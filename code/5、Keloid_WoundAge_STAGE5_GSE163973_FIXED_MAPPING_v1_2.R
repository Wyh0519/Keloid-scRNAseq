# =============================================================================
# 项目：瘢痕疙瘩“伤口状态持续与修复终止失败”
# 第五阶段：GSE163973独立瘢痕疙瘩队列QC、固定签名映射与患者级比较（V1.2修正版）
#
# 根目录：
#   C:/Users/33652/Desktop/LZZ/1、yssj
#
# 本阶段目的：
#   1. 读取GSE163973的6个单细胞矩阵（3例keloid，3例normal scar）；
#   2. 将作者all-cell metadata与原始矩阵barcode逐样本匹配；
#   3. 保守QC并按“患者×主细胞类型”生成pseudobulk；
#   4. 只使用第四阶段已固定、且在3次LOODO均稳定的伤口阶段基因；
#   5. 在共同基因背景上重新计算正常参考与疾病样本的秩评分；
#   6. 计算Skin/D1/D7/D30样状态、序数位置、早期持续和偏轨距离；
#   7. 以患者为统计单位比较keloid与normal scar；
#   8. 自动生成第五阶段检查包供复核。
#
# 统计边界：
#   - n=3 vs 3，不训练机器学习模型；
#   - 不使用细胞作为独立重复；
#   - 不在疾病队列中重新筛选签名基因；
#   - exact permutation P值只作描述，重点看效应量与方向一致性；
#   - 不将序数位置解释为精确伤口日龄。
#   - V1.2按样本内10x barcode核心序列匹配，完全去除整合对象追加的样本后缀。
#   - V1.2增加条码调试表，并修复后续PCA投影列名兼容问题。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANG = "en")

ROOT_DIR <- "C:/Users/33652/Desktop/LZZ/1、yssj"

RAW_TAR <- file.path(ROOT_DIR, "GSE163973_RAW.tar")
META_FILE <- file.path(
  ROOT_DIR,
  "GSE163973_integrate.all.NS.all.KL_cell.meta.data.csv.gz"
)

REFERENCE_PB_RDS <- file.path(
  ROOT_DIR,
  "03_STAGE3_REFERENCE_QC",
  "objects",
  "GSE241132_mainCellType_pseudobulk_counts.rds"
)

REFERENCE_MODEL_RDS <- file.path(
  ROOT_DIR,
  "04_STAGE4_WOUND_STATE_MODEL",
  "objects",
  "GSE241132_fixed_wound_state_reference.rds"
)

STAGE_DIR <- file.path(ROOT_DIR, "05_STAGE5_GSE163973_MAPPING")
REPORT_DIR <- file.path(STAGE_DIR, "report")
FIG_DIR <- file.path(STAGE_DIR, "figures")
OBJECT_DIR <- file.path(STAGE_DIR, "objects")
WORK_DIR <- file.path(STAGE_DIR, "_temporary_work")

PACKAGE_ZIP <- file.path(
  ROOT_DIR,
  "第五阶段_GSE163973固定伤口状态映射检查包.zip"
)
PACKAGE_TARGZ <- file.path(
  ROOT_DIR,
  "第五阶段_GSE163973固定伤口状态映射检查包.tar.gz"
)

if (dir.exists(STAGE_DIR)) {
  unlink(STAGE_DIR, recursive = TRUE, force = TRUE)
}

dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(REPORT_DIR, "00_stage5_log.txt")
WARNING_FILE <- file.path(REPORT_DIR, "00_stage5_warnings.txt")

log_msg <- function(...) {
  msg <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    paste(..., collapse = " ")
  )
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

warn_msg <- function(...) {
  msg <- paste(..., collapse = " ")
  log_msg("WARNING:", msg)
  cat(msg, "\n", file = WARNING_FILE, append = TRUE)
}

safe_write_csv <- function(x, path) {
  tryCatch(
    utils::write.csv(
      x,
      path,
      row.names = FALSE,
      fileEncoding = "UTF-8"
    ),
    error = function(e) {
      warn_msg(
        "CSV写入失败:",
        basename(path),
        "|",
        conditionMessage(e)
      )
    }
  )
}

safe_write_lines <- function(x, path) {
  tryCatch(
    writeLines(x, path, useBytes = TRUE),
    error = function(e) {
      warn_msg(
        "文本写入失败:",
        basename(path),
        "|",
        conditionMessage(e)
      )
    }
  )
}

normalize_slash <- function(x) {
  normalizePath(x, winslash = "/", mustWork = FALSE)
}

# ----------------------------- 依赖包 ------------------------------------------
if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages(
    "Matrix",
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("无法安装或加载Matrix包。")
}

# ----------------------------- 基础函数 ----------------------------------------
read_mtx_gz <- function(path) {
  con <- if (
    grepl("\\.gz$", path, ignore.case = TRUE)
  ) {
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }

  on.exit(close(con), add = TRUE)

  x <- Matrix::readMM(con)

  if (!inherits(x, "CsparseMatrix")) {
    x <- methods::as(x, "CsparseMatrix")
  }

  x
}

read_noheader_table <- function(
  path,
  sep = "\t"
) {
  con <- if (
    grepl("\\.gz$", path, ignore.case = TRUE)
  ) {
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }

  on.exit(close(con), add = TRUE)

  utils::read.table(
    con,
    header = FALSE,
    sep = sep,
    quote = "",
    comment.char = "",
    fill = TRUE,
    stringsAsFactors = FALSE
  )
}

read_metadata_csv <- function(path) {
  con <- if (
    grepl("\\.gz$", path, ignore.case = TRUE)
  ) {
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }

  on.exit(close(con), add = TRUE)

  d <- utils::read.csv(
    con,
    header = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # GEO文件首列是细胞ID，表头为空；R可能将其读入首列，也可能自动作为rownames。
  nm <- names(d)

  if (
    length(nm) &&
      (
        is.na(nm[1]) ||
          !nzchar(nm[1]) ||
          nm[1] %in% c(
            "X",
            "...1",
            "row.names"
          )
      )
  ) {
    names(d)[1] <- "cell_id_candidate"
  }

  d
}

normalize_10x_barcode <- function(x) {
  y <- toupper(
    trimws(
      as.character(x)
    )
  )

  out <- rep(
    NA_character_,
    length(y)
  )

  # GSE163973整合metadata中的_1、_2……是合并对象追加的样本序号，
  # 而每个独立10x矩阵中的条码后缀通常重新为-1。
  # 本脚本已经先按KL1-3、NS1-3分样本，因此只使用ACGT条码核心匹配，
  # 完全删除-1、_4、.6等数字后缀。
  hit <- grepl(
    "([ACGTN]{12,})(?:[-_.][0-9]+)?$",
    y,
    perl = TRUE
  )

  if (any(hit)) {
    out[hit] <- sub(
      "^.*?([ACGTN]{12,})(?:[-_.][0-9]+)?$",
      "\\1",
      y[hit],
      perl = TRUE
    )
  }

  out
}

extract_10x_barcode <- normalize_10x_barcode

select_metadata_barcode <- function(md) {
  candidates <- list()

  rn <- rownames(md)

  if (
    !is.null(rn) &&
      length(rn) == nrow(md) &&
      !identical(
        rn,
        as.character(seq_len(nrow(md)))
      )
  ) {
    candidates[[".rownames"]] <- rn
  }

  for (nm in names(md)) {
    v <- md[[nm]]

    if (
      is.character(v) ||
        is.factor(v)
    ) {
      candidates[[nm]] <- as.character(v)
    }
  }

  if (!length(candidates)) {
    return(
      list(
        name = NA_character_,
        values = rep(
          NA_character_,
          nrow(md)
        ),
        rate = 0
      )
    )
  }

  rates <- vapply(
    candidates,
    function(v) {
      mean(
        !is.na(
          extract_10x_barcode(v)
        )
      )
    },
    numeric(1)
  )

  best <- which.max(rates)

  list(
    name = names(candidates)[best],
    values = candidates[[best]],
    rate = rates[best]
  )
}

normalize_sample_id <- function(x) {
  raw <- toupper(
    trimws(
      as.character(x)
    )
  )

  out <- rep(
    NA_character_,
    length(raw)
  )

  # Keloid标签兼容：
  # KF1 / KL1 / K1 / keloid_1 / human_keloid_sample1
  k_hit <- grepl(
    "KELOID|(^|[^A-Z])KF[^A-Z0-9]*[0-9]+|(^|[^A-Z])KL[^A-Z0-9]*[0-9]+|^K[^A-Z0-9]*[0-9]+$",
    raw,
    perl = TRUE
  )

  if (any(k_hit)) {
    k_num <- sub(
      "^.*?(?:KELOID|KF|KL|K)[^0-9]*([0-9]+).*$",
      "\\1",
      raw[k_hit],
      perl = TRUE
    )

    valid_num <- grepl(
      "^[0-9]+$",
      k_num
    )

    idx <- which(k_hit)[valid_num]

    out[idx] <- paste0(
      "KL",
      k_num[valid_num]
    )
  }

  # Normal scar标签兼容：
  # NF1 / NS1 / N1 / scar_1 / normal_scar_1 / human_normal_scar_sample1
  n_hit <- is.na(out) & grepl(
    "NORMAL[^A-Z0-9]*SCAR|NORMALSCAR|(^|[^A-Z])NF[^A-Z0-9]*[0-9]+|(^|[^A-Z])NS[^A-Z0-9]*[0-9]+|(^|[^A-Z])SCAR[^A-Z0-9]*[0-9]+|^N[^A-Z0-9]*[0-9]+$",
    raw,
    perl = TRUE
  )

  if (any(n_hit)) {
    n_num <- sub(
      "^.*?(?:NORMAL[^A-Z0-9]*SCAR|NORMALSCAR|NF|NS|SCAR|N)[^0-9]*([0-9]+).*$",
      "\\1",
      raw[n_hit],
      perl = TRUE
    )

    valid_num <- grepl(
      "^[0-9]+$",
      n_num
    )

    idx <- which(n_hit)[valid_num]

    out[idx] <- paste0(
      "NS",
      n_num[valid_num]
    )
  }

  # 最后处理已经是标准形式但可能夹有分隔符的标签。
  compact <- gsub(
    "[^A-Z0-9]",
    "",
    raw
  )

  std_k <- is.na(out) & grepl(
    "^KL[0-9]+$",
    compact
  )

  out[std_k] <- compact[std_k]

  std_n <- is.na(out) & grepl(
    "^NS[0-9]+$",
    compact
  )

  out[std_n] <- compact[std_n]

  out
}

select_metadata_sample_column <- function(
  md,
  raw_sample_ids
) {
  preferred <- c(
    "orig.ident",
    "dataset",
    "sample",
    "sample_id",
    "Sample",
    "SampleID",
    "condition",
    "group"
  )

  character_columns <- names(md)[
    vapply(
      md,
      function(v) {
        is.character(v) ||
          is.factor(v)
      },
      logical(1)
    )
  ]

  candidate_names <- unique(
    c(
      intersect(
        preferred,
        names(md)
      ),
      character_columns
    )
  )

  if (!length(candidate_names)) {
    stop(
      "元数据中没有可用于恢复样本来源的字符字段。"
    )
  }

  normalized_list <- lapply(
    candidate_names,
    function(nm) {
      normalize_sample_id(
        md[[nm]]
      )
    }
  )

  names(normalized_list) <- candidate_names

  rates <- vapply(
    normalized_list,
    function(z) {
      mean(
        z %in% raw_sample_ids
      )
    },
    numeric(1)
  )

  best <- which.max(
    rates
  )

  selected <- normalized_list[[best]]

  # 使用其他字段逐行补充未映射样本，避免orig.ident与dataset命名不一致。
  for (nm in candidate_names) {
    z <- normalized_list[[nm]]

    fill <- (
      is.na(selected) &
        z %in% raw_sample_ids
    )

    selected[fill] <- z[fill]
  }

  list(
    name = candidate_names[best],
    normalized = selected,
    rate = mean(
      selected %in% raw_sample_ids
    ),
    field_rates = data.frame(
      field = candidate_names,
      mapping_rate = as.numeric(
        rates
      ),
      stringsAsFactors = FALSE
    )
  )
}

build_metadata_sample_value_map <- function(
  md,
  selected_field,
  normalized_values
) {
  original <- trimws(
    as.character(
      md[[selected_field]]
    )
  )

  x <- data.frame(
    original_value = original,
    normalized_sample_id = normalized_values,
    stringsAsFactors = FALSE
  )

  x <- unique(
    x
  )

  x <- x[
    order(
      x$normalized_sample_id,
      x$original_value,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  rownames(x) <- NULL
  x
}


harmonize_cell_type <- function(x) {
  y <- tolower(
    trimws(
      as.character(x)
    )
  )

  out <- rep(
    NA_character_,
    length(y)
  )

  out[
    grepl("fibro", y)
  ] <- "Fibroblast"

  out[
    is.na(out) &
      grepl("keratino|epithelial", y)
  ] <- "Keratinocyte"

  out[
    is.na(out) &
      grepl("endothelial", y)
  ] <- "Endothelial"

  out[
    is.na(out) &
      grepl(
        "myeloid|macroph|monocyte|dendritic|(^|[_ -])dc([_ -]|$)",
        y
      )
  ] <- "Myeloid"

  out[
    is.na(out) &
      grepl(
        "lymph|t[_ -]?cell|b[_ -]?cell|nk|natural[_ -]?killer",
        y
      )
  ] <- "Lymphoid"

  out
}

sample_id_from_nested_tar <- function(path) {
  x <- basename(path)

  x <- sub(
    "\\.tar\\.gz$",
    "",
    x,
    ignore.case = TRUE
  )

  x <- sub(
    "^GSM[0-9]+_",
    "",
    x
  )

  x <- sub(
    "_matrix$",
    "",
    x,
    ignore.case = TRUE
  )

  toupper(x)
}

extract_nested_triplet <- function(
  nested_tar,
  sample_tmp
) {
  if (dir.exists(sample_tmp)) {
    unlink(
      sample_tmp,
      recursive = TRUE,
      force = TRUE
    )
  }

  dir.create(
    sample_tmp,
    recursive = TRUE,
    showWarnings = FALSE
  )

  utils::untar(
    nested_tar,
    exdir = sample_tmp
  )

  files <- list.files(
    sample_tmp,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  )

  locate <- function(regex) {
    hit <- files[
      grepl(
        regex,
        basename(files),
        ignore.case = TRUE,
        perl = TRUE
      )
    ]

    if (!length(hit)) {
      return(NA_character_)
    }

    hit[1]
  }

  out <- list(
    matrix = locate(
      "matrix\\.mtx(\\.gz)?$"
    ),
    barcodes = locate(
      "barcodes?\\.tsv(\\.gz)?$"
    ),
    features = locate(
      "(features|genes)\\.tsv(\\.gz)?$"
    )
  )

  if (any(vapply(out, is.na, logical(1)))) {
    stop(
      "嵌套压缩包缺少matrix/barcodes/features：",
      basename(nested_tar)
    )
  }

  out
}

first_unique_nonempty <- function(
  x,
  label,
  sample_id
) {
  x <- trimws(
    as.character(x)
  )

  x <- unique(
    x[
      !is.na(x) &
        nzchar(x)
    ]
  )

  if (!length(x)) {
    warn_msg(
      sample_id,
      label,
      "为空，写入NA"
    )

    return(NA_character_)
  }

  if (length(x) > 1L) {
    warn_msg(
      sample_id,
      label,
      "存在多个取值：",
      paste(
        x,
        collapse = " | "
      ),
      "；使用第一个值"
    )
  }

  x[1]
}

collapse_duplicate_symbols <- function(
  counts,
  gene_symbol
) {
  gene_symbol <- trimws(
    as.character(gene_symbol)
  )

  valid <- (
    !is.na(gene_symbol) &
      nzchar(gene_symbol)
  )

  counts <- counts[
    valid,
    ,
    drop = FALSE
  ]

  gene_symbol <- gene_symbol[valid]

  if (!is.matrix(counts)) {
    counts <- as.matrix(counts)
  }

  collapsed <- rowsum(
    counts,
    group = gene_symbol,
    reorder = FALSE
  )

  storage.mode(collapsed) <- "double"

  collapsed
}

log_cpm <- function(counts) {
  lib <- colSums(counts)

  if (any(lib <= 0)) {
    stop("发现文库总计数≤0。")
  }

  log2(
    sweep(
      counts + 0.5,
      2,
      lib + 1,
      "/"
    ) * 1e6
  )
}

exclude_background_gene <- function(gene) {
  g <- toupper(
    as.character(gene)
  )

  grepl(
    "^MT-|^RPS[0-9]|^RPL[0-9]|^HB[ABDEGMQZ][0-9A-Z]*$",
    g
  )
}

score_samples_by_rank <- function(
  expr,
  signatures
) {
  if (!all(
    c(
      "Skin",
      "Wound1",
      "Wound7",
      "Wound30"
    ) %in% names(signatures)
  )) {
    stop("签名缺少四个阶段。")
  }

  out <- matrix(
    NA_real_,
    nrow = ncol(expr),
    ncol = 4,
    dimnames = list(
      colnames(expr),
      c(
        "Skin",
        "Wound1",
        "Wound7",
        "Wound30"
      )
    )
  )

  for (j in seq_len(ncol(expr))) {
    ranks <- rank(
      expr[, j],
      ties.method = "average",
      na.last = "keep"
    )

    ranks <- ranks /
      sum(is.finite(ranks))

    names(ranks) <- rownames(expr)

    for (st in colnames(out)) {
      genes <- intersect(
        signatures[[st]],
        names(ranks)
      )

      if (length(genes) < 10L) {
        stop(
          "阶段",
          st,
          "在共同基因背景中不足10个，实际为：",
          length(genes)
        )
      }

      out[j, st] <- mean(
        ranks[genes],
        na.rm = TRUE
      )
    }
  }

  out
}

get_reference_design <- function(
  counts,
  sample_info
) {
  sample_info$key <- paste(
    sample_info$sample_id,
    sample_info$main_cell_type,
    sep = "||"
  )

  idx <- match(
    colnames(counts),
    sample_info$key
  )

  if (anyNA(idx)) {
    stop(
      "参考pseudobulk列名无法全部匹配sample_info。"
    )
  }

  out <- sample_info[
    idx,
    ,
    drop = FALSE
  ]

  if (!identical(
    colnames(counts),
    out$key
  )) {
    stop(
      "参考pseudobulk设计顺序匹配失败。"
    )
  }

  out
}

get_stable_signatures <- function(
  model,
  cell_type
) {
  st <- model$signature_stability

  required_cols <- c(
    "cell_type",
    "stage",
    "gene",
    "stable_all_3_folds"
  )

  if (!all(
    required_cols %in% names(st)
  )) {
    stop(
      "固定模型中的signature_stability结构不完整。"
    )
  }

  x <- st[
    st$cell_type == cell_type &
      st$stable_all_3_folds,
    ,
    drop = FALSE
  ]

  stages <- c(
    "Skin",
    "Wound1",
    "Wound7",
    "Wound30"
  )

  sig <- setNames(
    lapply(
      stages,
      function(s) {
        unique(
          as.character(
            x$gene[
              x$stage == s
            ]
          )
        )
      }
    ),
    stages
  )

  n_genes <- vapply(
    sig,
    length,
    integer(1)
  )

  if (any(n_genes < 20L)) {
    stop(
      cell_type,
      "存在稳定签名少于20个基因：",
      paste(
        names(n_genes),
        n_genes,
        sep = "=",
        collapse = ", "
      )
    )
  }

  sig
}

standardize_scores <- function(
  ref_scores,
  new_scores
) {
  mu <- colMeans(
    ref_scores
  )

  sig <- apply(
    ref_scores,
    2,
    stats::sd
  )

  sig[
    !is.finite(sig) |
      sig == 0
  ] <- 1

  ref_z <- sweep(
    sweep(
      ref_scores,
      2,
      mu,
      "-"
    ),
    2,
    sig,
    "/"
  )

  new_z <- sweep(
    sweep(
      new_scores,
      2,
      mu,
      "-"
    ),
    2,
    sig,
    "/"
  )

  list(
    ref_z = ref_z,
    new_z = new_z,
    center = mu,
    scale = sig
  )
}

calculate_reference_centroids <- function(
  ref_z,
  ref_design
) {
  stages <- c(
    "Skin",
    "Wound1",
    "Wound7",
    "Wound30"
  )

  centroids <- do.call(
    rbind,
    lapply(
      stages,
      function(st) {
        colMeans(
          ref_z[
            ref_design$condition == st,
            ,
            drop = FALSE
          ]
        )
      }
    )
  )

  rownames(centroids) <- stages

  centroids
}

calculate_lodo_reference_distance <- function(
  ref_z,
  ref_design
) {
  out <- numeric(
    nrow(ref_z)
  )

  for (i in seq_len(nrow(ref_z))) {
    st <- ref_design$condition[i]
    donor <- ref_design$donor[i]

    train_idx <- which(
      ref_design$condition == st &
        ref_design$donor != donor
    )

    centroid <- colMeans(
      ref_z[
        train_idx,
        ,
        drop = FALSE
      ]
    )

    out[i] <- sqrt(
      sum(
        (
          ref_z[i, ] -
            centroid
        ) ^ 2
      )
    )
  }

  out
}

map_to_reference <- function(
  z_scores,
  centroids,
  reference_distance_limit
) {
  stage_position <- c(
    Skin = 0,
    Wound1 = 1,
    Wound7 = 2,
    Wound30 = 3
  )

  rows <- vector(
    "list",
    nrow(z_scores)
  )

  for (i in seq_len(nrow(z_scores))) {
    d <- apply(
      centroids,
      1,
      function(cen) {
        sqrt(
          sum(
            (
              z_scores[i, ] -
                cen
            ) ^ 2
          )
        )
      }
    )

    nearest <- names(
      which.min(d)
    )

    w <- exp(
      -0.5 * d ^ 2
    )

    if (
      !all(is.finite(w)) ||
        sum(w) == 0
    ) {
      w <- rep(
        1 / length(d),
        length(d)
      )
    } else {
      w <- w / sum(w)
    }

    entropy <- -sum(
      w * log(
        pmax(w, 1e-12)
      )
    ) / log(length(w))

    rows[[i]] <- data.frame(
      sample_key = rownames(z_scores)[i],
      nearest_stage = nearest,
      ordinal_position = sum(
        w *
          stage_position[
            names(w)
          ]
      ),
      off_trajectory_distance = min(d),
      off_trajectory_ratio =
        min(d) /
        reference_distance_limit,
      mixed_state_entropy = entropy,
      weight_Skin = w["Skin"],
      weight_Wound1 = w["Wound1"],
      weight_Wound7 = w["Wound7"],
      weight_Wound30 = w["Wound30"],
      stringsAsFactors = FALSE
    )
  }

  do.call(
    rbind,
    rows
  )
}

add_score_contrasts <- function(
  mapping,
  z_scores
) {
  idx <- match(
    mapping$sample_key,
    rownames(z_scores)
  )

  if (anyNA(idx)) {
    stop(
      "状态映射结果与z-score无法匹配。"
    )
  }

  z <- z_scores[
    idx,
    ,
    drop = FALSE
  ]

  mapping$z_Skin <- z[, "Skin"]
  mapping$z_Wound1 <- z[, "Wound1"]
  mapping$z_Wound7 <- z[, "Wound7"]
  mapping$z_Wound30 <- z[, "Wound30"]

  mapping$early_state_persistence <-
    (
      z[, "Wound1"] +
        z[, "Wound7"]
    ) / 2 -
    (
      z[, "Skin"] +
        z[, "Wound30"]
    ) / 2

  mapping$remodeling_completion <-
    z[, "Wound30"] -
    (
      z[, "Wound1"] +
        z[, "Wound7"]
    ) / 2

  mapping$wound_activation <-
    (
      z[, "Wound1"] +
        z[, "Wound7"] +
        z[, "Wound30"]
    ) / 3 -
    z[, "Skin"]

  mapping
}

exact_permutation_test <- function(
  value,
  group,
  group_a = "keloid",
  group_b = "normal_scar"
) {
  ok <- (
    is.finite(value) &
      group %in% c(
        group_a,
        group_b
      )
  )

  value <- value[ok]
  group <- group[ok]

  x <- value[
    group == group_a
  ]

  y <- value[
    group == group_b
  ]

  n1 <- length(x)
  n2 <- length(y)

  if (
    n1 < 2L ||
      n2 < 2L
  ) {
    return(
      list(
        mean_a = NA_real_,
        mean_b = NA_real_,
        median_a = NA_real_,
        median_b = NA_real_,
        mean_difference = NA_real_,
        hodges_lehmann_shift = NA_real_,
        hedges_g = NA_real_,
        cliffs_delta = NA_real_,
        exact_permutation_p = NA_real_
      )
    )
  }

  obs <- mean(x) - mean(y)

  all_values <- c(
    x,
    y
  )

  combos <- utils::combn(
    seq_along(all_values),
    n1
  )

  perm_diff <- apply(
    combos,
    2,
    function(id_a) {
      mean(
        all_values[id_a]
      ) -
        mean(
          all_values[-id_a]
        )
    }
  )

  p_exact <- mean(
    abs(perm_diff) >=
      abs(obs) - 1e-12
  )

  df <- n1 + n2 - 2

  pooled_sd <- sqrt(
    (
      (n1 - 1) * stats::var(x) +
        (n2 - 1) * stats::var(y)
    ) / df
  )

  hedges_g <- NA_real_

  if (
    is.finite(pooled_sd) &&
      pooled_sd > 0
  ) {
    cohen_d <- obs / pooled_sd
    correction <- 1 - 3 / (4 * df - 1)
    hedges_g <- correction * cohen_d
  }

  cliffs_delta <- mean(
    outer(
      x,
      y,
      "-"
    ) > 0
  ) -
    mean(
      outer(
        x,
        y,
        "-"
      ) < 0
    )

  list(
    mean_a = mean(x),
    mean_b = mean(y),
    median_a = stats::median(x),
    median_b = stats::median(y),
    mean_difference = obs,
    hodges_lehmann_shift = stats::median(
      as.vector(
        outer(
          x,
          y,
          "-"
        )
      )
    ),
    hedges_g = hedges_g,
    cliffs_delta = cliffs_delta,
    exact_permutation_p = p_exact
  )
}

effect_strength <- function(
  hedges_g,
  cliffs_delta
) {
  if (
    !is.finite(hedges_g) &&
      !is.finite(cliffs_delta)
  ) {
    return("not_estimable")
  }

  if (
    (
      is.finite(hedges_g) &&
        abs(hedges_g) >= 0.8
    ) ||
      (
        is.finite(cliffs_delta) &&
          abs(cliffs_delta) >= 0.56
      )
  ) {
    return("large")
  }

  if (
    (
      is.finite(hedges_g) &&
        abs(hedges_g) >= 0.5
    ) ||
      (
        is.finite(cliffs_delta) &&
          abs(cliffs_delta) >= 0.33
      )
  ) {
    return("moderate")
  }

  "small"
}

# ----------------------------- 输入检查 ----------------------------------------
log_msg("第五阶段开始")

required_files <- c(
  RAW_TAR,
  META_FILE,
  REFERENCE_PB_RDS,
  REFERENCE_MODEL_RDS
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files)) {
  stop(
    "缺少必要文件：",
    paste(
      missing_files,
      collapse = " | "
    )
  )
}

# ----------------------------- 解包与元数据 ------------------------------------
outer_dir <- file.path(
  WORK_DIR,
  "outer"
)

dir.create(
  outer_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

log_msg(
  "解包GSE163973_RAW.tar"
)

utils::untar(
  RAW_TAR,
  exdir = outer_dir
)

nested_tars <- list.files(
  outer_dir,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.tar\\.gz$",
  ignore.case = TRUE
)

raw_sample_ids <- vapply(
  nested_tars,
  sample_id_from_nested_tar,
  character(1)
)

expected_ids <- c(
  "KL1",
  "KL2",
  "KL3",
  "NS1",
  "NS2",
  "NS3"
)

if (
  length(nested_tars) != 6L ||
    !setequal(
      raw_sample_ids,
      expected_ids
    )
) {
  stop(
    "GSE163973样本不是预期的KL1-3与NS1-3。实际：",
    paste(
      raw_sample_ids,
      collapse = ", "
    )
  )
}

nested_tars <- nested_tars[
  match(
    expected_ids,
    raw_sample_ids
  )
]

raw_sample_ids <- expected_ids

log_msg(
  "读取all-cell metadata"
)

metadata <- read_metadata_csv(
  META_FILE
)

required_meta_cols <- c(
  "orig.ident",
  "dataset",
  "cellType"
)

missing_meta_cols <- setdiff(
  required_meta_cols,
  names(metadata)
)

if (length(missing_meta_cols)) {
  stop(
    "元数据缺少字段：",
    paste(
      missing_meta_cols,
      collapse = ", "
    )
  )
}

bc_info <- select_metadata_barcode(
  metadata
)

metadata$barcode_original <- as.character(
  bc_info$values
)

metadata$barcode_norm <- extract_10x_barcode(
  metadata$barcode_original
)

sample_info_detected <- select_metadata_sample_column(
  metadata,
  raw_sample_ids
)

metadata$sample_id_norm <-
  sample_info_detected$normalized

safe_write_csv(
  sample_info_detected$field_rates,
  file.path(
    REPORT_DIR,
    "01A_metadata_sample_field_mapping_rates.csv"
  )
)

sample_value_map <- build_metadata_sample_value_map(
  metadata,
  sample_info_detected$name,
  metadata$sample_id_norm
)

safe_write_csv(
  sample_value_map,
  file.path(
    REPORT_DIR,
    "01B_metadata_sample_value_mapping.csv"
  )
)

mapped_sample_counts <- as.data.frame(
  table(
    metadata$sample_id_norm,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

names(mapped_sample_counts) <- c(
  "sample_id_norm",
  "n_metadata_cells"
)

safe_write_csv(
  mapped_sample_counts,
  file.path(
    REPORT_DIR,
    "01C_metadata_cells_by_normalized_sample.csv"
  )
)

missing_expected_metadata <- setdiff(
  raw_sample_ids,
  unique(
    metadata$sample_id_norm[
      !is.na(
        metadata$sample_id_norm
      )
    ]
  )
)

if (length(missing_expected_metadata)) {
  stop(
    "元数据样本映射仍缺少：",
    paste(
      missing_expected_metadata,
      collapse = ", "
    ),
    "。请查看01A、01B、01C审计表。"
  )
}

metadata$cell_type_original <-
  trimws(
    as.character(
      metadata$cellType
    )
  )

metadata$main_cell_type <-
  harmonize_cell_type(
    metadata$cell_type_original
  )

if (bc_info$rate < 0.95) {
  warn_msg(
    "barcode字段自动识别率低于95%：",
    bc_info$name,
    "=",
    round(
      bc_info$rate,
      4
    )
  )
}

if (sample_info_detected$rate < 0.95) {
  warn_msg(
    "样本字段映射率低于95%：",
    sample_info_detected$name,
    "=",
    round(
      sample_info_detected$rate,
      4
    )
  )
}

metadata_audit <- data.frame(
  metadata_rows = nrow(metadata),
  barcode_source = bc_info$name,
  barcode_detection_rate = bc_info$rate,
  sample_source = sample_info_detected$name,
  sample_mapping_rate = sample_info_detected$rate,
  n_unmapped_barcodes = sum(
    is.na(
      metadata$barcode_norm
    )
  ),
  n_unmapped_samples = sum(
    is.na(
      metadata$sample_id_norm
    )
  ),
  n_unmapped_main_celltypes = sum(
    is.na(
      metadata$main_cell_type
    )
  ),
  all_six_expected_samples_present =
    length(
      setdiff(
        raw_sample_ids,
        unique(
          metadata$sample_id_norm[
            !is.na(
              metadata$sample_id_norm
            )
          ]
        )
      )
    ) == 0L,
  stringsAsFactors = FALSE
)

safe_write_csv(
  metadata_audit,
  file.path(
    REPORT_DIR,
    "01_metadata_mapping_audit.csv"
  )
)

original_celltype_counts <- as.data.frame(
  table(
    metadata$sample_id_norm,
    metadata$cell_type_original,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

names(original_celltype_counts) <- c(
  "sample_id",
  "cell_type_original",
  "n_cells"
)

original_celltype_counts <-
  original_celltype_counts[
    original_celltype_counts$n_cells > 0,
    ,
    drop = FALSE
  ]

safe_write_csv(
  original_celltype_counts,
  file.path(
    REPORT_DIR,
    "02_author_celltype_counts.csv"
  )
)

# ----------------------------- 逐样本QC与pseudobulk ----------------------------
qc_rows <- list()
barcode_rows <- list()
celltype_rows <- list()
pseudobulk_list <- list()
pb_info_rows <- list()
feature_rows <- list()
barcode_debug_rows <- list()

gene_reference_id <- NULL
gene_reference_symbol <- NULL

pdf(
  file.path(
    FIG_DIR,
    "Figure_GSE163973_QC_distributions.pdf"
  ),
  width = 10,
  height = 8,
  onefile = TRUE
)

for (nested_tar in nested_tars) {
  sid <- sample_id_from_nested_tar(
    nested_tar
  )

  log_msg(
    "处理样本:",
    sid
  )

  sample_tmp <- file.path(
    WORK_DIR,
    paste0(
      "sample_",
      sid
    )
  )

  triplet <- extract_nested_triplet(
    nested_tar,
    sample_tmp
  )

  barcodes_df <- read_noheader_table(
    triplet$barcodes,
    sep = "\t"
  )

  features_df <- read_noheader_table(
    triplet$features,
    sep = "\t"
  )

  counts <- read_mtx_gz(
    triplet$matrix
  )

  raw_barcode_original <- as.character(
    barcodes_df[[1]]
  )

  raw_barcodes <- normalize_10x_barcode(
    raw_barcode_original
  )

  if (
    anyNA(raw_barcodes) ||
      anyDuplicated(raw_barcodes)
  ) {
    stop(
      sid,
      "的矩阵barcode标准化失败或出现重复：",
      " NA=",
      sum(
        is.na(
          raw_barcodes
        )
      ),
      "; duplicated=",
      anyDuplicated(
        raw_barcodes
      )
    )
  }

  gene_id <- as.character(
    features_df[[1]]
  )

  gene_symbol <- if (
    ncol(features_df) >= 2L
  ) {
    as.character(
      features_df[[2]]
    )
  } else {
    gene_id
  }

  if (
    nrow(counts) != length(gene_id) ||
      ncol(counts) != length(raw_barcodes)
  ) {
    stop(
      "矩阵维度与features/barcodes不匹配：",
      sid
    )
  }

  if (is.null(gene_reference_id)) {
    gene_reference_id <- gene_id
    gene_reference_symbol <- gene_symbol
    identical_id <- TRUE
    identical_symbol <- TRUE
  } else {
    identical_id <- identical(
      gene_id,
      gene_reference_id
    )

    identical_symbol <- identical(
      gene_symbol,
      gene_reference_symbol
    )

    if (!identical_id) {
      if (
        anyDuplicated(gene_id) ||
          anyDuplicated(gene_reference_id)
      ) {
        stop(
          sid,
          "gene_id顺序不同且存在重复，无法安全对齐。"
        )
      }

      reorder_idx <- match(
        gene_reference_id,
        gene_id
      )

      if (anyNA(reorder_idx)) {
        stop(
          sid,
          "缺少参考gene_id，无法安全对齐。"
        )
      }

      counts <- counts[
        reorder_idx,
        ,
        drop = FALSE
      ]

      gene_id <- gene_id[
        reorder_idx
      ]

      gene_symbol <- gene_symbol[
        reorder_idx
      ]

      identical_id <- identical(
        gene_id,
        gene_reference_id
      )

      identical_symbol <- identical(
        gene_symbol,
        gene_reference_symbol
      )
    }
  }

  rownames(counts) <- make.unique(
    gene_reference_symbol
  )

  colnames(counts) <- raw_barcodes

  feature_rows[[length(feature_rows) + 1L]] <-
    data.frame(
      sample_id = sid,
      n_features = length(gene_id),
      identical_gene_ids_to_first = identical_id,
      identical_gene_symbols_to_first = identical_symbol,
      stringsAsFactors = FALSE
    )

  md_s <- metadata[
    metadata$sample_id_norm == sid,
    ,
    drop = FALSE
  ]

  md_s <- md_s[
    !is.na(
      md_s$barcode_norm
    ),
    ,
    drop = FALSE
  ]

  md_s <- md_s[
    !duplicated(
      md_s$barcode_norm
    ),
    ,
    drop = FALSE
  ]

  if (!nrow(md_s)) {
    available_ids <- sort(
      unique(
        metadata$sample_id_norm[
          !is.na(
            metadata$sample_id_norm
          )
        ]
      )
    )

    stop(
      "元数据中没有映射到样本：",
      sid,
      "；当前成功映射的样本为：",
      paste(
        available_ids,
        collapse = ", "
      ),
      "。请查看01A、01B、01C审计表。"
    )
  }

  n_debug_meta <- min(
    20L,
    nrow(md_s)
  )

  n_debug_raw <- min(
    20L,
    length(raw_barcodes)
  )

  barcode_debug_rows[[length(barcode_debug_rows) + 1L]] <-
    rbind(
      data.frame(
        sample_id = sid,
        source = "metadata",
        original_barcode = md_s$barcode_original[
          seq_len(n_debug_meta)
        ],
        normalized_core = md_s$barcode_norm[
          seq_len(n_debug_meta)
        ],
        stringsAsFactors = FALSE
      ),
      data.frame(
        sample_id = sid,
        source = "raw_matrix",
        original_barcode = raw_barcode_original[
          seq_len(n_debug_raw)
        ],
        normalized_core = raw_barcodes[
          seq_len(n_debug_raw)
        ],
        stringsAsFactors = FALSE
      )
    )

  idx <- match(
    md_s$barcode_norm,
    raw_barcodes
  )

  matched <- !is.na(idx)

  match_rate <- mean(
    matched
  )

  barcode_rows[[length(barcode_rows) + 1L]] <-
    data.frame(
      sample_id = sid,
      n_matrix_barcodes = length(
        raw_barcodes
      ),
      n_metadata_rows = nrow(md_s),
      n_metadata_matched = sum(
        matched
      ),
      metadata_match_rate = match_rate,
      n_matrix_not_in_metadata =
        length(raw_barcodes) -
        length(
          unique(
            idx[matched]
          )
        ),
      stringsAsFactors = FALSE
    )

  if (match_rate < 0.90) {
    mismatch_audit <- data.frame(
      sample_id = sid,
      metadata_original = head(
        md_s$barcode_original,
        100L
      ),
      metadata_core = head(
        md_s$barcode_norm,
        100L
      ),
      metadata_core_found_in_raw = head(
        md_s$barcode_norm %in% raw_barcodes,
        100L
      ),
      stringsAsFactors = FALSE
    )

    safe_write_csv(
      mismatch_audit,
      file.path(
        REPORT_DIR,
        paste0(
          "DEBUG_barcode_mismatch_",
          sid,
          ".csv"
        )
      )
    )

    raw_audit <- data.frame(
      sample_id = sid,
      raw_original = head(
        raw_barcode_original,
        100L
      ),
      raw_core = head(
        raw_barcodes,
        100L
      ),
      stringsAsFactors = FALSE
    )

    safe_write_csv(
      raw_audit,
      file.path(
        REPORT_DIR,
        paste0(
          "DEBUG_raw_barcodes_",
          sid,
          ".csv"
        )
      )
    )

    warn_msg(
      sid,
      "metadata与矩阵barcode匹配率低于90%：",
      round(
        match_rate,
        4
      ),
      "；已输出DEBUG条码审计表。"
    )
  }

  if (!any(matched)) {
    stop(
      sid,
      "没有任何metadata barcode匹配表达矩阵。"
    )
  }

  md_s <- md_s[
    matched,
    ,
    drop = FALSE
  ]

  counts <- counts[
    ,
    idx[matched],
    drop = FALSE
  ]

  colnames(counts) <-
    md_s$barcode_norm

  n_count <- Matrix::colSums(
    counts
  )

  n_feature <- Matrix::colSums(
    counts > 0
  )

  gene_upper <- toupper(
    rownames(counts)
  )

  mt_idx <- grepl(
    "^MT-",
    gene_upper
  )

  ribo_idx <- grepl(
    "^RP[SL][0-9]",
    gene_upper
  )

  hb_idx <- grepl(
    "^HB[ABDEGMQZ][0-9A-Z]*$",
    gene_upper
  )

  pct_mt <- if (any(mt_idx)) {
    Matrix::colSums(
      counts[
        mt_idx,
        ,
        drop = FALSE
      ]
    ) /
      pmax(
        n_count,
        1
      ) * 100
  } else {
    rep(
      0,
      ncol(counts)
    )
  }

  pct_ribo <- if (any(ribo_idx)) {
    Matrix::colSums(
      counts[
        ribo_idx,
        ,
        drop = FALSE
      ]
    ) /
      pmax(
        n_count,
        1
      ) * 100
  } else {
    rep(
      0,
      ncol(counts)
    )
  }

  pct_hb <- if (any(hb_idx)) {
    Matrix::colSums(
      counts[
        hb_idx,
        ,
        drop = FALSE
      ]
    ) /
      pmax(
        n_count,
        1
      ) * 100
  } else {
    rep(
      0,
      ncol(counts)
    )
  }

  md_s$nCount_calc <- as.numeric(
    n_count
  )

  md_s$nFeature_calc <- as.numeric(
    n_feature
  )

  md_s$percent_mt_calc <- as.numeric(
    pct_mt
  )

  md_s$percent_ribo_calc <- as.numeric(
    pct_ribo
  )

  md_s$percent_hb_calc <- as.numeric(
    pct_hb
  )

  # 作者已经提供过滤后的细胞。这里只做保守二次QC，不新增双细胞算法。
  md_s$qc_keep <- (
    is.finite(
      md_s$nFeature_calc
    ) &
      is.finite(
        md_s$nCount_calc
      ) &
      is.finite(
        md_s$percent_mt_calc
      ) &
      md_s$nFeature_calc >= 200 &
      md_s$nFeature_calc <= 7500 &
      md_s$nCount_calc >= 500 &
      md_s$percent_mt_calc <= 20
  )

  md_s$qc_keep[
    is.na(
      md_s$qc_keep
    )
  ] <- FALSE

  group_s <- if (
    startsWith(
      sid,
      "KL"
    )
  ) {
    "keloid"
  } else {
    "normal_scar"
  }

  qc_rows[[length(qc_rows) + 1L]] <-
    data.frame(
      sample_id = sid,
      group = group_s,
      n_matrix_barcodes = length(
        raw_barcodes
      ),
      n_author_metadata_cells = nrow(
        md_s
      ),
      n_qc_keep = sum(
        md_s$qc_keep
      ),
      qc_retention_rate = mean(
        md_s$qc_keep
      ),
      median_nFeature = stats::median(
        md_s$nFeature_calc
      ),
      median_nCount = stats::median(
        md_s$nCount_calc
      ),
      median_percent_mt = stats::median(
        md_s$percent_mt_calc
      ),
      stringsAsFactors = FALSE
    )

  op <- par(
    no.readonly = TRUE
  )

  par(
    mfrow = c(2, 2),
    mar = c(4, 4, 3, 1)
  )

  hist(
    md_s$nFeature_calc,
    breaks = 50,
    main = paste0(
      sid,
      " nFeature"
    ),
    xlab = "Detected genes"
  )

  abline(
    v = c(
      200,
      7500
    ),
    lty = 2
  )

  hist(
    log10(
      md_s$nCount_calc + 1
    ),
    breaks = 50,
    main = paste0(
      sid,
      " log10 nCount"
    ),
    xlab = "log10 counts + 1"
  )

  abline(
    v = log10(
      501
    ),
    lty = 2
  )

  hist(
    md_s$percent_mt_calc,
    breaks = 50,
    main = paste0(
      sid,
      " mitochondrial %"
    ),
    xlab = "percent.mt"
  )

  abline(
    v = 20,
    lty = 2
  )

  plot(
    log10(
      md_s$nCount_calc + 1
    ),
    md_s$nFeature_calc,
    pch = 16,
    cex = 0.35,
    xlab = "log10 counts + 1",
    ylab = "Detected genes",
    main = paste0(
      sid,
      " QC scatter"
    )
  )

  par(op)

  keep_idx <- which(
    md_s$qc_keep
  )

  counts_qc <- counts[
    ,
    keep_idx,
    drop = FALSE
  ]

  md_qc <- md_s[
    keep_idx,
    ,
    drop = FALSE
  ]

  if (!nrow(md_qc)) {
    stop(
      sid,
      "在保守QC后没有保留任何细胞。"
    )
  }

  ct_input <- data.frame(
    sample_id = rep(
      sid,
      nrow(md_qc)
    ),
    group = rep(
      group_s,
      nrow(md_qc)
    ),
    main_cell_type = trimws(
      as.character(
        md_qc$main_cell_type
      )
    ),
    stringsAsFactors = FALSE
  )

  ct_input <- ct_input[
    !is.na(
      ct_input$main_cell_type
    ) &
      nzchar(
        ct_input$main_cell_type
      ),
    ,
    drop = FALSE
  ]

  if (nrow(ct_input)) {
    ct_df <- stats::aggregate(
      x = list(
        n_cells = rep(
          1L,
          nrow(ct_input)
        )
      ),
      by = ct_input,
      FUN = sum
    )

    celltype_rows[[length(celltype_rows) + 1L]] <-
      ct_df
  }

  cell_types <- sort(
    unique(
      ct_input$main_cell_type
    )
  )

  for (ct in cell_types) {
    j <- which(
      trimws(
        as.character(
          md_qc$main_cell_type
        )
      ) == ct
    )

    if (length(j) < 20L) {
      next
    }

    pb <- Matrix::rowSums(
      counts_qc[
        ,
        j,
        drop = FALSE
      ]
    )

    names(pb) <- make.unique(
      gene_reference_symbol
    )

    key <- paste(
      sid,
      ct,
      sep = "||"
    )

    pseudobulk_list[[key]] <- pb

    pb_info_rows[[length(pb_info_rows) + 1L]] <-
      data.frame(
        sample_key = key,
        sample_id = sid,
        group = group_s,
        main_cell_type = ct,
        n_cells = length(j),
        pseudobulk_library_size = sum(
          pb
        ),
        pseudobulk_detected_genes = sum(
          pb > 0
        ),
        stringsAsFactors = FALSE
      )
  }

  rm(
    counts,
    counts_qc,
    md_s,
    md_qc
  )

  gc(
    verbose = FALSE
  )

  unlink(
    sample_tmp,
    recursive = TRUE,
    force = TRUE
  )
}

dev.off()

if (length(barcode_debug_rows)) {
  barcode_debug_summary <- do.call(
    rbind,
    barcode_debug_rows
  )

  safe_write_csv(
    barcode_debug_summary,
    file.path(
      REPORT_DIR,
      "02A_barcode_normalization_examples.csv"
    )
  )
}

if (!length(qc_rows)) {
  stop("未生成任何样本QC结果。")
}

if (!length(barcode_rows)) {
  stop("未生成任何barcode匹配结果。")
}

if (!length(celltype_rows)) {
  stop("未生成任何细胞类型计数。")
}

if (!length(pb_info_rows)) {
  stop("未生成任何pseudobulk信息。")
}

if (!length(feature_rows)) {
  stop("未生成features一致性结果。")
}

qc_summary <- do.call(
  rbind,
  qc_rows
)

barcode_summary <- do.call(
  rbind,
  barcode_rows
)

celltype_counts <- do.call(
  rbind,
  celltype_rows
)

pb_info <- do.call(
  rbind,
  pb_info_rows
)

feature_consistency <- do.call(
  rbind,
  feature_rows
)

safe_write_csv(
  qc_summary,
  file.path(
    REPORT_DIR,
    "03_sample_QC_summary.csv"
  )
)

safe_write_csv(
  barcode_summary,
  file.path(
    REPORT_DIR,
    "04_barcode_matching_summary.csv"
  )
)

safe_write_csv(
  celltype_counts,
  file.path(
    REPORT_DIR,
    "05_harmonized_celltype_counts.csv"
  )
)

safe_write_csv(
  pb_info,
  file.path(
    REPORT_DIR,
    "06_disease_pseudobulk_information.csv"
  )
)

safe_write_csv(
  feature_consistency,
  file.path(
    REPORT_DIR,
    "07_feature_consistency.csv"
  )
)

if (!length(pseudobulk_list)) {
  stop(
    "没有生成任何GSE163973 pseudobulk。"
  )
}

pb_lengths <- vapply(
  pseudobulk_list,
  length,
  integer(1)
)

if (
  length(
    unique(
      pb_lengths
    )
  ) != 1L
) {
  stop(
    "GSE163973 pseudobulk向量长度不一致。"
  )
}

disease_pb_raw <- do.call(
  cbind,
  pseudobulk_list
)

rownames(disease_pb_raw) <-
  make.unique(
    gene_reference_symbol
  )

disease_pb <- collapse_duplicate_symbols(
  disease_pb_raw,
  gene_reference_symbol
)

saveRDS(
  list(
    counts = disease_pb,
    sample_info = pb_info,
    gene_symbol = rownames(
      disease_pb
    ),
    source = "GSE163973"
  ),
  file = file.path(
    OBJECT_DIR,
    "GSE163973_patient_celltype_pseudobulk.rds"
  ),
  compress = "gzip"
)

# ----------------------------- 固定签名映射 ------------------------------------
log_msg(
  "读取正常伤口参考pseudobulk与固定模型"
)

ref_pb_obj <- readRDS(
  REFERENCE_PB_RDS
)

ref_model <- readRDS(
  REFERENCE_MODEL_RDS
)

ref_counts_raw <- ref_pb_obj$counts
ref_gene_symbol <- ref_pb_obj$gene_symbol
ref_sample_info <- ref_pb_obj$sample_info

ref_design <- get_reference_design(
  ref_counts_raw,
  ref_sample_info
)

ref_counts <- collapse_duplicate_symbols(
  ref_counts_raw,
  ref_gene_symbol
)

available_disease_celltypes <- unique(
  pb_info$main_cell_type
)

available_reference_celltypes <- unique(
  ref_design$main_cell_type
)

candidate_celltypes <- intersect(
  c(
    "Fibroblast",
    "Keratinocyte",
    "Myeloid",
    "Endothelial",
    "Lymphoid"
  ),
  intersect(
    available_disease_celltypes,
    available_reference_celltypes
  )
)

if (!"Fibroblast" %in% candidate_celltypes) {
  stop(
    "GSE163973未形成可用Fibroblast pseudobulk。"
  )
}

mapping_rows <- list()
gene_coverage_rows <- list()
comparison_rows <- list()
reference_mapping_rows <- list()
calibration_objects <- list()

for (ct in candidate_celltypes) {
  log_msg(
    "固定签名映射：",
    ct
  )

  disease_keys <- pb_info$sample_key[
    pb_info$main_cell_type == ct
  ]

  # 只有6名患者全部有该细胞类型时，才进行正式组间比较。
  if (length(disease_keys) < 6L) {
    warn_msg(
      ct,
      "只有",
      length(disease_keys),
      "个患者pseudobulk；保留评分，但不作为完整6样本比较。"
    )
  }

  ref_idx <- which(
    ref_design$main_cell_type == ct
  )

  disease_idx <- match(
    disease_keys,
    colnames(
      disease_pb
    )
  )

  if (anyNA(disease_idx)) {
    stop(
      ct,
      "疾病pseudobulk列匹配失败。"
    )
  }

  ref_counts_ct <- ref_counts[
    ,
    ref_idx,
    drop = FALSE
  ]

  ref_design_ct <- ref_design[
    ref_idx,
    ,
    drop = FALSE
  ]

  disease_counts_ct <- disease_pb[
    ,
    disease_idx,
    drop = FALSE
  ]

  disease_info_ct <- pb_info[
    match(
      disease_keys,
      pb_info$sample_key
    ),
    ,
    drop = FALSE
  ]

  common_genes <- intersect(
    rownames(
      ref_counts_ct
    ),
    rownames(
      disease_counts_ct
    )
  )

  common_genes <- common_genes[
    !exclude_background_gene(
      common_genes
    )
  ]

  if (length(common_genes) < 5000L) {
    stop(
      ct,
      "参考与疾病共同背景基因少于5000个，实际：",
      length(common_genes)
    )
  }

  stable_signatures <- get_stable_signatures(
    ref_model,
    ct
  )

  signature_coverage <- vapply(
    stable_signatures,
    function(g) {
      sum(
        g %in% common_genes
      )
    },
    integer(1)
  )

  gene_coverage_rows[[ct]] <- data.frame(
    cell_type = ct,
    n_common_background_genes = length(
      common_genes
    ),
    stage = names(
      signature_coverage
    ),
    n_stable_signature_genes =
      vapply(
        stable_signatures,
        length,
        integer(1)
      ),
    n_signature_genes_in_common =
      as.integer(
        signature_coverage
      ),
    coverage_fraction =
      as.integer(
        signature_coverage
      ) /
      vapply(
        stable_signatures,
        length,
        integer(1)
      ),
    stringsAsFactors = FALSE
  )

  if (any(signature_coverage < 15L)) {
    stop(
      ct,
      "至少一个阶段在共同背景中少于15个稳定基因。"
    )
  }

  ref_expr <- log_cpm(
    ref_counts_ct[
      common_genes,
      ,
      drop = FALSE
    ]
  )

  disease_expr <- log_cpm(
    disease_counts_ct[
      common_genes,
      ,
      drop = FALSE
    ]
  )

  ref_scores <- score_samples_by_rank(
    ref_expr,
    stable_signatures
  )

  disease_scores <- score_samples_by_rank(
    disease_expr,
    stable_signatures
  )

  standardized <- standardize_scores(
    ref_scores,
    disease_scores
  )

  ref_z <- standardized$ref_z
  disease_z <- standardized$new_z

  centroids <- calculate_reference_centroids(
    ref_z,
    ref_design_ct
  )

  reference_lodo_distance <-
    calculate_lodo_reference_distance(
      ref_z,
      ref_design_ct
    )

  reference_distance_limit <-
    as.numeric(
      stats::quantile(
        reference_lodo_distance,
        probs = 0.95,
        na.rm = TRUE,
        names = FALSE,
        type = 8
      )
    )

  if (
    !is.finite(
      reference_distance_limit
    ) ||
      reference_distance_limit <= 0
  ) {
    reference_distance_limit <- max(
      reference_lodo_distance,
      na.rm = TRUE
    )
  }

  if (
    !is.finite(
      reference_distance_limit
    ) ||
      reference_distance_limit <= 0
  ) {
    reference_distance_limit <- 1
  }

  disease_map <- map_to_reference(
    disease_z,
    centroids,
    reference_distance_limit
  )

  disease_map <- add_score_contrasts(
    disease_map,
    disease_z
  )

  disease_map$cell_type <- ct

  disease_map <- merge(
    disease_map,
    disease_info_ct[
      ,
      c(
        "sample_key",
        "sample_id",
        "group",
        "n_cells",
        "pseudobulk_library_size",
        "pseudobulk_detected_genes"
      )
    ],
    by = "sample_key",
    all.x = TRUE,
    sort = FALSE
  )

  disease_map <- disease_map[
    match(
      rownames(
        disease_z
      ),
      disease_map$sample_key
    ),
    ,
    drop = FALSE
  ]

  mapping_rows[[ct]] <- disease_map

  # 正常参考样本也用同一坐标体系映射，用于图形与阈值审计。
  ref_map <- map_to_reference(
    ref_z,
    centroids,
    reference_distance_limit
  )

  ref_map <- add_score_contrasts(
    ref_map,
    ref_z
  )

  ref_map$cell_type <- ct
  ref_map$sample_id <- ref_design_ct$sample_id
  ref_map$donor <- ref_design_ct$donor
  ref_map$true_stage <- ref_design_ct$condition
  ref_map$reference_lodo_true_stage_distance <-
    reference_lodo_distance

  reference_mapping_rows[[ct]] <- ref_map

  calibration_objects[[ct]] <- list(
    stable_signatures = stable_signatures,
    common_background_genes = common_genes,
    score_center = standardized$center,
    score_scale = standardized$scale,
    stage_centroids = centroids,
    reference_distance_limit =
      reference_distance_limit
  )

  # 患者级比较
  metrics_to_compare <- c(
    "z_Skin",
    "z_Wound1",
    "z_Wound7",
    "z_Wound30",
    "early_state_persistence",
    "remodeling_completion",
    "wound_activation",
    "ordinal_position",
    "off_trajectory_distance",
    "off_trajectory_ratio",
    "mixed_state_entropy"
  )

  for (metric in metrics_to_compare) {
    tst <- exact_permutation_test(
      disease_map[[metric]],
      disease_map$group
    )

    comparison_rows[[length(comparison_rows) + 1L]] <-
      data.frame(
        cell_type = ct,
        metric = metric,
        primary_endpoint =
          ct == "Fibroblast" &
          metric %in% c(
            "early_state_persistence",
            "off_trajectory_ratio"
          ),
        n_keloid = sum(
          disease_map$group == "keloid" &
            is.finite(
              disease_map[[metric]]
            )
        ),
        n_normal_scar = sum(
          disease_map$group == "normal_scar" &
            is.finite(
              disease_map[[metric]]
            )
        ),
        mean_keloid = tst$mean_a,
        mean_normal_scar = tst$mean_b,
        median_keloid = tst$median_a,
        median_normal_scar = tst$median_b,
        mean_difference = tst$mean_difference,
        hodges_lehmann_shift =
          tst$hodges_lehmann_shift,
        hedges_g = tst$hedges_g,
        cliffs_delta = tst$cliffs_delta,
        exact_permutation_p =
          tst$exact_permutation_p,
        effect_strength = effect_strength(
          tst$hedges_g,
          tst$cliffs_delta
        ),
        stringsAsFactors = FALSE
      )
  }
}

mapping_results <- do.call(
  rbind,
  mapping_rows
)

reference_mapping <- do.call(
  rbind,
  reference_mapping_rows
)

gene_coverage <- do.call(
  rbind,
  gene_coverage_rows
)

comparison_results <- do.call(
  rbind,
  comparison_rows
)

rownames(mapping_results) <- NULL
rownames(reference_mapping) <- NULL
rownames(gene_coverage) <- NULL
rownames(comparison_results) <- NULL

comparison_results$BH_FDR <- stats::p.adjust(
  comparison_results$exact_permutation_p,
  method = "BH"
)

safe_write_csv(
  gene_coverage,
  file.path(
    REPORT_DIR,
    "08_stable_signature_gene_coverage.csv"
  )
)

safe_write_csv(
  mapping_results,
  file.path(
    REPORT_DIR,
    "09_GSE163973_patient_wound_state_scores.csv"
  )
)

safe_write_csv(
  reference_mapping,
  file.path(
    REPORT_DIR,
    "10_reference_sample_calibration_scores.csv"
  )
)

safe_write_csv(
  comparison_results,
  file.path(
    REPORT_DIR,
    "11_patient_level_group_comparisons.csv"
  )
)

saveRDS(
  list(
    disease_scores = mapping_results,
    reference_scores = reference_mapping,
    comparisons = comparison_results,
    calibration = calibration_objects,
    disease_pseudobulk = list(
      counts = disease_pb,
      sample_info = pb_info
    )
  ),
  file = file.path(
    OBJECT_DIR,
    "GSE163973_fixed_wound_state_mapping.rds"
  ),
  compress = "gzip"
)

# ----------------------------- 图形 --------------------------------------------
fib <- mapping_results[
  mapping_results$cell_type == "Fibroblast",
  ,
  drop = FALSE
]

ref_fib <- reference_mapping[
  reference_mapping$cell_type == "Fibroblast",
  ,
  drop = FALSE
]

if (nrow(fib) != 6L) {
  warn_msg(
    "Fibroblast映射样本不是6个，实际：",
    nrow(fib)
  )
}

# 1. Fibroblast患者级指标
pdf(
  file.path(
    FIG_DIR,
    "Figure_fibroblast_patient_level_metrics.pdf"
  ),
  width = 11,
  height = 8
)

op <- par(
  no.readonly = TRUE
)

par(
  mfrow = c(2, 3),
  mar = c(6, 4, 3, 1)
)

metrics_plot <- c(
  "early_state_persistence",
  "remodeling_completion",
  "off_trajectory_ratio",
  "ordinal_position",
  "z_Wound1",
  "z_Wound7"
)

metric_titles <- c(
  "Early-state persistence",
  "Remodeling completion",
  "Off-trajectory ratio",
  "Ordinal wound-state position",
  "D1-like score",
  "D7-like score"
)

for (k in seq_along(metrics_plot)) {
  metric <- metrics_plot[k]

  grp_num <- ifelse(
    fib$group == "normal_scar",
    1,
    2
  )

  plot(
    jitter(
      grp_num,
      amount = 0.06
    ),
    fib[[metric]],
    pch = ifelse(
      fib$group == "normal_scar",
      1,
      16
    ),
    xaxt = "n",
    xlim = c(0.5, 2.5),
    xlab = "",
    ylab = metric,
    main = metric_titles[k]
  )

  axis(
    1,
    at = c(1, 2),
    labels = c(
      "Normal scar",
      "Keloid"
    ),
    las = 2
  )

  segments(
    x0 = 0.85,
    y0 = stats::median(
      fib[[metric]][
        fib$group == "normal_scar"
      ]
    ),
    x1 = 1.15,
    y1 = stats::median(
      fib[[metric]][
        fib$group == "normal_scar"
      ]
    ),
    lwd = 2
  )

  segments(
    x0 = 1.85,
    y0 = stats::median(
      fib[[metric]][
        fib$group == "keloid"
      ]
    ),
    x1 = 2.15,
    y1 = stats::median(
      fib[[metric]][
        fib$group == "keloid"
      ]
    ),
    lwd = 2
  )

  text(
    jitter(
      grp_num,
      amount = 0.06
    ),
    fib[[metric]],
    labels = fib$sample_id,
    pos = 3,
    cex = 0.7
  )
}

par(op)
dev.off()

# 2. Fibroblast四阶段z-score热图
fib_z <- as.matrix(
  fib[
    ,
    c(
      "z_Skin",
      "z_Wound1",
      "z_Wound7",
      "z_Wound30"
    )
  ]
)

rownames(fib_z) <- fib$sample_id
colnames(fib_z) <- c(
  "Skin",
  "D1",
  "D7",
  "D30"
)

pdf(
  file.path(
    FIG_DIR,
    "Figure_fibroblast_stage_score_heatmap.pdf"
  ),
  width = 7,
  height = 6
)

par(
  mar = c(5, 7, 3, 2)
)

image(
  x = seq_len(
    ncol(fib_z)
  ),
  y = seq_len(
    nrow(fib_z)
  ),
  z = t(fib_z),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "GSE163973 fibroblast fixed wound-state scores"
)

axis(
  1,
  at = seq_len(
    ncol(fib_z)
  ),
  labels = colnames(
    fib_z
  )
)

axis(
  2,
  at = seq_len(
    nrow(fib_z)
  ),
  labels = paste(
    fib$sample_id,
    fib$group,
    sep = " | "
  ),
  las = 2,
  cex.axis = 0.8
)

box()
dev.off()

# 3. Fibroblast四维状态空间的参考PCA投影
ref_z_matrix <- as.matrix(
  ref_fib[
    ,
    c(
      "z_Skin",
      "z_Wound1",
      "z_Wound7",
      "z_Wound30"
    )
  ]
)

colnames(ref_z_matrix) <- c(
  "Skin",
  "D1",
  "D7",
  "D30"
)

disease_z_matrix <- as.matrix(
  fib[
    ,
    c(
      "z_Skin",
      "z_Wound1",
      "z_Wound7",
      "z_Wound30"
    )
  ]
)

colnames(disease_z_matrix) <- c(
  "Skin",
  "D1",
  "D7",
  "D30"
)

pca <- stats::prcomp(
  ref_z_matrix,
  center = TRUE,
  scale. = FALSE
)

ref_pc <- predict(
  pca,
  newdata = ref_z_matrix
)

disease_pc <- predict(
  pca,
  newdata = disease_z_matrix
)

pdf(
  file.path(
    FIG_DIR,
    "Figure_fibroblast_reference_trajectory_projection.pdf"
  ),
  width = 8,
  height = 7
)

stage_pch <- c(
  Skin = 1,
  Wound1 = 2,
  Wound7 = 3,
  Wound30 = 4
)

plot(
  ref_pc[, 1],
  ref_pc[, 2],
  pch = stage_pch[
    ref_fib$true_stage
  ],
  cex = 1.1,
  xlab = "Reference PC1",
  ylab = "Reference PC2",
  main = "Fibroblast projection onto normal wound-state space"
)

points(
  disease_pc[
    fib$group == "normal_scar",
    1
  ],
  disease_pc[
    fib$group == "normal_scar",
    2
  ],
  pch = 15,
  cex = 1.2
)

points(
  disease_pc[
    fib$group == "keloid",
    1
  ],
  disease_pc[
    fib$group == "keloid",
    2
  ],
  pch = 16,
  cex = 1.2
)

text(
  disease_pc[, 1],
  disease_pc[, 2],
  labels = fib$sample_id,
  pos = 3,
  cex = 0.75
)

legend(
  "topright",
  legend = c(
    "Reference Skin",
    "Reference D1",
    "Reference D7",
    "Reference D30",
    "Normal scar",
    "Keloid"
  ),
  pch = c(
    1,
    2,
    3,
    4,
    15,
    16
  ),
  bty = "n",
  cex = 0.8
)

dev.off()

# ----------------------------- 决策闸门 ----------------------------------------
primary <- comparison_results[
  comparison_results$cell_type == "Fibroblast" &
    comparison_results$metric %in% c(
      "early_state_persistence",
      "off_trajectory_ratio"
    ),
  ,
  drop = FALSE
]

get_primary <- function(metric) {
  x <- primary[
    primary$metric == metric,
    ,
    drop = FALSE
  ]

  if (nrow(x) != 1L) {
    stop(
      "无法唯一获得主要终点：",
      metric
    )
  }

  x
}

early <- get_primary(
  "early_state_persistence"
)

offtraj <- get_primary(
  "off_trajectory_ratio"
)

early_positive <- (
  early$mean_difference > 0 &&
    (
      (
        is.finite(
          early$hedges_g
        ) &&
          early$hedges_g >= 0.8
      ) ||
        (
          is.finite(
            early$cliffs_delta
          ) &&
            early$cliffs_delta >= 0.56
        )
    )
)

offtraj_positive <- (
  offtraj$mean_difference > 0 &&
    (
      (
        is.finite(
          offtraj$hedges_g
        ) &&
          offtraj$hedges_g >= 0.8
      ) ||
        (
          is.finite(
            offtraj$cliffs_delta
          ) &&
            offtraj$cliffs_delta >= 0.56
        )
    )
)

fib_count_pass <- all(
  fib$n_cells >= 100
)

signature_coverage_pass <- all(
  gene_coverage$n_signature_genes_in_common >= 15
)

if (
  fib_count_pass &&
    signature_coverage_pass &&
    early_positive &&
    offtraj_positive
) {
  project_gate <- "PASS_INDEPENDENT_DISEASE_SIGNAL"
  next_step <- paste(
    "进入第二独立瘢痕疙瘩队列GSE181316验证；",
    "主要复制early-state persistence与off-trajectory ratio。"
  )
} else if (
  fib_count_pass &&
    signature_coverage_pass &&
    (
      early_positive ||
        offtraj_positive
    )
) {
  project_gate <- "PASS_PARTIAL_DISEASE_SIGNAL"
  next_step <- paste(
    "进入GSE181316验证，但只保留本队列方向一致的终点为候选；",
    "不能提前宣称完整修复终止失败。"
  )
} else {
  project_gate <- "NO_CLEAR_SIGNAL_IN_GSE163973"
  next_step <- paste(
    "仍可运行GSE181316作最终独立核验；",
    "若第二队列同样无方向，则停止伤口状态持续主线。"
  )
}

decision_lines <- c(
  "第五阶段：GSE163973固定伤口状态映射结论（V1.2修正版）",
  "============================================================",
  paste0(
    "运行时间：",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  paste0(
    "项目闸门：",
    project_gate
  ),
  "",
  "主要终点（Fibroblast，keloid减normal scar）：",
  paste0(
    "- early-state persistence mean difference = ",
    round(
      early$mean_difference,
      4
    ),
    "; Hedges g = ",
    round(
      early$hedges_g,
      3
    ),
    "; Cliff delta = ",
    round(
      early$cliffs_delta,
      3
    ),
    "; exact P = ",
    round(
      early$exact_permutation_p,
      3
    )
  ),
  paste0(
    "- off-trajectory ratio mean difference = ",
    round(
      offtraj$mean_difference,
      4
    ),
    "; Hedges g = ",
    round(
      offtraj$hedges_g,
      3
    ),
    "; Cliff delta = ",
    round(
      offtraj$cliffs_delta,
      3
    ),
    "; exact P = ",
    round(
      offtraj$exact_permutation_p,
      3
    )
  ),
  "",
  paste0(
    "6例Fibroblast是否均≥100细胞：",
    fib_count_pass
  ),
  paste0(
    "稳定签名共同基因覆盖是否合格：",
    signature_coverage_pass
  ),
  "",
  "下一步：",
  next_step,
  "",
  "解释边界：",
  "- n=3 vs 3时exact permutation双侧P值分辨率有限，不以P<0.05作为继续条件。",
  "- 主要依据患者级效应量、方向一致性及下一独立队列复制。",
  "- normal scar并非未受伤皮肤，因此结果解释为keloid相对成熟瘢痕的状态偏移。",
  "- ordinal position只代表相对正常伤口状态，不代表真实术后日数。",
  "- 本阶段没有使用GSE163973重新选择任何签名基因。"
)

safe_write_lines(
  decision_lines,
  file.path(
    REPORT_DIR,
    "12_STAGE5_DECISION.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    REPORT_DIR,
    "13_SESSION_INFO.txt"
  )
)

safe_write_lines(
  c(
    paste0(
      "STAGE5_COMPLETED=",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    paste0(
      "PROJECT_GATE=",
      project_gate
    ),
    paste0(
      "FIBROBLAST_COUNT_PASS=",
      fib_count_pass
    ),
    paste0(
      "SIGNATURE_COVERAGE_PASS=",
      signature_coverage_pass
    ),
    paste0(
      "EARLY_PERSISTENCE_POSITIVE=",
      early_positive
    ),
    paste0(
      "OFF_TRAJECTORY_POSITIVE=",
      offtraj_positive
    )
  ),
  file.path(
    REPORT_DIR,
    "STAGE5_COMPLETED.txt"
  )
)

# 删除临时解包内容；本地objects保留。
if (dir.exists(WORK_DIR)) {
  unlink(
    WORK_DIR,
    recursive = TRUE,
    force = TRUE
  )
}

# ----------------------------- 自动打包 ----------------------------------------
package_dir <- file.path(
  STAGE_DIR,
  "_package_for_review"
)

if (dir.exists(package_dir)) {
  unlink(
    package_dir,
    recursive = TRUE,
    force = TRUE
  )
}

dir.create(
  package_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

file.copy(
  REPORT_DIR,
  package_dir,
  recursive = TRUE
)

file.copy(
  FIG_DIR,
  package_dir,
  recursive = TRUE
)

create_zip_windows <- function(
  source_dir,
  zip_path
) {
  if (file.exists(zip_path)) {
    file.remove(zip_path)
  }

  esc <- function(x) {
    gsub(
      "'",
      "''",
      normalize_slash(x),
      fixed = TRUE
    )
  }

  cmd <- paste0(
    "$ErrorActionPreference='Stop'; ",
    "$items=Get-ChildItem -LiteralPath '",
    esc(source_dir),
    "'; ",
    "Compress-Archive -Path $items.FullName ",
    "-DestinationPath '",
    esc(zip_path),
    "' -Force"
  )

  try(
    system2(
      "powershell.exe",
      c(
        "-NoProfile",
        "-ExecutionPolicy",
        "Bypass",
        "-Command",
        cmd
      ),
      stdout = TRUE,
      stderr = TRUE
    ),
    silent = TRUE
  )

  file.exists(zip_path) &&
    file.info(zip_path)$size > 0
}

zip_ok <- create_zip_windows(
  package_dir,
  PACKAGE_ZIP
)

if (!zip_ok) {
  if (file.exists(PACKAGE_TARGZ)) {
    file.remove(PACKAGE_TARGZ)
  }

  old_wd <- getwd()
  setwd(STAGE_DIR)

  try(
    utils::tar(
      tarfile = PACKAGE_TARGZ,
      files = basename(package_dir),
      compression = "gzip",
      tar = "internal"
    ),
    silent = TRUE
  )

  setwd(old_wd)
}

unlink(
  package_dir,
  recursive = TRUE,
  force = TRUE
)

log_msg("第五阶段完成")

cat("\n============================================================\n")
cat("第五阶段运行完成。\n")
cat(
  "本地映射对象保存在：",
  file.path(
    OBJECT_DIR,
    "GSE163973_fixed_wound_state_mapping.rds"
  ),
  "\n",
  sep = ""
)

if (file.exists(PACKAGE_ZIP)) {
  cat(
    "请上传：",
    PACKAGE_ZIP,
    "\n",
    sep = ""
  )
} else if (file.exists(PACKAGE_TARGZ)) {
  cat(
    "请上传备用包：",
    PACKAGE_TARGZ,
    "\n",
    sep = ""
  )
} else {
  cat(
    "自动压缩失败，请手动压缩report和figures文件夹。\n"
  )
}

cat("============================================================\n")
