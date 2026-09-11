# =============================================================================
# 项目：瘢痕疙瘩“伤口状态持续与修复终止失败”
# 第六阶段A：GSE181316未过滤10x矩阵的cell calling、QC与高特异性成纤维细胞提取
#
# 根目录：
#   C:/Users/33652/Desktop/LZZ/1、yssj
#
# 重要说明：
#   GSE181316每个样本有6,794,880个barcode，这是未过滤的10x全barcode矩阵，
#   绝不能把全部barcode当作细胞。本阶段先完成emptyDrops cell calling，
#   再做保守QC和预设标志物规则的成纤维细胞提取。
#
# 本阶段不进行：
#   - 不计算任何伤口状态评分；
#   - 不比较keloid与normal scar；
#   - 不根据结果修改成纤维细胞规则；
#   - 不把keloid_3L与keloid_3R当作两个独立患者进行统计。
#
# 本阶段目的：
#   1. 使用DropletUtils::emptyDrops识别真实细胞；
#   2. 进行nCount、nFeature和线粒体比例QC；
#   3. 使用预先固定的高特异性标志物规则提取fibroblast-like细胞；
#   4. 同时生成较宽松的broad fibroblast-like敏感性集合；
#   5. 生成每个样本的成纤维细胞pseudobulk；
#   6. 自动生成第六阶段A检查包，供人工复核后再运行伤口状态映射。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANG = "en")

ROOT_DIR <- "C:/Users/33652/Desktop/LZZ/1、yssj"

RAW_TAR <- file.path(
  ROOT_DIR,
  "GSE181316_RAW.tar"
)

STAGE_DIR <- file.path(
  ROOT_DIR,
  "06A_STAGE6A_GSE181316_CELL_CALLING"
)

REPORT_DIR <- file.path(
  STAGE_DIR,
  "report"
)

FIG_DIR <- file.path(
  STAGE_DIR,
  "figures"
)

OBJECT_DIR <- file.path(
  STAGE_DIR,
  "objects"
)

WORK_DIR <- file.path(
  STAGE_DIR,
  "_temporary_work"
)

PACKAGE_ZIP <- file.path(
  ROOT_DIR,
  "第六阶段A_GSE181316细胞识别与成纤维细胞提取检查包.zip"
)

PACKAGE_TARGZ <- file.path(
  ROOT_DIR,
  "第六阶段A_GSE181316细胞识别与成纤维细胞提取检查包.tar.gz"
)

if (dir.exists(STAGE_DIR)) {
  unlink(
    STAGE_DIR,
    recursive = TRUE,
    force = TRUE
  )
}

dir.create(
  REPORT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  FIG_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  OBJECT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  WORK_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

LOG_FILE <- file.path(
  REPORT_DIR,
  "00_stage6A_log.txt"
)

WARNING_FILE <- file.path(
  REPORT_DIR,
  "00_stage6A_warnings.txt"
)

log_msg <- function(...) {
  msg <- paste0(
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    " | ",
    paste(
      ...,
      collapse = " "
    )
  )

  cat(
    msg,
    "\n"
  )

  cat(
    msg,
    "\n",
    file = LOG_FILE,
    append = TRUE
  )
}

warn_msg <- function(...) {
  msg <- paste(
    ...,
    collapse = " "
  )

  log_msg(
    "WARNING:",
    msg
  )

  cat(
    msg,
    "\n",
    file = WARNING_FILE,
    append = TRUE
  )
}

safe_write_csv <- function(
  x,
  path
) {
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

safe_write_lines <- function(
  x,
  path
) {
  tryCatch(
    writeLines(
      x,
      path,
      useBytes = TRUE
    ),
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
  normalizePath(
    x,
    winslash = "/",
    mustWork = FALSE
  )
}

# ----------------------------- 依赖包 ------------------------------------------
if (!requireNamespace(
  "Matrix",
  quietly = TRUE
)) {
  install.packages(
    "Matrix",
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace(
  "BiocManager",
  quietly = TRUE
)) {
  install.packages(
    "BiocManager",
    repos = "https://cloud.r-project.org"
  )
}

if (!requireNamespace(
  "DropletUtils",
  quietly = TRUE
)) {
  BiocManager::install(
    "DropletUtils",
    ask = FALSE,
    update = FALSE
  )
}

required_packages <- c(
  "Matrix",
  "DropletUtils",
  "S4Vectors"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages)) {
  stop(
    "无法加载必要R包：",
    paste(
      missing_packages,
      collapse = ", "
    ),
    "。本阶段不使用未经验证的替代cell calling方法。"
  )
}

# ----------------------------- 固定样本设计 ------------------------------------
sample_design <- data.frame(
  sample_id = c(
    "skin_7",
    "keloid_1",
    "keloid_2",
    "keloid_3L",
    "keloid_3R",
    "scar_1",
    "scar_2",
    "scar_3"
  ),
  group = c(
    "healthy_skin",
    "keloid",
    "keloid",
    "keloid",
    "keloid",
    "normal_scar",
    "normal_scar",
    "normal_scar"
  ),
  patient_id = c(
    "H7",
    "K1",
    "K2",
    "K3",
    "K3",
    "S1",
    "S2",
    "S3"
  ),
  independent_patient = c(
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    TRUE,
    TRUE,
    TRUE
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  sample_design,
  file.path(
    REPORT_DIR,
    "01_fixed_sample_design.csv"
  )
)

# ----------------------------- 固定标志物 --------------------------------------
fibro_core <- c(
  "COL1A1",
  "COL1A2",
  "COL3A1",
  "DCN",
  "LUM",
  "COL6A1",
  "COL6A2",
  "SPARC"
)

fibro_support <- c(
  "PDGFRA",
  "DPT",
  "C7",
  "COL14A1",
  "COL12A1",
  "FBLN1",
  "FBLN2",
  "C1S",
  "C1R"
)

immune_markers <- c(
  "PTPRC",
  "TYROBP",
  "LST1",
  "FCER1G",
  "CTSS",
  "CD3D",
  "CD3E",
  "MS4A1",
  "NKG7"
)

epithelial_markers <- c(
  "EPCAM",
  "KRT5",
  "KRT14",
  "KRT1",
  "KRT10",
  "KRT17",
  "KRT19",
  "KRT8",
  "KRT18",
  "SFN"
)

endothelial_markers <- c(
  "PECAM1",
  "VWF",
  "KDR",
  "EMCN",
  "ENG",
  "RAMP2",
  "CLDN5",
  "ESAM"
)

neural_markers <- c(
  "SOX10",
  "S100B",
  "PLP1",
  "MPZ",
  "SLC1A3",
  "NGFR",
  "FABP7",
  "PMP22"
)

marker_sets <- list(
  fibro_core = fibro_core,
  fibro_support = fibro_support,
  immune = immune_markers,
  epithelial = epithelial_markers,
  endothelial = endothelial_markers,
  neural = neural_markers
)

marker_definition <- do.call(
  rbind,
  lapply(
    names(marker_sets),
    function(nm) {
      data.frame(
        marker_set = nm,
        gene = marker_sets[[nm]],
        stringsAsFactors = FALSE
      )
    }
  )
)

safe_write_csv(
  marker_definition,
  file.path(
    REPORT_DIR,
    "02_fixed_marker_definition.csv"
  )
)

# ----------------------------- 基础函数 ----------------------------------------
read_mtx_gz <- function(path) {
  con <- if (
    grepl(
      "\\.gz$",
      path,
      ignore.case = TRUE
    )
  ) {
    gzfile(
      path,
      "rt"
    )
  } else {
    file(
      path,
      "rt"
    )
  }

  on.exit(
    close(con),
    add = TRUE
  )

  x <- Matrix::readMM(
    con
  )

  if (!inherits(
    x,
    "CsparseMatrix"
  )) {
    x <- methods::as(
      x,
      "CsparseMatrix"
    )
  }

  x
}

read_noheader_table <- function(
  path,
  sep = "\t"
) {
  con <- if (
    grepl(
      "\\.gz$",
      path,
      ignore.case = TRUE
    )
  ) {
    gzfile(
      path,
      "rt"
    )
  } else {
    file(
      path,
      "rt"
    )
  }

  on.exit(
    close(con),
    add = TRUE
  )

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

read_selected_lines <- function(
  path,
  selected_indices,
  chunk_size = 200000L
) {
  selected_indices <- as.integer(
    selected_indices
  )

  if (!length(selected_indices)) {
    return(character(0))
  }

  if (
    anyNA(selected_indices) ||
      any(selected_indices < 1L) ||
      anyDuplicated(selected_indices)
  ) {
    stop(
      "selected_indices必须是唯一的正整数。"
    )
  }

  ord <- order(
    selected_indices
  )

  sorted_idx <- selected_indices[
    ord
  ]

  result_sorted <- character(
    length(sorted_idx)
  )

  con <- if (
    grepl(
      "\\.gz$",
      path,
      ignore.case = TRUE
    )
  ) {
    gzfile(
      path,
      "rt"
    )
  } else {
    file(
      path,
      "rt"
    )
  }

  on.exit(
    close(con),
    add = TRUE
  )

  global_start <- 1L
  target_pointer <- 1L

  repeat {
    lines <- readLines(
      con,
      n = chunk_size,
      warn = FALSE
    )

    if (!length(lines)) {
      break
    }

    global_end <- global_start +
      length(lines) -
      1L

    while (
      target_pointer <= length(sorted_idx) &&
        sorted_idx[target_pointer] <= global_end
    ) {
      if (
        sorted_idx[target_pointer] >= global_start
      ) {
        local_position <-
          sorted_idx[target_pointer] -
          global_start +
          1L

        result_sorted[target_pointer] <-
          lines[local_position]
      }

      target_pointer <- target_pointer + 1L
    }

    if (
      target_pointer > length(sorted_idx)
    ) {
      break
    }

    global_start <- global_end + 1L
  }

  if (
    target_pointer <= length(sorted_idx) ||
      any(!nzchar(result_sorted))
  ) {
    stop(
      "未能从barcodes文件中读取全部选定行。"
    )
  }

  result <- character(
    length(result_sorted)
  )

  result[ord] <- result_sorted

  result
}

extract_sample_id <- function(
  matrix_filename
) {
  x <- basename(
    matrix_filename
  )

  x <- sub(
    "^GSM[0-9]+_",
    "",
    x
  )

  x <- sub(
    "_matrix\\.mtx\\.gz$",
    "",
    x,
    ignore.case = TRUE
  )

  x
}

extract_triplet_from_outer <- function(
  outer_dir,
  sample_id
) {
  files <- list.files(
    outer_dir,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  )

  escape_regex <- function(x) {
    gsub(
      "([][{}()+*^$|\\\\?.])",
      "\\\\\\1",
      x
    )
  }

  sid_regex <- escape_regex(
    sample_id
  )

  matrix_path <- files[
    grepl(
      paste0(
        "_",
        sid_regex,
        "_matrix\\.mtx\\.gz$"
      ),
      basename(files),
      ignore.case = TRUE,
      perl = TRUE
    )
  ]

  barcode_path <- files[
    grepl(
      paste0(
        "_",
        sid_regex,
        "_barcodes\\.tsv\\.gz$"
      ),
      basename(files),
      ignore.case = TRUE,
      perl = TRUE
    )
  ]

  feature_path <- files[
    grepl(
      paste0(
        "_",
        sid_regex,
        "_features\\.tsv\\.gz$"
      ),
      basename(files),
      ignore.case = TRUE,
      perl = TRUE
    )
  ]

  if (
    length(matrix_path) != 1L ||
      length(barcode_path) != 1L ||
      length(feature_path) != 1L
  ) {
    stop(
      "无法唯一定位样本",
      sample_id,
      "的matrix/barcodes/features文件。"
    )
  }

  list(
    matrix = matrix_path,
    barcodes = barcode_path,
    features = feature_path
  )
}

collapse_vector_by_symbol <- function(
  value,
  gene_symbol
) {
  gene_symbol <- trimws(
    as.character(
      gene_symbol
    )
  )

  valid <- (
    !is.na(
      gene_symbol
    ) &
      nzchar(
        gene_symbol
      ) &
      is.finite(
        value
      )
  )

  sums <- rowsum(
    matrix(
      value[valid],
      ncol = 1
    ),
    group = gene_symbol[valid],
    reorder = FALSE
  )

  out <- as.numeric(
    sums[, 1]
  )

  names(out) <- rownames(
    sums
  )

  out
}

marker_count_matrix <- function(
  counts,
  gene_symbol,
  markers
) {
  gene_upper <- toupper(
    trimws(
      as.character(
        gene_symbol
      )
    )
  )

  markers <- toupper(
    markers
  )

  out <- matrix(
    0,
    nrow = length(markers),
    ncol = ncol(counts),
    dimnames = list(
      markers,
      colnames(counts)
    )
  )

  for (i in seq_along(markers)) {
    idx <- which(
      gene_upper == markers[i]
    )

    if (length(idx)) {
      out[i, ] <- Matrix::colSums(
        counts[
          idx,
          ,
          drop = FALSE
        ]
      )
    }
  }

  out
}

marker_metrics <- function(
  marker_counts,
  library_size
) {
  detected <- colSums(
    marker_counts > 0
  )

  normalized <- log1p(
    sweep(
      marker_counts,
      2,
      pmax(
        library_size,
        1
      ),
      "/"
    ) * 10000
  )

  score <- colMeans(
    normalized
  )

  list(
    detected = as.numeric(
      detected
    ),
    score = as.numeric(
      score
    )
  )
}

safe_quantile <- function(
  x,
  probs
) {
  if (
    !length(x) ||
      all(
        is.na(x)
      )
  ) {
    return(
      rep(
        NA_real_,
        length(probs)
      )
    )
  }

  as.numeric(
    stats::quantile(
      x,
      probs = probs,
      na.rm = TRUE,
      names = FALSE,
      type = 8
    )
  )
}

# ----------------------------- 输入检查 ----------------------------------------
log_msg(
  "第六阶段A开始"
)

if (!file.exists(RAW_TAR)) {
  stop(
    "未找到：",
    RAW_TAR
  )
}

outer_dir <- file.path(
  WORK_DIR,
  "GSE181316_RAW"
)

dir.create(
  outer_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

log_msg(
  "解包GSE181316_RAW.tar"
)

utils::untar(
  RAW_TAR,
  exdir = outer_dir
)

matrix_files <- list.files(
  outer_dir,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "_matrix\\.mtx\\.gz$",
  ignore.case = TRUE
)

detected_sample_ids <- vapply(
  matrix_files,
  extract_sample_id,
  character(1)
)

if (
  length(matrix_files) != 8L ||
    !setequal(
      detected_sample_ids,
      sample_design$sample_id
    )
) {
  stop(
    "GSE181316样本结构与预期不一致。实际：",
    paste(
      detected_sample_ids,
      collapse = ", "
    )
  )
}

# ----------------------------- 逐样本处理 --------------------------------------
cell_call_rows <- list()
marker_summary_rows <- list()
fibro_pb_high <- list()
fibro_pb_broad <- list()
object_info_rows <- list()

pdf(
  file.path(
    FIG_DIR,
    "Figure_cell_calling_and_QC_distributions.pdf"
  ),
  width = 10,
  height = 8,
  onefile = TRUE
)

pdf_device_open <- TRUE

for (
  sample_id in sample_design$sample_id
) {
  log_msg(
    "处理样本:",
    sample_id
  )

  design_row <- sample_design[
    sample_design$sample_id == sample_id,
    ,
    drop = FALSE
  ]

  triplet <- extract_triplet_from_outer(
    outer_dir,
    sample_id
  )

  features <- read_noheader_table(
    triplet$features,
    sep = "\t"
  )

  gene_id <- as.character(
    features[[1]]
  )

  gene_symbol <- if (
    ncol(features) >= 2L
  ) {
    as.character(
      features[[2]]
    )
  } else {
    gene_id
  }

  counts_all <- read_mtx_gz(
    triplet$matrix
  )

  if (
    nrow(counts_all) != length(gene_id)
  ) {
    stop(
      sample_id,
      "的features行数与矩阵行数不一致。"
    )
  }

  n_all_barcodes <- ncol(
    counts_all
  )

  totals_all <- Matrix::colSums(
    counts_all
  )

  nonzero_idx <- which(
    totals_all > 0
  )

  if (
    length(nonzero_idx) < 1000L
  ) {
    stop(
      sample_id,
      "非零barcode少于1000，无法合理进行emptyDrops。"
    )
  }

  counts_nonzero <- counts_all[
    ,
    nonzero_idx,
    drop = FALSE
  ]

  totals_nonzero <- totals_all[
    nonzero_idx
  ]

  rm(
    counts_all,
    totals_all
  )

  gc(
    verbose = FALSE
  )

  # 标准emptyDrops cell calling。
  set.seed(
    20260618L +
      match(
        sample_id,
        sample_design$sample_id
      )
  )

  barcode_ranks <- DropletUtils::barcodeRanks(
    counts_nonzero,
    lower = 100
  )

  rank_meta <- S4Vectors::metadata(
    barcode_ranks
  )

  knee <- suppressWarnings(
    as.numeric(
      rank_meta$knee
    )
  )

  inflection <- suppressWarnings(
    as.numeric(
      rank_meta$inflection
    )
  )

  empty_result <- tryCatch(
    DropletUtils::emptyDrops(
      counts_nonzero,
      lower = 100,
      niters = 10000
    ),
    error = function(e) {
      stop(
        sample_id,
        " emptyDrops失败：",
        conditionMessage(e)
      )
    }
  )

  empty_fdr <- as.numeric(
    empty_result$FDR
  )

  called_by_fdr <- (
    !is.na(
      empty_fdr
    ) &
      empty_fdr <= 0.01
  )

  count_retain_threshold <- knee

  if (
    !is.finite(
      count_retain_threshold
    ) ||
      count_retain_threshold <= 0
  ) {
    count_retain_threshold <- inflection
  }

  if (
    !is.finite(
      count_retain_threshold
    ) ||
      count_retain_threshold <= 0
  ) {
    count_retain_threshold <- as.numeric(
      stats::quantile(
        totals_nonzero,
        probs = 0.995,
        na.rm = TRUE,
        names = FALSE,
        type = 8
      )
    )

    warn_msg(
      sample_id,
      "barcodeRanks未返回有效knee/inflection，使用99.5%计数分位数作为高计数保留阈值。"
    )
  }

  called_by_high_count <- (
    totals_nonzero >=
      count_retain_threshold
  )

  called <- (
    called_by_fdr |
      called_by_high_count
  )

  called_idx_nonzero <- which(
    called
  )

  if (
    length(called_idx_nonzero) < 100L
  ) {
    stop(
      sample_id,
      "识别到的细胞少于100个，停止。"
    )
  }

  if (
    length(called_idx_nonzero) > 100000L
  ) {
    stop(
      sample_id,
      "识别到的细胞超过100,000个，疑似cell calling异常。"
    )
  }

  counts_called <- counts_nonzero[
    ,
    called_idx_nonzero,
    drop = FALSE
  ]

  called_original_indices <-
    nonzero_idx[
      called_idx_nonzero
    ]

  rm(
    counts_nonzero,
    empty_result,
    barcode_ranks
  )

  gc(
    verbose = FALSE
  )

  n_count <- Matrix::colSums(
    counts_called
  )

  n_feature <- Matrix::colSums(
    counts_called > 0
  )

  gene_upper <- toupper(
    gene_symbol
  )

  mt_idx <- grepl(
    "^MT-",
    gene_upper
  )

  pct_mt <- if (
    any(mt_idx)
  ) {
    Matrix::colSums(
      counts_called[
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
      ncol(
        counts_called
      )
    )
  }

  qc_keep <- (
    is.finite(
      n_count
    ) &
      is.finite(
        n_feature
      ) &
      is.finite(
        pct_mt
      ) &
      n_count >= 500 &
      n_feature >= 200 &
      n_feature <= 7500 &
      pct_mt <= 20
  )

  qc_keep[
    is.na(
      qc_keep
    )
  ] <- FALSE

  qc_idx_called <- which(
    qc_keep
  )

  if (
    length(qc_idx_called) < 100L
  ) {
    stop(
      sample_id,
      "QC后细胞少于100个，停止。"
    )
  }

  counts_qc <- counts_called[
    ,
    qc_idx_called,
    drop = FALSE
  ]

  qc_original_indices <-
    called_original_indices[
      qc_idx_called
    ]

  selected_barcodes <- read_selected_lines(
    triplet$barcodes,
    qc_original_indices
  )

  if (
    length(selected_barcodes) !=
      ncol(counts_qc)
  ) {
    stop(
      sample_id,
      "读取到的QC barcode数量与矩阵列数不一致。"
    )
  }

  colnames(
    counts_qc
  ) <- selected_barcodes

  n_count_qc <- as.numeric(
    Matrix::colSums(
      counts_qc
    )
  )

  n_feature_qc <- as.numeric(
    Matrix::colSums(
      counts_qc > 0
    )
  )

  pct_mt_qc <- as.numeric(
    pct_mt[
      qc_idx_called
    ]
  )

  marker_results <- list()

  for (
    marker_name in names(marker_sets)
  ) {
    marker_counts <- marker_count_matrix(
      counts_qc,
      gene_symbol,
      marker_sets[[marker_name]]
    )

    marker_results[[marker_name]] <-
      marker_metrics(
        marker_counts,
        n_count_qc
      )
  }

  max_exclusion_score <- pmax(
    marker_results$immune$score,
    marker_results$epithelial$score,
    marker_results$endothelial$score,
    marker_results$neural$score
  )

  broad_fibroblast <- (
    marker_results$fibro_core$detected >= 2 &
      marker_results$fibro_core$score >
      max_exclusion_score
  )

  high_specificity_fibroblast <- (
    marker_results$fibro_core$detected >= 3 &
      (
        marker_results$fibro_support$detected >= 1 |
          marker_results$fibro_core$detected >= 4
      ) &
      marker_results$fibro_core$score >
      max_exclusion_score + 0.15 &
      marker_results$immune$detected <= 1 &
      marker_results$epithelial$detected <= 1 &
      marker_results$endothelial$detected <= 1 &
      marker_results$neural$detected <= 1
  )

  high_specificity_fibroblast[
    is.na(
      high_specificity_fibroblast
    )
  ] <- FALSE

  broad_fibroblast[
    is.na(
      broad_fibroblast
    )
  ] <- FALSE

  # 高特异性集合必须是宽松集合的子集。
  high_specificity_fibroblast <-
    high_specificity_fibroblast &
    broad_fibroblast

  n_fib_high <- sum(
    high_specificity_fibroblast
  )

  n_fib_broad <- sum(
    broad_fibroblast
  )

  if (
    n_fib_high < 50L
  ) {
    warn_msg(
      sample_id,
      "高特异性fibroblast-like细胞少于50：",
      n_fib_high
    )
  }

  if (
    n_fib_broad < 100L
  ) {
    warn_msg(
      sample_id,
      "宽松fibroblast-like细胞少于100：",
      n_fib_broad
    )
  }

  qc_metadata <- data.frame(
    barcode = selected_barcodes,
    original_matrix_column = qc_original_indices,
    nCount = n_count_qc,
    nFeature = n_feature_qc,
    percent_mt = pct_mt_qc,
    fibro_core_detected =
      marker_results$fibro_core$detected,
    fibro_support_detected =
      marker_results$fibro_support$detected,
    immune_detected =
      marker_results$immune$detected,
    epithelial_detected =
      marker_results$epithelial$detected,
    endothelial_detected =
      marker_results$endothelial$detected,
    neural_detected =
      marker_results$neural$detected,
    fibro_core_score =
      marker_results$fibro_core$score,
    fibro_support_score =
      marker_results$fibro_support$score,
    immune_score =
      marker_results$immune$score,
    epithelial_score =
      marker_results$epithelial$score,
    endothelial_score =
      marker_results$endothelial$score,
    neural_score =
      marker_results$neural$score,
    broad_fibroblast_like =
      broad_fibroblast,
    high_specificity_fibroblast_like =
      high_specificity_fibroblast,
    stringsAsFactors = FALSE
  )

  pb_high_raw <- Matrix::rowSums(
    counts_qc[
      ,
      high_specificity_fibroblast,
      drop = FALSE
    ]
  )

  pb_broad_raw <- Matrix::rowSums(
    counts_qc[
      ,
      broad_fibroblast,
      drop = FALSE
    ]
  )

  pb_high <- collapse_vector_by_symbol(
    pb_high_raw,
    gene_symbol
  )

  pb_broad <- collapse_vector_by_symbol(
    pb_broad_raw,
    gene_symbol
  )

  fibro_pb_high[[sample_id]] <- pb_high
  fibro_pb_broad[[sample_id]] <- pb_broad

  saveRDS(
    list(
      sample_id = sample_id,
      group = design_row$group,
      patient_id = design_row$patient_id,
      cell_calling_method =
        "DropletUtils::emptyDrops(FDR<=0.01) OR barcodeRanks knee retention",
      qc_rule =
        "nCount>=500; nFeature>=200 and <=7500; percent.mt<=20",
      gene_id = gene_id,
      gene_symbol = gene_symbol,
      qc_metadata = qc_metadata,
      fibroblast_high_counts = pb_high,
      fibroblast_broad_counts = pb_broad
    ),
    file = file.path(
      OBJECT_DIR,
      paste0(
        "GSE181316_",
        sample_id,
        "_cellcalling_fibroblast.rds"
      )
    ),
    compress = "gzip"
  )

  q_count <- safe_quantile(
    n_count_qc,
    c(
      0.01,
      0.25,
      0.50,
      0.75,
      0.99
    )
  )

  q_feature <- safe_quantile(
    n_feature_qc,
    c(
      0.01,
      0.25,
      0.50,
      0.75,
      0.99
    )
  )

  cell_call_rows[[length(cell_call_rows) + 1L]] <-
    data.frame(
      sample_id = sample_id,
      group = design_row$group,
      patient_id = design_row$patient_id,
      n_all_barcodes = n_all_barcodes,
      n_nonzero_barcodes = length(
        nonzero_idx
      ),
      knee_count_threshold = knee,
      inflection_count_threshold = inflection,
      retained_count_threshold =
        count_retain_threshold,
      n_called_by_FDR01 = sum(
        called_by_fdr
      ),
      n_called_by_high_count = sum(
        called_by_high_count
      ),
      n_called_union = length(
        called_idx_nonzero
      ),
      n_QC_keep = ncol(
        counts_qc
      ),
      QC_retention_from_called =
        ncol(
          counts_qc
        ) /
        length(
          called_idx_nonzero
        ),
      n_fibroblast_high = n_fib_high,
      n_fibroblast_broad = n_fib_broad,
      fibroblast_high_fraction =
        n_fib_high /
        ncol(
          counts_qc
        ),
      fibroblast_broad_fraction =
        n_fib_broad /
        ncol(
          counts_qc
        ),
      nCount_p01 = q_count[1],
      nCount_p25 = q_count[2],
      nCount_median = q_count[3],
      nCount_p75 = q_count[4],
      nCount_p99 = q_count[5],
      nFeature_p01 = q_feature[1],
      nFeature_p25 = q_feature[2],
      nFeature_median = q_feature[3],
      nFeature_p75 = q_feature[4],
      nFeature_p99 = q_feature[5],
      median_percent_mt = stats::median(
        pct_mt_qc
      ),
      stringsAsFactors = FALSE
    )

  subset_names <- c(
    "all_QC_cells",
    "broad_fibroblast_like",
    "high_specificity_fibroblast_like"
  )

  subset_indices <- list(
    rep(
      TRUE,
      ncol(
        counts_qc
      )
    ),
    broad_fibroblast,
    high_specificity_fibroblast
  )

  for (
    i in seq_along(
      subset_names
    )
  ) {
    idx_subset <- subset_indices[[i]]

    if (!any(idx_subset)) {
      next
    }

    marker_summary_rows[[length(marker_summary_rows) + 1L]] <-
      data.frame(
        sample_id = sample_id,
        group = design_row$group,
        subset = subset_names[i],
        n_cells = sum(
          idx_subset
        ),
        median_fibro_core_detected =
          stats::median(
            marker_results$fibro_core$detected[
              idx_subset
            ]
          ),
        median_fibro_support_detected =
          stats::median(
            marker_results$fibro_support$detected[
              idx_subset
            ]
          ),
        median_fibro_core_score =
          stats::median(
            marker_results$fibro_core$score[
              idx_subset
            ]
          ),
        median_max_exclusion_score =
          stats::median(
            max_exclusion_score[
              idx_subset
            ]
          ),
        pct_PTPRC_positive = mean(
          marker_count_matrix(
            counts_qc[
              ,
              idx_subset,
              drop = FALSE
            ],
            gene_symbol,
            "PTPRC"
          ) > 0
        ) * 100,
        pct_EPCAM_positive = mean(
          marker_count_matrix(
            counts_qc[
              ,
              idx_subset,
              drop = FALSE
            ],
            gene_symbol,
            "EPCAM"
          ) > 0
        ) * 100,
        pct_PECAM1_positive = mean(
          marker_count_matrix(
            counts_qc[
              ,
              idx_subset,
              drop = FALSE
            ],
            gene_symbol,
            "PECAM1"
          ) > 0
        ) * 100,
        pct_SOX10_positive = mean(
          marker_count_matrix(
            counts_qc[
              ,
              idx_subset,
              drop = FALSE
            ],
            gene_symbol,
            "SOX10"
          ) > 0
        ) * 100,
        stringsAsFactors = FALSE
      )
  }

  object_info_rows[[length(object_info_rows) + 1L]] <-
    data.frame(
      sample_id = sample_id,
      object_file = paste0(
        "GSE181316_",
        sample_id,
        "_cellcalling_fibroblast.rds"
      ),
      n_QC_cells = ncol(
        counts_qc
      ),
      n_high_specificity_fibroblast =
        n_fib_high,
      n_broad_fibroblast =
        n_fib_broad,
      stringsAsFactors = FALSE
    )

  # 每个样本4页审计图。
  op <- par(
    no.readonly = TRUE
  )

  par(
    mfrow = c(2, 2),
    mar = c(4, 4, 3, 1)
  )

  plot(
    seq_along(
      sort(
        totals_nonzero,
        decreasing = TRUE
      )
    ),
    sort(
      totals_nonzero,
      decreasing = TRUE
    ),
    log = "xy",
    pch = 16,
    cex = 0.2,
    xlab = "Barcode rank",
    ylab = "Total UMI",
    main = paste0(
      sample_id,
      " barcode rank"
    )
  )

  abline(
    h = count_retain_threshold,
    lty = 2
  )

  hist(
    log10(
      n_count_qc + 1
    ),
    breaks = 50,
    xlab = "log10 nCount + 1",
    main = paste0(
      sample_id,
      " called-cell counts"
    )
  )

  hist(
    n_feature_qc,
    breaks = 50,
    xlab = "nFeature",
    main = paste0(
      sample_id,
      " called-cell genes"
    )
  )

  plot(
    marker_results$fibro_core$score,
    max_exclusion_score,
    pch = ifelse(
      high_specificity_fibroblast,
      16,
      1
    ),
    cex = 0.35,
    xlab = "Fibroblast core score",
    ylab = "Maximum exclusion score",
    main = paste0(
      sample_id,
      " marker classification"
    )
  )

  abline(
    a = -0.15,
    b = 1,
    lty = 2
  )

  par(op)

  rm(
    counts_called,
    counts_qc,
    marker_results,
    qc_metadata
  )

  gc(
    verbose = FALSE
  )
}

if (pdf_device_open) {
  dev.off()
}

# ----------------------------- 汇总输出 ----------------------------------------
cell_call_summary <- do.call(
  rbind,
  cell_call_rows
)

marker_summary <- do.call(
  rbind,
  marker_summary_rows
)

object_summary <- do.call(
  rbind,
  object_info_rows
)

rownames(
  cell_call_summary
) <- NULL

rownames(
  marker_summary
) <- NULL

rownames(
  object_summary
) <- NULL

safe_write_csv(
  cell_call_summary,
  file.path(
    REPORT_DIR,
    "03_cell_calling_QC_fibroblast_summary.csv"
  )
)

safe_write_csv(
  marker_summary,
  file.path(
    REPORT_DIR,
    "04_marker_identity_summary.csv"
  )
)

safe_write_csv(
  object_summary,
  file.path(
    REPORT_DIR,
    "05_local_object_summary.csv"
  )
)

saveRDS(
  list(
    high_specificity = fibro_pb_high,
    broad_sensitivity = fibro_pb_broad,
    sample_design = sample_design,
    marker_definition = marker_definition,
    cell_call_summary = cell_call_summary,
    classification_rule = list(
      broad = paste(
        "fibro_core_detected>=2 and",
        "fibro_core_score>max exclusion score"
      ),
      high_specificity = paste(
        "fibro_core_detected>=3;",
        "fibro_support_detected>=1 or fibro_core_detected>=4;",
        "fibro_core_score>max exclusion+0.15;",
        "immune/epithelial/endothelial/neural detected<=1"
      )
    )
  ),
  file = file.path(
    OBJECT_DIR,
    "GSE181316_fibroblast_like_pseudobulk_list.rds"
  ),
  compress = "gzip"
)

# ----------------------------- 决策闸门 ----------------------------------------
main_sample_ids <- sample_design$sample_id[
  sample_design$group %in% c(
    "keloid",
    "normal_scar"
  )
]

main_summary <- cell_call_summary[
  match(
    main_sample_ids,
    cell_call_summary$sample_id
  ),
  ,
  drop = FALSE
]

cell_calling_pass <- all(
  main_summary$n_called_union >= 500 &
    main_summary$n_called_union <= 50000
)

qc_pass <- all(
  main_summary$n_QC_keep >= 300 &
    main_summary$QC_retention_from_called >= 0.50
)

fibro_high_pass <- all(
  main_summary$n_fibroblast_high >= 100
)

fibro_broad_pass <- all(
  main_summary$n_fibroblast_broad >= 200
)

patient_structure_pass <- (
  length(
    unique(
      sample_design$patient_id[
        sample_design$group == "keloid"
      ]
    )
  ) == 3L &&
    length(
      unique(
        sample_design$patient_id[
          sample_design$group == "normal_scar"
        ]
      )
    ) == 3L
)

if (
  cell_calling_pass &&
    qc_pass &&
    fibro_high_pass &&
    fibro_broad_pass &&
    patient_structure_pass
) {
  project_gate <-
    "PASS_CELL_CALLING_AND_FIBROBLAST_EXTRACTION"

  next_step <- paste(
    "进入第六阶段B固定伤口状态复制；",
    "主要终点预先冻结为late-remodeling state，",
    "关键次要终点为ordinal wound-state position。"
  )
} else if (
  cell_calling_pass &&
    qc_pass &&
    fibro_broad_pass &&
    patient_structure_pass
) {
  project_gate <-
    "PASS_BROAD_FIBROBLAST_ONLY"

  next_step <- paste(
    "高特异性集合数量不足；",
    "需人工复核标志物分布，",
    "不能直接修改规则后追求阳性结果。"
  )
} else {
  project_gate <-
    "FAIL_GSE181316_EXTRACTION"

  next_step <- paste(
    "暂停GSE181316伤口状态复制；",
    "先检查cell calling和标志物识别失败原因。"
  )
}

decision_lines <- c(
  "第六阶段A：GSE181316 cell calling与成纤维细胞提取结论",
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
  paste0(
    "7个keloid/normal-scar样本cell calling是否合格：",
    cell_calling_pass
  ),
  paste0(
    "7个样本QC是否合格：",
    qc_pass
  ),
  paste0(
    "7个样本高特异性fibroblast-like是否均≥100：",
    fibro_high_pass
  ),
  paste0(
    "7个样本宽松fibroblast-like是否均≥200：",
    fibro_broad_pass
  ),
  paste0(
    "患者结构是否为3例keloid vs 3例normal scar：",
    patient_structure_pass
  ),
  "",
  "下一步：",
  next_step,
  "",
  "已冻结的第六阶段B复制终点：",
  "- Primary: fibroblast late-remodeling state，定义为D30 score减D1/D7平均score。",
  "- Key secondary: fibroblast ordinal wound-state position。",
  "- Supporting: D30-like score、wound activation、Skin-like score。",
  "- Exploratory: off-trajectory ratio与原early-state persistence。",
  "",
  "解释边界：",
  "- keloid_3L和keloid_3R属于同一患者K3，后续必须先合并到患者层面。",
  "- 本阶段的fibroblast-like为高特异性标志物定义，不等同于完整无监督细胞图谱。",
  "- 不允许根据第六阶段B结果重新调整本阶段的标志物规则。",
  "- healthy skin仅用于描述，不参与3 vs 3主要统计。"
)

safe_write_lines(
  decision_lines,
  file.path(
    REPORT_DIR,
    "06_STAGE6A_DECISION.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    REPORT_DIR,
    "07_SESSION_INFO.txt"
  )
)

safe_write_lines(
  c(
    paste0(
      "STAGE6A_COMPLETED=",
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
      "CELL_CALLING_PASS=",
      cell_calling_pass
    ),
    paste0(
      "QC_PASS=",
      qc_pass
    ),
    paste0(
      "FIBRO_HIGH_PASS=",
      fibro_high_pass
    ),
    paste0(
      "FIBRO_BROAD_PASS=",
      fibro_broad_pass
    )
  ),
  file.path(
    REPORT_DIR,
    "STAGE6A_COMPLETED.txt"
  )
)

# 删除临时解包数据，保留objects。
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

  command <- paste0(
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
        command
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

  setwd(
    STAGE_DIR
  )

  try(
    utils::tar(
      tarfile = PACKAGE_TARGZ,
      files = basename(
        package_dir
      ),
      compression = "gzip",
      tar = "internal"
    ),
    silent = TRUE
  )

  setwd(
    old_wd
  )
}

unlink(
  package_dir,
  recursive = TRUE,
  force = TRUE
)

log_msg(
  "第六阶段A完成"
)

cat("\n============================================================\n")
cat("第六阶段A运行完成。\n")
cat(
  "本地pseudobulk对象保存在：",
  file.path(
    OBJECT_DIR,
    "GSE181316_fibroblast_like_pseudobulk_list.rds"
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
