# =============================================================================
# 项目：瘢痕疙瘩“晚期重塑状态持续与伤口状态终止失败”
# 第九阶段：GSE265972静脉溃疡疾病对照与修复失败模式特异性分析
#
# 根目录：
#   C:/Users/33652/Desktop/LZZ/1、yssj
#
# 数据：
#   Normal skin: NS23, NS24, NS36, NS62E, NS63D（n=5）
#   Venous ulcer: VU1, VU4, VU5, VU6（n=4）
#
# 分析定位：
#   - 本阶段不是瘢痕疙瘩的第三/第四个验证队列；
#   - 它是慢性不愈合伤口疾病对照，用于检验不同修复失败模式；
#   - 主要分析作者注释的FB1+FB2成纤维细胞；
#   - FB1和FB2分别作为预设亚型敏感性分析；
#   - 不在GSE265972中重新选择伤口阶段签名；
#   - 不根据结果修改前面已经固定的瘢痕疙瘩结论。
#
# 主要描述指标：
#   late_remodeling_state
#   early_state_persistence
#   late_vs_early_balance
#   ordinal_wound_state_position
#   off_trajectory_ratio
#
# 解释原则：
#   - 若静脉溃疡以early-state persistence或off-trajectory为主，
#     而非稳定D30-like late-remodeling dominance，则支持不同修复失败模式；
#   - 若静脉溃疡同样表现为强late-remodeling dominance，
#     则晚期重塑持续可能是慢性伤口共有现象，不能称为瘢痕疙瘩特异机制；
#   - GSE265972正常对照为正常皮肤，不是成熟瘢痕，跨疾病效应量只比较方向，
#     不作直接数值等价解释。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANG = "en")

ROOT_DIR <- "C:/Users/33652/Desktop/LZZ/1、yssj"

PROCESSED_TAR <- file.path(
  ROOT_DIR,
  "GSE265972_processed.tar.gz"
)

METADATA_FILE <- file.path(
  ROOT_DIR,
  "GSE265972_VU_anno_metadata.txt.gz"
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

STAGE5_RDS <- file.path(
  ROOT_DIR,
  "05_STAGE5_GSE163973_MAPPING",
  "objects",
  "GSE163973_fixed_wound_state_mapping.rds"
)

STAGE6B_RDS <- file.path(
  ROOT_DIR,
  "06B_STAGE6B_GSE181316_REPLICATION",
  "objects",
  "GSE181316_fixed_wound_state_replication.rds"
)

STAGE7_RDS <- file.path(
  ROOT_DIR,
  "07_STAGE7_GSE220300_ACTIVITY_REGION",
  "objects",
  "GSE220300_activity_region_wound_state_support.rds"
)

STAGE_DIR <- file.path(
  ROOT_DIR,
  "09_STAGE9_GSE265972_REVERSE_CONTROL"
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
  "第九阶段_GSE265972静脉溃疡疾病对照检查包.zip"
)

PACKAGE_TARGZ <- file.path(
  ROOT_DIR,
  "第九阶段_GSE265972静脉溃疡疾病对照检查包.tar.gz"
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
  "00_stage9_log.txt"
)

WARNING_FILE <- file.path(
  REPORT_DIR,
  "00_stage9_warnings.txt"
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
  "Matrix",
  quietly = TRUE
)) {
  stop(
    "无法安装或加载Matrix包。"
  )
}

# ----------------------------- 固定样本设计 ------------------------------------
sample_design <- data.frame(
  sample_id = c(
    "NS23",
    "NS24",
    "NS36",
    "NS62E",
    "NS63D",
    "VU1",
    "VU4",
    "VU5",
    "VU6"
  ),
  group = c(
    rep(
      "normal_skin",
      5
    ),
    rep(
      "venous_ulcer",
      4
    )
  ),
  group_code = c(
    rep(
      0L,
      5
    ),
    rep(
      1L,
      4
    )
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

read_metadata <- function(path) {
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

  utils::read.delim(
    con,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

normalize_barcode_core <- function(x) {
  y <- toupper(
    trimws(
      as.character(
        x
      )
    )
  )

  out <- rep(
    NA_character_,
    length(
      y
    )
  )

  hit <- grepl(
    "([ACGTN]{12,})(?:[-_.][0-9]+)?$",
    y,
    perl = TRUE
  )

  if (any(
    hit
  )) {
    out[
      hit
    ] <- sub(
      "^.*?([ACGTN]{12,})(?:[-_.][0-9]+)?$",
      "\\1",
      y[
        hit
      ],
      perl = TRUE
    )
  }

  out
}

extract_zip_sample_id <- function(path) {
  toupper(
    sub(
      "\\.zip$",
      "",
      basename(
        path
      ),
      ignore.case = TRUE
    )
  )
}

extract_triplet_from_zip <- function(
  zip_path,
  output_dir
) {
  if (dir.exists(
    output_dir
  )) {
    unlink(
      output_dir,
      recursive = TRUE,
      force = TRUE
    )
  }

  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  utils::unzip(
    zip_path,
    exdir = output_dir
  )

  files <- list.files(
    output_dir,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  )

  files <- files[
    !grepl(
      "(^|[/\\\\])__MACOSX([/\\\\]|$)",
      files,
      ignore.case = TRUE,
      perl = TRUE
    ) &
      !startsWith(
        basename(
          files
        ),
        "._"
      )
  ]

  locate_one <- function(pattern) {
    candidates <- files[
      grepl(
        pattern,
        basename(
          files
        ),
        ignore.case = TRUE,
        perl = TRUE
      )
    ]

    if (length(
      candidates
    ) !=
      1L) {
      stop(
        "无法唯一定位文件：",
        pattern,
        " | ZIP=",
        basename(
          zip_path
        ),
        " | n=",
        length(
          candidates
        )
      )
    }

    candidates[1]
  }

  list(
    matrix = locate_one(
      "^matrix\\.mtx(\\.gz)?$"
    ),
    barcodes = locate_one(
      "^barcodes\\.tsv(\\.gz)?$"
    ),
    features = locate_one(
      "^(features|genes)\\.tsv(\\.gz)?$"
    )
  )
}

collapse_vector_by_symbol <- function(
  value,
  gene_symbol
) {
  gene <- toupper(
    trimws(
      as.character(
        gene_symbol
      )
    )
  )

  valid <- (
    !is.na(
      gene
    ) &
      nzchar(
        gene
      ) &
      is.finite(
        value
      )
  )

  if (!any(
    valid
  )) {
    stop(
      "没有可用于基因符号合并的有效值。"
    )
  }

  summed <- rowsum(
    matrix(
      value[
        valid
      ],
      ncol = 1
    ),
    group = gene[
      valid
    ],
    reorder = FALSE
  )

  result <- as.numeric(
    summed[, 1]
  )

  names(
    result
  ) <- rownames(
    summed
  )

  result
}

canonicalize_reference_counts <- function(
  counts,
  gene_symbol
) {
  gene <- toupper(
    trimws(
      as.character(
        gene_symbol
      )
    )
  )

  valid <- (
    !is.na(
      gene
    ) &
      nzchar(
        gene
      )
  )

  counts <- counts[
    valid,
    ,
    drop = FALSE
  ]

  gene <- gene[
    valid
  ]

  counts <- as.matrix(
    counts
  )

  collapsed <- rowsum(
    counts,
    group = gene,
    reorder = FALSE
  )

  storage.mode(
    collapsed
  ) <- "double"

  collapsed
}

merge_named_vector_list <- function(
  vector_list,
  required_names,
  minimum_common_genes = 10000L
) {
  missing_samples <- setdiff(
    required_names,
    names(
      vector_list
    )
  )

  if (length(
    missing_samples
  )) {
    stop(
      "pseudobulk列表缺少样本：",
      paste(
        missing_samples,
        collapse = ", "
      )
    )
  }

  vector_list <- vector_list[
    required_names
  ]

  common_genes <- Reduce(
    intersect,
    lapply(
      vector_list,
      names
    )
  )

  if (length(
    common_genes
  ) <
      minimum_common_genes
  ) {
    warn_msg(
      "样本间共同基因数低于",
      minimum_common_genes,
      "：",
      length(
        common_genes
      )
    )
  }

  result <- do.call(
    cbind,
    lapply(
      vector_list,
      function(x) {
        x[
          common_genes
        ]
      }
    )
  )

  rownames(
    result
  ) <- common_genes

  colnames(
    result
  ) <- names(
    vector_list
  )

  storage.mode(
    result
  ) <- "double"

  result
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
    colnames(
      counts
    ),
    sample_info$key
  )

  if (anyNA(
    idx
  )) {
    stop(
      "正常伤口pseudobulk列名无法匹配sample_info。"
    )
  }

  result <- sample_info[
    idx,
    ,
    drop = FALSE
  ]

  if (!identical(
    colnames(
      counts
    ),
    result$key
  )) {
    stop(
      "正常伤口pseudobulk顺序匹配失败。"
    )
  }

  result
}

get_stable_fibroblast_signatures <- function(
  model
) {
  stability <- model$signature_stability

  required_columns <- c(
    "cell_type",
    "stage",
    "gene",
    "stable_all_3_folds"
  )

  if (!all(
    required_columns %in%
      names(
        stability
      )
  )) {
    stop(
      "第四阶段模型缺少稳定签名字段。"
    )
  }

  x <- stability[
    stability$cell_type ==
      "Fibroblast" &
      stability$stable_all_3_folds,
    ,
    drop = FALSE
  ]

  stages <- c(
    "Skin",
    "Wound1",
    "Wound7",
    "Wound30"
  )

  signatures <- setNames(
    lapply(
      stages,
      function(stage) {
        unique(
          toupper(
            trimws(
              as.character(
                x$gene[
                  x$stage ==
                    stage
                ]
              )
            )
          )
        )
      }
    ),
    stages
  )

  n_genes <- vapply(
    signatures,
    length,
    integer(1)
  )

  if (any(
    n_genes <
      20L
  )) {
    stop(
      "稳定Fibroblast签名不足：",
      paste(
        names(
          n_genes
        ),
        n_genes,
        sep = "=",
        collapse = ", "
      )
    )
  }

  signatures
}

exclude_background_gene <- function(gene) {
  g <- toupper(
    as.character(
      gene
    )
  )

  grepl(
    "^MT-|^RPS[0-9]|^RPL[0-9]|^HB[ABDEGMQZ][0-9A-Z]*$",
    g
  )
}

log_cpm <- function(counts) {
  library_size <- colSums(
    counts
  )

  if (any(
    library_size <=
      0
  )) {
    stop(
      "存在文库总计数≤0。"
    )
  }

  log2(
    sweep(
      counts +
        0.5,
      2,
      library_size +
        1,
      "/"
    ) *
      1e6
  )
}

score_samples_by_rank <- function(
  expression_matrix,
  signatures
) {
  stages <- c(
    "Skin",
    "Wound1",
    "Wound7",
    "Wound30"
  )

  result <- matrix(
    NA_real_,
    nrow = ncol(
      expression_matrix
    ),
    ncol = length(
      stages
    ),
    dimnames = list(
      colnames(
        expression_matrix
      ),
      stages
    )
  )

  for (
    j in seq_len(
      ncol(
        expression_matrix
      )
    )
  ) {
    ranks <- rank(
      expression_matrix[, j],
      ties.method = "average",
      na.last = "keep"
    )

    ranks <- ranks /
      sum(
        is.finite(
          ranks
        )
      )

    names(
      ranks
    ) <- rownames(
      expression_matrix
    )

    for (
      stage in stages
    ) {
      genes <- intersect(
        signatures[[stage]],
        names(
          ranks
        )
      )

      if (length(
        genes
      ) <
          15L) {
        stop(
          "阶段",
          stage,
          "共同稳定签名少于15：",
          length(
            genes
          )
        )
      }

      result[j, stage] <- mean(
        ranks[
          genes
        ],
        na.rm = TRUE
      )
    }
  }

  result
}

standardize_by_reference <- function(
  reference_scores,
  new_scores
) {
  center <- colMeans(
    reference_scores
  )

  scale_value <- apply(
    reference_scores,
    2,
    stats::sd
  )

  scale_value[
    !is.finite(
      scale_value
    ) |
      scale_value ==
      0
  ] <- 1

  reference_z <- sweep(
    sweep(
      reference_scores,
      2,
      center,
      "-"
    ),
    2,
    scale_value,
    "/"
  )

  new_z <- sweep(
    sweep(
      new_scores,
      2,
      center,
      "-"
    ),
    2,
    scale_value,
    "/"
  )

  list(
    reference_z =
      reference_z,
    new_z =
      new_z,
    center =
      center,
    scale =
      scale_value
  )
}

calculate_reference_centroids <- function(
  reference_z,
  reference_design
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
      function(stage) {
        colMeans(
          reference_z[
            reference_design$condition ==
              stage,
            ,
            drop = FALSE
          ]
        )
      }
    )
  )

  rownames(
    centroids
  ) <- stages

  centroids
}

calculate_reference_lodo_distance <- function(
  reference_z,
  reference_design
) {
  output <- numeric(
    nrow(
      reference_z
    )
  )

  for (
    i in seq_len(
      nrow(
        reference_z
      )
    )
  ) {
    stage <- reference_design$condition[i]
    donor <- reference_design$donor[i]

    train_idx <- which(
      reference_design$condition ==
        stage &
        reference_design$donor !=
        donor
    )

    centroid <- colMeans(
      reference_z[
        train_idx,
        ,
        drop = FALSE
      ]
    )

    output[i] <- sqrt(
      sum(
        (
          reference_z[i, ] -
            centroid
        ) ^
            2
      )
    )
  }

  output
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
    nrow(
      z_scores
    )
  )

  for (
    i in seq_len(
      nrow(
        z_scores
      )
    )
  ) {
    distance <- apply(
      centroids,
      1,
      function(centroid) {
        sqrt(
          sum(
            (
              z_scores[i, ] -
                centroid
            ) ^
              2
          )
        )
      }
    )

    nearest_stage <- names(
      which.min(
        distance
      )
    )

    weights <- exp(
      -0.5 *
        distance ^
        2
    )

    if (
      !all(
        is.finite(
          weights
        )
      ) ||
        sum(
          weights
        ) ==
        0
    ) {
      weights <- rep(
        1 /
          length(
            distance
          ),
        length(
          distance
        )
      )
    } else {
      weights <- weights /
        sum(
          weights
        )
    }

    names(
      weights
    ) <- names(
      distance
    )

    entropy <- -sum(
      weights *
        log(
          pmax(
            weights,
            1e-12
          )
        )
    ) /
      log(
        length(
          weights
        )
      )

    rows[[i]] <- data.frame(
      sample_key = rownames(
        z_scores
      )[i],
      nearest_stage =
        nearest_stage,
      ordinal_wound_state_position =
        sum(
          weights *
            stage_position[
              names(
                weights
              )
            ]
        ),
      off_trajectory_distance =
        min(
          distance
        ),
      off_trajectory_ratio =
        min(
          distance
        ) /
        reference_distance_limit,
      mixed_state_entropy =
        entropy,
      weight_Skin =
        weights[
          "Skin"
        ],
      weight_Wound1 =
        weights[
          "Wound1"
        ],
      weight_Wound7 =
        weights[
          "Wound7"
        ],
      weight_Wound30 =
        weights[
          "Wound30"
        ],
      stringsAsFactors = FALSE
    )
  }

  do.call(
    rbind,
    rows
  )
}

add_fixed_contrasts <- function(
  mapping,
  z_scores
) {
  idx <- match(
    mapping$sample_key,
    rownames(
      z_scores
    )
  )

  if (anyNA(
    idx
  )) {
    stop(
      "映射结果无法匹配z-score。"
    )
  }

  z <- z_scores[
    idx,
    ,
    drop = FALSE
  ]

  mapping$z_Skin <-
    z[, "Skin"]

  mapping$z_Wound1 <-
    z[, "Wound1"]

  mapping$z_Wound7 <-
    z[, "Wound7"]

  mapping$z_Wound30 <-
    z[, "Wound30"]

  mapping$late_remodeling_state <-
    z[, "Wound30"] -
    (
      z[, "Wound1"] +
        z[, "Wound7"]
    ) /
    2

  mapping$early_state_persistence <-
    (
      z[, "Wound1"] +
        z[, "Wound7"]
    ) /
    2 -
    (
      z[, "Skin"] +
        z[, "Wound30"]
    ) /
    2

  mapping$late_vs_early_balance <-
    mapping$late_remodeling_state -
    mapping$early_state_persistence

  mapping$wound_activation <-
    (
      z[, "Wound1"] +
        z[, "Wound7"] +
        z[, "Wound30"]
    ) /
    3 -
    z[, "Skin"]

  mapping$skin_return_state <-
    z[, "Skin"] -
    (
      z[, "Wound1"] +
        z[, "Wound7"] +
        z[, "Wound30"]
    ) /
    3

  mapping
}

exact_group_comparison <- function(
  value,
  group,
  group_a = "venous_ulcer",
  group_b = "normal_skin"
) {
  ok <- (
    is.finite(
      value
    ) &
      group %in%
      c(
        group_a,
        group_b
      )
  )

  value <- value[
    ok
  ]

  group <- group[
    ok
  ]

  x <- value[
    group ==
      group_a
  ]

  y <- value[
    group ==
      group_b
  ]

  n1 <- length(
    x
  )

  n2 <- length(
    y
  )

  if (
    n1 <
      2L ||
      n2 <
      2L
  ) {
    return(
      list(
        n_a = n1,
        n_b = n2,
        mean_a = NA_real_,
        mean_b = NA_real_,
        median_a = NA_real_,
        median_b = NA_real_,
        mean_difference = NA_real_,
        hodges_lehmann_shift = NA_real_,
        hedges_g = NA_real_,
        cliffs_delta = NA_real_,
        exact_two_sided_p = NA_real_,
        exact_one_sided_greater_p = NA_real_,
        complete_separation_greater = NA,
        complete_separation_less = NA
      )
    )
  }

  observed_difference <- mean(
    x
  ) -
    mean(
      y
    )

  all_values <- c(
    x,
    y
  )

  combinations <- utils::combn(
    seq_along(
      all_values
    ),
    n1
  )

  permutation_difference <- apply(
    combinations,
    2,
    function(group_a_indices) {
      mean(
        all_values[
          group_a_indices
        ]
      ) -
        mean(
          all_values[
            -group_a_indices
        ]
      )
    }
  )

  exact_two_sided_p <- mean(
    abs(
      permutation_difference
    ) >=
      abs(
        observed_difference
      ) -
      1e-12
  )

  exact_one_sided_greater_p <- mean(
    permutation_difference >=
      observed_difference -
      1e-12
  )

  degrees_freedom <- n1 +
    n2 -
    2

  pooled_sd <- sqrt(
    (
      (n1 -
        1) *
        stats::var(
          x
        ) +
        (n2 -
          1) *
        stats::var(
          y
        )
    ) /
      degrees_freedom
  )

  hedges_g <- NA_real_

  if (
    is.finite(
      pooled_sd
    ) &&
      pooled_sd >
      0
  ) {
    cohen_d <-
      observed_difference /
      pooled_sd

    correction <-
      1 -
      3 /
      (
        4 *
          degrees_freedom -
          1
      )

    hedges_g <-
      correction *
      cohen_d
  }

  pairwise_difference <- outer(
    x,
    y,
    "-"
  )

  cliffs_delta <- mean(
    pairwise_difference >
      0
  ) -
    mean(
      pairwise_difference <
        0
    )

  list(
    n_a =
      n1,
    n_b =
      n2,
    mean_a =
      mean(
        x
      ),
    mean_b =
      mean(
        y
      ),
    median_a =
      stats::median(
        x
      ),
    median_b =
      stats::median(
        y
      ),
    mean_difference =
      observed_difference,
    hodges_lehmann_shift =
      stats::median(
        as.vector(
          pairwise_difference
        )
      ),
    hedges_g =
      hedges_g,
    cliffs_delta =
      cliffs_delta,
    exact_two_sided_p =
      exact_two_sided_p,
    exact_one_sided_greater_p =
      exact_one_sided_greater_p,
    complete_separation_greater =
      min(
        x
      ) >
      max(
        y
      ),
    complete_separation_less =
      max(
        x
      ) <
      min(
        y
      )
  )
}

effect_strength <- function(
  hedges_g,
  cliffs_delta
) {
  if (
    !is.finite(
      hedges_g
    ) &&
      !is.finite(
        cliffs_delta
      )
  ) {
    return(
      "not_estimable"
    )
  }

  if (
    (
      is.finite(
        hedges_g
      ) &&
        abs(
          hedges_g
        ) >=
        0.8
    ) ||
      (
        is.finite(
          cliffs_delta
        ) &&
          abs(
            cliffs_delta
          ) >=
          0.56
      )
  ) {
    return(
      "large"
    )
  }

  if (
    (
      is.finite(
        hedges_g
      ) &&
        abs(
          hedges_g
        ) >=
        0.5
    ) ||
      (
        is.finite(
          cliffs_delta
        ) &&
          abs(
            cliffs_delta
          ) >=
          0.33
      )
  ) {
    return(
      "moderate"
    )
  }

  "small"
}

analyse_fibroblast_set <- function(
  set_name,
  vector_list,
  required_samples,
  sample_design,
  reference_counts,
  reference_design,
  signatures
) {
  log_msg(
    "固定伤口状态映射：",
    set_name
  )

  counts <- merge_named_vector_list(
    vector_list,
    required_names =
      required_samples
  )

  common_genes <- intersect(
    rownames(
      reference_counts
    ),
    rownames(
      counts
    )
  )

  common_genes <- common_genes[
    !exclude_background_gene(
      common_genes
    )
  ]

  if (
    length(
      common_genes
    ) <
      5000L
  ) {
    stop(
      set_name,
      "与正常伤口参考共同背景基因少于5000：",
      length(
        common_genes
      )
    )
  }

  signature_coverage <- vapply(
    signatures,
    function(genes) {
      sum(
        genes %in%
          common_genes
      )
    },
    integer(1)
  )

  if (any(
    signature_coverage <
      15L
  )) {
    stop(
      set_name,
      "至少一个阶段共同稳定签名少于15。"
    )
  }

  reference_expression <- log_cpm(
    reference_counts[
      common_genes,
      ,
      drop = FALSE
    ]
  )

  new_expression <- log_cpm(
    counts[
      common_genes,
      ,
      drop = FALSE
    ]
  )

  reference_scores <- score_samples_by_rank(
    reference_expression,
    signatures
  )

  new_scores <- score_samples_by_rank(
    new_expression,
    signatures
  )

  standardized <- standardize_by_reference(
    reference_scores,
    new_scores
  )

  reference_z <-
    standardized$reference_z

  new_z <-
    standardized$new_z

  centroids <- calculate_reference_centroids(
    reference_z,
    reference_design
  )

  reference_lodo_distance <-
    calculate_reference_lodo_distance(
      reference_z,
      reference_design
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
      reference_distance_limit <=
      0
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
      reference_distance_limit <=
      0
  ) {
    reference_distance_limit <- 1
  }

  mapping <- map_to_reference(
    new_z,
    centroids,
    reference_distance_limit
  )

  mapping <- add_fixed_contrasts(
    mapping,
    new_z
  )

  mapping$analysis_set <-
    set_name

  mapping <- merge(
    mapping,
    sample_design[
      ,
      c(
        "sample_id",
        "group",
        "group_code"
      )
    ],
    by.x = "sample_key",
    by.y = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )

  mapping <- mapping[
    match(
      rownames(
        new_z
      ),
      mapping$sample_key
    ),
    ,
    drop = FALSE
  ]

  metrics <- c(
    "late_remodeling_state",
    "early_state_persistence",
    "late_vs_early_balance",
    "ordinal_wound_state_position",
    "off_trajectory_ratio",
    "wound_activation",
    "z_Skin",
    "z_Wound1",
    "z_Wound7",
    "z_Wound30",
    "skin_return_state",
    "mixed_state_entropy"
  )

  comparison_rows <- list()

  for (
    metric in metrics
  ) {
    test <- exact_group_comparison(
      mapping[[metric]],
      mapping$group
    )

    comparison_rows[[length(
      comparison_rows
    ) + 1L]] <- data.frame(
      analysis_set =
        set_name,
      metric =
        metric,
      n_venous_ulcer =
        test$n_a,
      n_normal_skin =
        test$n_b,
      mean_venous_ulcer =
        test$mean_a,
      mean_normal_skin =
        test$mean_b,
      median_venous_ulcer =
        test$median_a,
      median_normal_skin =
        test$median_b,
      mean_difference =
        test$mean_difference,
      hodges_lehmann_shift =
        test$hodges_lehmann_shift,
      hedges_g =
        test$hedges_g,
      cliffs_delta =
        test$cliffs_delta,
      exact_two_sided_p =
        test$exact_two_sided_p,
      exact_one_sided_greater_p =
        test$exact_one_sided_greater_p,
      complete_separation_greater =
        test$complete_separation_greater,
      complete_separation_less =
        test$complete_separation_less,
      effect_strength =
        effect_strength(
          test$hedges_g,
          test$cliffs_delta
        ),
      stringsAsFactors = FALSE
    )
  }

  comparisons <- do.call(
    rbind,
    comparison_rows
  )

  coverage <- data.frame(
    analysis_set =
      set_name,
    stage =
      names(
        signature_coverage
      ),
    n_stable_signature_genes =
      vapply(
        signatures,
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
        signatures,
        length,
        integer(1)
      ),
    n_common_background_genes =
      length(
        common_genes
      ),
    reference_distance_limit =
      reference_distance_limit,
    stringsAsFactors = FALSE
  )

  list(
    counts =
      counts,
    mapping =
      mapping,
    comparisons =
      comparisons,
    coverage =
      coverage,
    calibration = list(
      common_genes =
        common_genes,
      signatures =
        signatures,
      reference_center =
        standardized$center,
      reference_scale =
        standardized$scale,
      centroids =
        centroids,
      reference_lodo_distance =
        reference_lodo_distance,
      reference_distance_limit =
        reference_distance_limit
    )
  )
}

# ----------------------------- 输入检查 ----------------------------------------
log_msg(
  "第九阶段开始"
)

required_files <- c(
  PROCESSED_TAR,
  METADATA_FILE,
  REFERENCE_PB_RDS,
  REFERENCE_MODEL_RDS
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (length(
  missing_files
)) {
  stop(
    "缺少必要文件：",
    paste(
      missing_files,
      collapse = " | "
    )
  )
}

metadata <- read_metadata(
  METADATA_FILE
)

required_metadata_columns <- c(
  "samples",
  "Condition",
  "barcode",
  "CellType"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  names(
    metadata
  )
)

if (length(
  missing_metadata_columns
)) {
  stop(
    "GSE265972 metadata缺少字段：",
    paste(
      missing_metadata_columns,
      collapse = ", "
    )
  )
}

metadata$sample_id <- toupper(
  trimws(
    as.character(
      metadata$samples
    )
  )
)

metadata$group <- ifelse(
  toupper(
    trimws(
      as.character(
        metadata$Condition
      )
    )
  ) ==
    "VU",
  "venous_ulcer",
  "normal_skin"
)

metadata$barcode_original <-
  as.character(
    metadata$barcode
  )

metadata$barcode_core <-
  normalize_barcode_core(
    metadata$barcode_original
  )

metadata$CellType <-
  toupper(
    trimws(
      as.character(
        metadata$CellType
      )
    )
  )

metadata$ct <-
  if ("ct" %in%
      names(
        metadata
      )) {
    toupper(
      trimws(
        as.character(
          metadata$ct
        )
      )
    )
  } else {
    NA_character_
  }

metadata_sample_counts <- as.data.frame(
  table(
    metadata$sample_id,
    metadata$group,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

names(
  metadata_sample_counts
) <- c(
  "sample_id",
  "group",
  "n_metadata_cells"
)

metadata_sample_counts <-
  metadata_sample_counts[
    metadata_sample_counts$n_metadata_cells >
      0,
    ,
    drop = FALSE
  ]

safe_write_csv(
  metadata_sample_counts,
  file.path(
    REPORT_DIR,
    "02_metadata_sample_counts.csv"
  )
)

metadata_celltype_counts <- as.data.frame(
  table(
    metadata$sample_id,
    metadata$group,
    metadata$CellType,
    useNA = "ifany"
  ),
  stringsAsFactors = FALSE
)

names(
  metadata_celltype_counts
) <- c(
  "sample_id",
  "group",
  "cell_type",
  "n_cells"
)

metadata_celltype_counts <-
  metadata_celltype_counts[
    metadata_celltype_counts$n_cells >
      0,
    ,
    drop = FALSE
  ]

safe_write_csv(
  metadata_celltype_counts,
  file.path(
    REPORT_DIR,
    "03_author_celltype_counts.csv"
  )
)

missing_metadata_samples <- setdiff(
  sample_design$sample_id,
  unique(
    metadata$sample_id
  )
)

if (length(
  missing_metadata_samples
)) {
  stop(
    "metadata缺少样本：",
    paste(
      missing_metadata_samples,
      collapse = ", "
    )
  )
}

# ----------------------------- 正常参考 ----------------------------------------
reference_pb_object <- readRDS(
  REFERENCE_PB_RDS
)

reference_model <- readRDS(
  REFERENCE_MODEL_RDS
)

reference_design_all <- get_reference_design(
  reference_pb_object$counts,
  reference_pb_object$sample_info
)

reference_counts_all <-
  canonicalize_reference_counts(
    reference_pb_object$counts,
    reference_pb_object$gene_symbol
  )

reference_fibroblast_idx <- which(
  reference_design_all$main_cell_type ==
    "Fibroblast"
)

if (
  length(
    reference_fibroblast_idx
  ) !=
    12L
) {
  stop(
    "正常伤口Fibroblast pseudobulk不是12个。"
  )
}

reference_counts <- reference_counts_all[
  ,
  reference_fibroblast_idx,
  drop = FALSE
]

reference_design <- reference_design_all[
  reference_fibroblast_idx,
  ,
  drop = FALSE
]

signatures <- get_stable_fibroblast_signatures(
  reference_model
)

signature_summary <- data.frame(
  stage =
    names(
      signatures
    ),
  n_stable_genes =
    vapply(
      signatures,
      length,
      integer(1)
    ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  signature_summary,
  file.path(
    REPORT_DIR,
    "04_fixed_signature_summary.csv"
  )
)

# ----------------------------- 解包九个样本ZIP ---------------------------------
outer_dir <- file.path(
  WORK_DIR,
  "GSE265972_outer"
)

dir.create(
  outer_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

log_msg(
  "解包GSE265972_processed.tar.gz"
)

utils::untar(
  PROCESSED_TAR,
  exdir = outer_dir
)

nested_zips <- list.files(
  outer_dir,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.zip$",
  ignore.case = TRUE
)

zip_sample_ids <- vapply(
  nested_zips,
  extract_zip_sample_id,
  character(1)
)

if (
  length(
    nested_zips
  ) !=
    9L ||
    !setequal(
      zip_sample_ids,
      sample_design$sample_id
    )
) {
  stop(
    "GSE265972嵌套ZIP结构与预期不一致。实际：",
    paste(
      zip_sample_ids,
      collapse = ", "
    )
  )
}

nested_zips <- nested_zips[
  match(
    sample_design$sample_id,
    zip_sample_ids
  )
]

# ----------------------------- 逐样本匹配与pseudobulk --------------------------
combined_pb_list <- list()
fb1_pb_list <- list()
fb2_pb_list <- list()
sample_qc_rows <- list()
barcode_rows <- list()
fibro_rows <- list()
feature_rows <- list()

gene_reference_id <- NULL
gene_reference_symbol <- NULL

pdf(
  file.path(
    FIG_DIR,
    "Figure_GSE265972_QC_distributions.pdf"
  ),
  width = 10,
  height = 8,
  onefile = TRUE
)

for (
  zip_path in nested_zips
) {
  sample_id <- extract_zip_sample_id(
    zip_path
  )

  log_msg(
    "处理样本：",
    sample_id
  )

  sample_tmp <- file.path(
    WORK_DIR,
    paste0(
      "sample_",
      sample_id
    )
  )

  triplet <- extract_triplet_from_zip(
    zip_path,
    sample_tmp
  )

  barcodes <- read_noheader_table(
    triplet$barcodes,
    sep = "\t"
  )

  features <- read_noheader_table(
    triplet$features,
    sep = "\t"
  )

  counts <- read_mtx_gz(
    triplet$matrix
  )

  raw_barcode_original <- as.character(
    barcodes[[1]]
  )

  raw_barcode_core <- normalize_barcode_core(
    raw_barcode_original
  )

  if (
    anyNA(
      raw_barcode_core
    ) ||
      anyDuplicated(
        raw_barcode_core
      )
  ) {
    stop(
      sample_id,
      "原始barcode标准化失败或重复：NA=",
      sum(
        is.na(
          raw_barcode_core
        )
      ),
      "; duplicated=",
      anyDuplicated(
        raw_barcode_core
      )
    )
  }

  gene_id <- as.character(
    features[[1]]
  )

  gene_symbol <- if (
    ncol(
      features
    ) >=
      2L
  ) {
    as.character(
      features[[2]]
    )
  } else {
    gene_id
  }

  if (
    nrow(
      counts
    ) !=
      length(
        gene_id
      ) ||
      ncol(
        counts
      ) !=
      length(
        raw_barcode_core
      )
  ) {
    stop(
      sample_id,
      "矩阵维度与features/barcodes不一致。"
    )
  }

  if (is.null(
    gene_reference_id
  )) {
    gene_reference_id <-
      gene_id

    gene_reference_symbol <-
      gene_symbol

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
        anyDuplicated(
          gene_id
        ) ||
          anyDuplicated(
            gene_reference_id
          )
      ) {
        stop(
          sample_id,
          "gene_id不一致且存在重复，无法安全对齐。"
        )
      }

      reorder_idx <- match(
        gene_reference_id,
        gene_id
      )

      if (anyNA(
        reorder_idx
      )) {
        stop(
          sample_id,
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

  rownames(
    counts
  ) <- make.unique(
    gene_reference_symbol
  )

  colnames(
    counts
  ) <- raw_barcode_core

  metadata_sample <- metadata[
    metadata$sample_id ==
      sample_id &
      !is.na(
        metadata$barcode_core
      ),
    ,
    drop = FALSE
  ]

  metadata_sample <- metadata_sample[
    !duplicated(
      metadata_sample$barcode_core
    ),
    ,
    drop = FALSE
  ]

  if (!nrow(
    metadata_sample
  )) {
    stop(
      sample_id,
      "没有可用metadata细胞。"
    )
  }

  match_index <- match(
    metadata_sample$barcode_core,
    raw_barcode_core
  )

  matched <- !is.na(
    match_index
  )

  match_rate <- mean(
    matched
  )

  barcode_rows[[length(
    barcode_rows
  ) + 1L]] <- data.frame(
    sample_id =
      sample_id,
    n_matrix_barcodes =
      length(
        raw_barcode_core
      ),
    n_metadata_rows =
      nrow(
        metadata_sample
      ),
    n_metadata_matched =
      sum(
        matched
      ),
    metadata_match_rate =
      match_rate,
    stringsAsFactors = FALSE
  )

  if (
    match_rate <
      0.95
  ) {
    warn_msg(
      sample_id,
      "metadata与矩阵barcode匹配率低于95%：",
      round(
        match_rate,
        4
      )
    )
  }

  if (!any(
    matched
  )) {
    stop(
      sample_id,
      "没有任何metadata barcode匹配表达矩阵。"
    )
  }

  metadata_sample <- metadata_sample[
    matched,
    ,
    drop = FALSE
  ]

  counts_matched <- counts[
    ,
    match_index[
      matched
    ],
    drop = FALSE
  ]

  colnames(
    counts_matched
  ) <- metadata_sample$barcode_core

  n_count <- Matrix::colSums(
    counts_matched
  )

  n_feature <- Matrix::colSums(
    counts_matched >
      0
  )

  gene_upper <- toupper(
    rownames(
      counts_matched
    )
  )

  mt_idx <- grepl(
    "^MT-",
    gene_upper
  )

  pct_mt <- if (
    any(
      mt_idx
    )
  ) {
    Matrix::colSums(
      counts_matched[
        mt_idx,
        ,
        drop = FALSE
      ]
    ) /
      pmax(
        n_count,
        1
      ) *
      100
  } else {
    rep(
      0,
      ncol(
        counts_matched
      )
    )
  }

  # metadata已经是作者质控和双细胞处理后的细胞；
  # 下列只作审计，不再以额外阈值删除作者保留细胞。
  audit_qc_pass <- (
    is.finite(
      n_count
    ) &
      is.finite(
        n_feature
      ) &
      is.finite(
        pct_mt
      ) &
      n_count >=
      500 &
      n_feature >=
      200 &
      pct_mt <=
      20
  )

  audit_qc_pass[
    is.na(
      audit_qc_pass
    )
  ] <- FALSE

  combined_idx <- which(
    metadata_sample$CellType %in%
      c(
        "FB1",
        "FB2"
      )
  )

  fb1_idx <- which(
    metadata_sample$CellType ==
      "FB1"
  )

  fb2_idx <- which(
    metadata_sample$CellType ==
      "FB2"
  )

  if (
    length(
      combined_idx
    ) <
      100L
  ) {
    warn_msg(
      sample_id,
      "作者注释FB1+FB2少于100：",
      length(
        combined_idx
      )
    )
  }

  if (
    length(
      combined_idx
    ) <
      50L
  ) {
    stop(
      sample_id,
      "作者注释FB1+FB2少于50，无法建立稳定pseudobulk。"
    )
  }

  if (
    length(
      fb1_idx
    ) <
      50L
  ) {
    warn_msg(
      sample_id,
      "FB1少于50：",
      length(
        fb1_idx
      )
    )
  }

  if (
    length(
      fb2_idx
    ) <
      50L
  ) {
    warn_msg(
      sample_id,
      "FB2少于50：",
      length(
        fb2_idx
      )
    )
  }

  combined_raw <- Matrix::rowSums(
    counts_matched[
      ,
      combined_idx,
      drop = FALSE
    ]
  )

  fb1_raw <- Matrix::rowSums(
    counts_matched[
      ,
      fb1_idx,
      drop = FALSE
    ]
  )

  fb2_raw <- Matrix::rowSums(
    counts_matched[
      ,
      fb2_idx,
      drop = FALSE
    ]
  )

  combined_pb_list[[sample_id]] <-
    collapse_vector_by_symbol(
      combined_raw,
      gene_reference_symbol
    )

  fb1_pb_list[[sample_id]] <-
    collapse_vector_by_symbol(
      fb1_raw,
      gene_reference_symbol
    )

  fb2_pb_list[[sample_id]] <-
    collapse_vector_by_symbol(
      fb2_raw,
      gene_reference_symbol
    )

  design_row <- sample_design[
    sample_design$sample_id ==
      sample_id,
    ,
    drop = FALSE
  ]

  sample_qc_rows[[length(
    sample_qc_rows
  ) + 1L]] <- data.frame(
    sample_id =
      sample_id,
    group =
      design_row$group,
    n_matrix_barcodes =
      ncol(
        counts
      ),
    n_author_metadata_cells =
      ncol(
        counts_matched
      ),
    author_cells_audit_QC_pass =
      sum(
        audit_qc_pass
      ),
    audit_QC_pass_fraction =
      mean(
        audit_qc_pass
      ),
    median_nCount =
      stats::median(
        n_count
      ),
    median_nFeature =
      stats::median(
        n_feature
      ),
    median_percent_mt =
      stats::median(
        pct_mt
      ),
    stringsAsFactors = FALSE
  )

  fibro_rows[[length(
    fibro_rows
  ) + 1L]] <- data.frame(
    sample_id =
      sample_id,
    group =
      design_row$group,
    n_FB1 =
      length(
        fb1_idx
      ),
    n_FB2 =
      length(
        fb2_idx
      ),
    n_FB1_FB2_combined =
      length(
        combined_idx
      ),
    combined_library_size =
      sum(
        combined_raw
      ),
    combined_detected_genes =
      sum(
        combined_raw >
          0
      ),
    stringsAsFactors = FALSE
  )

  feature_rows[[length(
    feature_rows
  ) + 1L]] <- data.frame(
    sample_id =
      sample_id,
    n_features =
      length(
        gene_id
      ),
    identical_gene_ids_to_first =
      identical_id,
    identical_gene_symbols_to_first =
      identical_symbol,
    stringsAsFactors = FALSE
  )

  op <- par(
    no.readonly = TRUE
  )

  par(
    mfrow = c(
      2,
      2
    ),
    mar = c(
      4,
      4,
      3,
      1
    )
  )

  hist(
    log10(
      n_count +
        1
    ),
    breaks = 50,
    xlab =
      "log10 nCount + 1",
    main = paste0(
      sample_id,
      " counts"
    )
  )

  hist(
    n_feature,
    breaks = 50,
    xlab =
      "Detected genes",
    main = paste0(
      sample_id,
      " features"
    )
  )

  hist(
    pct_mt,
    breaks = 50,
    xlab =
      "Mitochondrial %",
    main = paste0(
      sample_id,
      " mitochondrial %"
    )
  )

  barplot(
    c(
      FB1 =
        length(
          fb1_idx
        ),
      FB2 =
        length(
          fb2_idx
        ),
      Other =
        ncol(
          counts_matched
        ) -
        length(
          combined_idx
        )
    ),
    ylab =
      "Author-annotated cells",
    main = paste0(
      sample_id,
      " cell composition"
    )
  )

  par(
    op
  )

  rm(
    counts,
    counts_matched
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

sample_qc_summary <- do.call(
  rbind,
  sample_qc_rows
)

barcode_summary <- do.call(
  rbind,
  barcode_rows
)

fibro_summary <- do.call(
  rbind,
  fibro_rows
)

feature_summary <- do.call(
  rbind,
  feature_rows
)

rownames(
  sample_qc_summary
) <- NULL

rownames(
  barcode_summary
) <- NULL

rownames(
  fibro_summary
) <- NULL

rownames(
  feature_summary
) <- NULL

safe_write_csv(
  sample_qc_summary,
  file.path(
    REPORT_DIR,
    "05_sample_QC_summary.csv"
  )
)

safe_write_csv(
  barcode_summary,
  file.path(
    REPORT_DIR,
    "06_barcode_matching_summary.csv"
  )
)

safe_write_csv(
  fibro_summary,
  file.path(
    REPORT_DIR,
    "07_fibroblast_cell_counts.csv"
  )
)

safe_write_csv(
  feature_summary,
  file.path(
    REPORT_DIR,
    "08_feature_consistency.csv"
  )
)

# ----------------------------- 固定映射 ----------------------------------------
combined_result <- analyse_fibroblast_set(
  set_name =
    "FB1_FB2_combined_primary",
  vector_list =
    combined_pb_list,
  required_samples =
    sample_design$sample_id,
  sample_design =
    sample_design,
  reference_counts =
    reference_counts,
  reference_design =
    reference_design,
  signatures =
    signatures
)

fb1_result <- analyse_fibroblast_set(
  set_name =
    "FB1_sensitivity",
  vector_list =
    fb1_pb_list,
  required_samples =
    sample_design$sample_id,
  sample_design =
    sample_design,
  reference_counts =
    reference_counts,
  reference_design =
    reference_design,
  signatures =
    signatures
)

fb2_result <- analyse_fibroblast_set(
  set_name =
    "FB2_sensitivity",
  vector_list =
    fb2_pb_list,
  required_samples =
    sample_design$sample_id,
  sample_design =
    sample_design,
  reference_counts =
    reference_counts,
  reference_design =
    reference_design,
  signatures =
    signatures
)

sample_scores <- rbind(
  combined_result$mapping,
  fb1_result$mapping,
  fb2_result$mapping
)

comparisons <- rbind(
  combined_result$comparisons,
  fb1_result$comparisons,
  fb2_result$comparisons
)

coverage <- rbind(
  combined_result$coverage,
  fb1_result$coverage,
  fb2_result$coverage
)

rownames(
  sample_scores
) <- NULL

rownames(
  comparisons
) <- NULL

rownames(
  coverage
) <- NULL

comparisons$BH_FDR_all_tests <-
  stats::p.adjust(
    comparisons$exact_two_sided_p,
    method = "BH"
  )

safe_write_csv(
  coverage,
  file.path(
    REPORT_DIR,
    "09_signature_gene_coverage.csv"
  )
)

safe_write_csv(
  sample_scores,
  file.path(
    REPORT_DIR,
    "10_GSE265972_wound_state_scores.csv"
  )
)

safe_write_csv(
  comparisons,
  file.path(
    REPORT_DIR,
    "11_GSE265972_group_comparisons.csv"
  )
)

# ----------------------------- 跨疾病方向汇总 ----------------------------------
cross_disease_rows <- list()

primary_metrics <- c(
  "late_remodeling_state",
  "early_state_persistence",
  "ordinal_wound_state_position",
  "off_trajectory_ratio",
  "wound_activation",
  "z_Skin",
  "z_Wound30"
)

vu_primary <- comparisons[
  comparisons$analysis_set ==
    "FB1_FB2_combined_primary" &
    comparisons$metric %in%
    primary_metrics,
  ,
  drop = FALSE
]

cross_disease_rows[[
  "GSE265972"
]] <- data.frame(
  cohort =
    "GSE265972",
  disease =
    "venous_ulcer",
  analysis_set =
    "FB1_FB2_combined_primary",
  metric =
    vu_primary$metric,
  mean_difference =
    vu_primary$mean_difference,
  hedges_g =
    vu_primary$hedges_g,
  cliffs_delta =
    vu_primary$cliffs_delta,
  stringsAsFactors = FALSE
)

if (file.exists(
  STAGE5_RDS
)) {
  stage5 <- readRDS(
    STAGE5_RDS
  )

  if ("comparisons" %in%
      names(
        stage5
      )) {
    x <- stage5$comparisons

    metric_map <- data.frame(
      original = c(
        "remodeling_completion",
        "early_state_persistence",
        "ordinal_position",
        "off_trajectory_ratio",
        "wound_activation",
        "z_Skin",
        "z_Wound30"
      ),
      harmonized = c(
        "late_remodeling_state",
        "early_state_persistence",
        "ordinal_wound_state_position",
        "off_trajectory_ratio",
        "wound_activation",
        "z_Skin",
        "z_Wound30"
      ),
      stringsAsFactors = FALSE
    )

    x <- x[
      x$cell_type ==
        "Fibroblast" &
        x$metric %in%
        metric_map$original,
      ,
      drop = FALSE
    ]

    if (nrow(
      x
    )) {
      harmonized_metric <-
        metric_map$harmonized[
          match(
            x$metric,
            metric_map$original
          )
        ]

      cross_disease_rows[[
        "GSE163973"
      ]] <- data.frame(
        cohort =
          "GSE163973",
        disease =
          "keloid",
        analysis_set =
          "author_annotated_fibroblast",
        metric =
          harmonized_metric,
        mean_difference =
          x$mean_difference,
        hedges_g =
          x$hedges_g,
        cliffs_delta =
          x$cliffs_delta,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (file.exists(
  STAGE6B_RDS
)) {
  stage6b <- readRDS(
    STAGE6B_RDS
  )

  if ("comparisons" %in%
      names(
        stage6b
      )) {
    x <- stage6b$comparisons

    x <- x[
      x$analysis_set ==
        "high_specificity" &
        x$metric %in%
        primary_metrics,
      ,
      drop = FALSE
    ]

    if (nrow(
      x
    )) {
      cross_disease_rows[[
        "GSE181316"
      ]] <- data.frame(
        cohort =
          "GSE181316",
        disease =
          "keloid",
        analysis_set =
          "high_specificity",
        metric =
          x$metric,
        mean_difference =
          x$mean_difference,
        hedges_g =
          x$hedges_g,
        cliffs_delta =
          x$cliffs_delta,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (file.exists(
  STAGE7_RDS
)) {
  stage7 <- readRDS(
    STAGE7_RDS
  )

  if ("keloid_scar_results" %in%
      names(
        stage7
      )) {
    x <- stage7$keloid_scar_results

    x <- x[
      x$analysis_set ==
        "high_specificity" &
        x$metric %in%
        primary_metrics,
      ,
      drop = FALSE
    ]

    if (nrow(
      x
    )) {
      cross_disease_rows[[
        "GSE220300"
      ]] <- data.frame(
        cohort =
          "GSE220300",
        disease =
          "keloid",
        analysis_set =
          "high_specificity",
        metric =
          x$metric,
        mean_difference =
          x$mean_difference,
        hedges_g =
          NA_real_,
        cliffs_delta =
          x$cliffs_delta,
        stringsAsFactors = FALSE
      )
    }
  }
}

cross_disease_summary <- do.call(
  rbind,
  cross_disease_rows
)

rownames(
  cross_disease_summary
) <- NULL

safe_write_csv(
  cross_disease_summary,
  file.path(
    REPORT_DIR,
    "12_cross_disease_failure_mode_summary.csv"
  )
)

# ----------------------------- 结果判定 ----------------------------------------
get_primary_row <- function(metric) {
  x <- comparisons[
    comparisons$analysis_set ==
      "FB1_FB2_combined_primary" &
      comparisons$metric ==
      metric,
    ,
    drop = FALSE
  ]

  if (nrow(
    x
  ) !=
      1L) {
    stop(
      "无法唯一获得主要比较：",
      metric
    )
  }

  x
}

late_row <- get_primary_row(
  "late_remodeling_state"
)

early_row <- get_primary_row(
  "early_state_persistence"
)

balance_row <- get_primary_row(
  "late_vs_early_balance"
)

offtraj_row <- get_primary_row(
  "off_trajectory_ratio"
)

ordinal_row <- get_primary_row(
  "ordinal_wound_state_position"
)

early_or_offtrajectory_signal <- (
  (
    early_row$mean_difference >
      0 &&
      (
        early_row$cliffs_delta >=
          0.33 ||
          early_row$hedges_g >=
          0.5
      )
  ) ||
    (
      offtraj_row$mean_difference >
        0 &&
        (
          offtraj_row$cliffs_delta >=
            0.33 ||
            offtraj_row$hedges_g >=
            0.5
        )
    )
)

late_dominant_signal <- (
  late_row$mean_difference >
    0 &&
    late_row$cliffs_delta >=
    0.56 &&
    late_row$mean_difference >=
    early_row$mean_difference
)

distinct_failure_balance <- (
  balance_row$mean_difference <
    0
)

if (
  early_or_offtrajectory_signal &&
    distinct_failure_balance &&
    !late_dominant_signal
) {
  project_gate <-
    "PASS_DISTINCT_CHRONIC_WOUND_FAILURE_MODE"

  interpretation <- paste(
    "静脉溃疡成纤维细胞主要表现为早期状态持续和/或偏离正常轨迹，",
    "而不是瘢痕疙瘩中稳定的D30-like晚期重塑优势；",
    "支持不同慢性伤口具有不同修复终止失败模式。"
  )
} else if (
  late_dominant_signal &&
    !distinct_failure_balance
) {
  project_gate <-
    "SHARED_LATE_REMODELING_PATTERN"

  interpretation <- paste(
    "静脉溃疡同样表现为明显late-remodeling dominance；",
    "晚期重塑持续可能不是瘢痕疙瘩特异，而是部分慢性伤口共有现象。"
  )
} else if (
  early_or_offtrajectory_signal ||
    late_dominant_signal
) {
  project_gate <-
    "MIXED_CHRONIC_WOUND_PATTERN"

  interpretation <- paste(
    "静脉溃疡呈现混合修复失败状态；",
    "可作为疾病对照描述，但不能支持明确的瘢痕疙瘩特异性。"
  )
} else {
  project_gate <-
    "NO_CLEAR_VENOUS_ULCER_PATTERN"

  interpretation <- paste(
    "静脉溃疡队列未呈现稳定的固定伤口状态差异；",
    "本阶段不提供明确的疾病特异性证据。"
  )
}

decision_criteria <- data.frame(
  criterion = c(
    "early_or_offtrajectory_signal",
    "late_dominant_signal",
    "distinct_failure_balance"
  ),
  passed = c(
    early_or_offtrajectory_signal,
    late_dominant_signal,
    distinct_failure_balance
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  decision_criteria,
  file.path(
    REPORT_DIR,
    "13_reverse_control_criteria.csv"
  )
)

# ----------------------------- 图形 --------------------------------------------
# 1. 主要患者级指标
primary_scores <- sample_scores[
  sample_scores$analysis_set ==
    "FB1_FB2_combined_primary",
  ,
  drop = FALSE
]

pdf(
  file.path(
    FIG_DIR,
    "Figure_patient_level_reverse_control_metrics.pdf"
  ),
  width = 11,
  height = 8
)

op <- par(
  no.readonly = TRUE
)

par(
  mfrow = c(
    2,
    2
  ),
  mar = c(
    6,
    4,
    3,
    1
  )
)

plot_metric <- function(
  data,
  metric,
  ylab
) {
  x_position <- ifelse(
    data$group ==
      "normal_skin",
    1,
    2
  )

  y_range <- range(
    data[[metric]],
    finite = TRUE
  )

  y_padding <- max(
    0.1,
    diff(
      y_range
    ) *
      0.15
  )

  plot(
    jitter(
      x_position,
      amount = 0.05
    ),
    data[[metric]],
    pch = ifelse(
      data$group ==
        "normal_skin",
      1,
      16
    ),
    xlim = c(
      0.5,
      2.5
    ),
    ylim = c(
      y_range[1] -
        y_padding,
      y_range[2] +
        y_padding
    ),
    xaxt = "n",
    xlab = "",
    ylab =
      ylab
  )

  axis(
    1,
    at = c(
      1,
      2
    ),
    labels = c(
      "Normal skin",
      "Venous ulcer"
    ),
    las = 2
  )

  text(
    jitter(
      x_position,
      amount = 0.05
    ),
    data[[metric]],
    labels =
      data$sample_key,
    pos = 3,
    cex = 0.7
  )

  segments(
    0.85,
    stats::median(
      data[[metric]][
        data$group ==
          "normal_skin"
      ]
    ),
    1.15,
    stats::median(
      data[[metric]][
        data$group ==
          "normal_skin"
      ]
    ),
    lwd = 2
  )

  segments(
    1.85,
    stats::median(
      data[[metric]][
        data$group ==
          "venous_ulcer"
      ]
    ),
    2.15,
    stats::median(
      data[[metric]][
        data$group ==
          "venous_ulcer"
      ]
    ),
    lwd = 2
  )
}

plot_metric(
  primary_scores,
  "late_remodeling_state",
  "Late-remodeling state"
)

plot_metric(
  primary_scores,
  "early_state_persistence",
  "Early-state persistence"
)

plot_metric(
  primary_scores,
  "late_vs_early_balance",
  "Late minus early balance"
)

plot_metric(
  primary_scores,
  "off_trajectory_ratio",
  "Off-trajectory ratio"
)

par(
  op
)

dev.off()

# 2. 早期-晚期失败模式平面
pdf(
  file.path(
    FIG_DIR,
    "Figure_early_vs_late_failure_mode_plane.pdf"
  ),
  width = 8,
  height = 7
)

plot(
  primary_scores$early_state_persistence,
  primary_scores$late_remodeling_state,
  pch = ifelse(
    primary_scores$group ==
      "normal_skin",
    1,
    16
  ),
  xlab =
    "Early-state persistence",
  ylab =
    "Late-remodeling state",
  main =
    "Venous ulcer fibroblast failure-mode plane"
)

text(
  primary_scores$early_state_persistence,
  primary_scores$late_remodeling_state,
  labels =
    primary_scores$sample_key,
  pos = 3,
  cex = 0.75
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

dev.off()

# 3. 四阶段固定评分热图
heat_order <- c(
  sample_design$sample_id[
    sample_design$group ==
      "normal_skin"
  ],
  sample_design$sample_id[
    sample_design$group ==
      "venous_ulcer"
  ]
)

heat_data <- primary_scores[
  match(
    heat_order,
    primary_scores$sample_key
  ),
  ,
  drop = FALSE
]

heat_matrix <- as.matrix(
  heat_data[
    ,
    c(
      "z_Skin",
      "z_Wound1",
      "z_Wound7",
      "z_Wound30"
    )
  ]
)

rownames(
  heat_matrix
) <- paste(
  heat_data$sample_key,
  heat_data$group,
  sep = " | "
)

colnames(
  heat_matrix
) <- c(
  "Skin",
  "D1",
  "D7",
  "D30"
)

pdf(
  file.path(
    FIG_DIR,
    "Figure_GSE265972_fixed_stage_score_heatmap.pdf"
  ),
  width = 7,
  height = 8
)

par(
  mar = c(
    5,
    10,
    3,
    2
  )
)

image(
  x = seq_len(
    ncol(
      heat_matrix
    )
  ),
  y = seq_len(
    nrow(
      heat_matrix
    )
  ),
  z = t(
    heat_matrix
  ),
  axes = FALSE,
  xlab = "",
  ylab = "",
  main =
    "GSE265972 fixed wound-state scores"
)

axis(
  1,
  at = seq_len(
    ncol(
      heat_matrix
    )
  ),
  labels =
    colnames(
      heat_matrix
    )
)

axis(
  2,
  at = seq_len(
    nrow(
      heat_matrix
    )
  ),
  labels =
    rownames(
      heat_matrix
    ),
  las = 2,
  cex.axis = 0.75
)

box()

dev.off()

# 4. 跨疾病效应方向图
late_cross <- cross_disease_summary[
  cross_disease_summary$metric ==
    "late_remodeling_state",
  ,
  drop = FALSE
]

early_cross <- cross_disease_summary[
  cross_disease_summary$metric ==
    "early_state_persistence",
  ,
  drop = FALSE
]

if (
  nrow(
    late_cross
  ) >
    0L &&
    nrow(
      early_cross
    ) >
    0L
) {
  plot_table <- merge(
    late_cross[
      ,
      c(
        "cohort",
        "disease",
        "mean_difference"
      )
    ],
    early_cross[
      ,
      c(
        "cohort",
        "mean_difference"
      )
    ],
    by = "cohort",
    suffixes = c(
      "_late",
      "_early"
    ),
    all = TRUE
  )

  pdf(
    file.path(
      FIG_DIR,
      "Figure_cross_disease_failure_mode_effects.pdf"
    ),
    width = 10,
    height = 6
  )

  effect_matrix <- rbind(
    Late_remodeling =
      plot_table$mean_difference_late,
    Early_persistence =
      plot_table$mean_difference_early
  )

  barplot(
    effect_matrix,
    beside = TRUE,
    names.arg =
      plot_table$cohort,
    ylab =
      "Disease-control mean difference",
    main =
      "Keloid and venous-ulcer wound-state effects",
    las = 2
  )

  abline(
    h = 0,
    lty = 2
  )

  legend(
    "topright",
    legend = rownames(
      effect_matrix
    ),
    pch = 15,
    bty = "n"
  )

  dev.off()
}

# ----------------------------- 保存对象 ----------------------------------------
saveRDS(
  list(
    project_gate =
      project_gate,
    interpretation =
      interpretation,
    sample_design =
      sample_design,
    sample_qc_summary =
      sample_qc_summary,
    barcode_summary =
      barcode_summary,
    fibro_summary =
      fibro_summary,
    combined_primary =
      combined_result,
    FB1_sensitivity =
      fb1_result,
    FB2_sensitivity =
      fb2_result,
    sample_scores =
      sample_scores,
    comparisons =
      comparisons,
    cross_disease_summary =
      cross_disease_summary,
    decision_criteria =
      decision_criteria,
    fixed_signatures =
      signatures
  ),
  file = file.path(
    OBJECT_DIR,
    "GSE265972_venous_ulcer_reverse_control.rds"
  ),
  compress = "gzip"
)

# ----------------------------- 最终结论 ----------------------------------------
decision_lines <- c(
  "第九阶段：GSE265972静脉溃疡疾病对照结论",
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
  "结果解释：",
  interpretation,
  "",
  "主要FB1+FB2结果：",
  paste0(
    "- late-remodeling mean difference = ",
    round(
      late_row$mean_difference,
      4
    )
  ),
  paste0(
    "- late-remodeling Hedges g = ",
    round(
      late_row$hedges_g,
      3
    )
  ),
  paste0(
    "- late-remodeling Cliff delta = ",
    round(
      late_row$cliffs_delta,
      3
    )
  ),
  paste0(
    "- early persistence mean difference = ",
    round(
      early_row$mean_difference,
      4
    )
  ),
  paste0(
    "- early persistence Hedges g = ",
    round(
      early_row$hedges_g,
      3
    )
  ),
  paste0(
    "- early persistence Cliff delta = ",
    round(
      early_row$cliffs_delta,
      3
    )
  ),
  paste0(
    "- late-vs-early balance mean difference = ",
    round(
      balance_row$mean_difference,
      4
    )
  ),
  paste0(
    "- off-trajectory mean difference = ",
    round(
      offtraj_row$mean_difference,
      4
    )
  ),
  paste0(
    "- ordinal position mean difference = ",
    round(
      ordinal_row$mean_difference,
      4
    )
  ),
  "",
  "解释边界：",
  "- 正常对照为正常皮肤，不是成熟瘢痕。",
  "- GSE265972是疾病对照，不参与瘢痕疙瘩主要效应合并。",
  "- 主要集合为作者注释FB1+FB2；FB1和FB2分别为敏感性分析。",
  "- metadata已经过作者QC和双细胞处理，本阶段不额外删除作者保留细胞。",
  "- 跨疾病只比较效应方向和失败模式，不比较绝对数值大小。",
  "- 若静脉溃疡也表现为late-remodeling dominance，必须降低瘢痕疙瘩特异性表述。",
  "- ordinal position不是精确伤口日龄。"
)

safe_write_lines(
  decision_lines,
  file.path(
    REPORT_DIR,
    "14_STAGE9_DECISION.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    REPORT_DIR,
    "15_SESSION_INFO.txt"
  )
)

safe_write_lines(
  c(
    paste0(
      "STAGE9_COMPLETED=",
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
      "EARLY_OR_OFFTRAJECTORY_SIGNAL=",
      early_or_offtrajectory_signal
    ),
    paste0(
      "LATE_DOMINANT_SIGNAL=",
      late_dominant_signal
    ),
    paste0(
      "DISTINCT_FAILURE_BALANCE=",
      distinct_failure_balance
    )
  ),
  file.path(
    REPORT_DIR,
    "STAGE9_COMPLETED.txt"
  )
)

if (dir.exists(
  WORK_DIR
)) {
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

if (dir.exists(
  package_dir
)) {
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
  if (file.exists(
    zip_path
  )) {
    file.remove(
      zip_path
    )
  }

  escape_path <- function(x) {
    gsub(
      "'",
      "''",
      normalize_slash(
        x
      ),
      fixed = TRUE
    )
  }

  command <- paste0(
    "$ErrorActionPreference='Stop'; ",
    "$items=Get-ChildItem -LiteralPath '",
    escape_path(
      source_dir
    ),
    "'; ",
    "Compress-Archive -Path $items.FullName ",
    "-DestinationPath '",
    escape_path(
      zip_path
    ),
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

  file.exists(
    zip_path
  ) &&
    file.info(
      zip_path
    )$size >
    0
}

zip_ok <- create_zip_windows(
  package_dir,
  PACKAGE_ZIP
)

if (!zip_ok) {
  if (file.exists(
    PACKAGE_TARGZ
  )) {
    file.remove(
      PACKAGE_TARGZ
    )
  }

  old_wd <- getwd()

  setwd(
    STAGE_DIR
  )

  try(
    utils::tar(
      tarfile =
        PACKAGE_TARGZ,
      files =
        basename(
          package_dir
        ),
      compression =
        "gzip",
      tar =
        "internal"
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
  "第九阶段完成"
)

cat(
  "\n============================================================\n"
)

cat(
  "第九阶段运行完成。\n"
)

cat(
  "本地分析对象保存在：",
  file.path(
    OBJECT_DIR,
    "GSE265972_venous_ulcer_reverse_control.rds"
  ),
  "\n",
  sep = ""
)

if (file.exists(
  PACKAGE_ZIP
)) {
  cat(
    "请上传：",
    PACKAGE_ZIP,
    "\n",
    sep = ""
  )
} else if (file.exists(
  PACKAGE_TARGZ
)) {
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

cat(
  "============================================================\n"
)
