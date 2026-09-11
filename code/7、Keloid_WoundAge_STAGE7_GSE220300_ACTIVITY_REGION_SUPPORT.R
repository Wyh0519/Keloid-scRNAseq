# =============================================================================
# 项目：瘢痕疙瘩“晚期重塑状态持续与伤口状态终止失败”
# 第七阶段：GSE220300活动度、中心/周边与成熟瘢痕患者内支持分析
#
# 根目录：
#   C:/Users/33652/Desktop/LZZ/1、yssj
#
# 数据结构：
#   Active keloid:
#     A1 = AC01 + AP01
#     A2 = AC02 + AP02
#
#   Inactive keloid:
#     I1 = IC01 + IP01
#     I2 = IC02 + IP02
#
#   Mature scar:
#     MS01来自I1患者
#     MS02来自I2患者
#
#   Normal scar:
#     NS01，独立瘢痕参照
#
# 分析定位：
#   - 本阶段是活动度、区域和患者内成熟瘢痕支持分析；
#   - 不是第三个正式独立验证队列；
#   - 主要关注固定的late-remodeling state；
#   - active vs inactive及center vs periphery仅作探索性描述；
#   - I1/I2病灶与同患者成熟瘢痕为关键患者内支持；
#   - 所有成纤维细胞标志物和伤口阶段签名均来自前面已经冻结的对象。
#
# 冻结终点：
#   Primary supportive:
#     late_remodeling_state
#     = D30 score - mean(D1 score, D7 score)
#
#   Key secondary:
#     ordinal_wound_state_position
#
#   Supporting:
#     D30-like score
#     wound activation
#     Skin-like score / skin-return state
#
#   Exploratory:
#     active vs inactive
#     center vs periphery
#     off-trajectory ratio
#     early-state persistence
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANG = "en")

ROOT_DIR <- "C:/Users/33652/Desktop/LZZ/1、yssj"

RAW_TAR <- file.path(
  ROOT_DIR,
  "GSE220300_RAW.tar"
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

STAGE6A_RDS <- file.path(
  ROOT_DIR,
  "06A_STAGE6A_GSE181316_CELL_CALLING",
  "objects",
  "GSE181316_fibroblast_like_pseudobulk_list.rds"
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

STAGE_DIR <- file.path(
  ROOT_DIR,
  "07_STAGE7_GSE220300_ACTIVITY_REGION"
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
  "第七阶段_GSE220300活动度区域与成熟瘢痕支持分析检查包.zip"
)

PACKAGE_TARGZ <- file.path(
  ROOT_DIR,
  "第七阶段_GSE220300活动度区域与成熟瘢痕支持分析检查包.tar.gz"
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
  "00_stage7_log.txt"
)

WARNING_FILE <- file.path(
  REPORT_DIR,
  "00_stage7_warnings.txt"
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
    "AC01",
    "AP01",
    "AC02",
    "AP02",
    "IC01",
    "IP01",
    "MS01",
    "IC02",
    "IP02",
    "MS02",
    "NS01"
  ),
  patient_id = c(
    "A1",
    "A1",
    "A2",
    "A2",
    "I1",
    "I1",
    "I1",
    "I2",
    "I2",
    "I2",
    "NS1"
  ),
  tissue_state = c(
    "active_keloid",
    "active_keloid",
    "active_keloid",
    "active_keloid",
    "inactive_keloid",
    "inactive_keloid",
    "mature_scar",
    "inactive_keloid",
    "inactive_keloid",
    "mature_scar",
    "normal_scar"
  ),
  activity = c(
    "active",
    "active",
    "active",
    "active",
    "inactive",
    "inactive",
    "resolved_scar",
    "inactive",
    "inactive",
    "resolved_scar",
    "resolved_scar"
  ),
  region = c(
    "center",
    "periphery",
    "center",
    "periphery",
    "center",
    "periphery",
    "mature_scar",
    "center",
    "periphery",
    "mature_scar",
    "normal_scar"
  ),
  lesion_patient = c(
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    TRUE,
    TRUE,
    FALSE,
    FALSE
  ),
  paired_mature_scar_patient = c(
    NA,
    NA,
    NA,
    NA,
    "I1",
    "I1",
    "I1",
    "I2",
    "I2",
    "I2",
    NA
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

extract_sample_id <- function(
  matrix_path
) {
  x <- basename(
    matrix_path
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

  toupper(
    x
  )
}

locate_triplet <- function(
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

  matrix_path <- files[
    grepl(
      paste0(
        "_",
        sample_id,
        "_matrix\\.mtx\\.gz$"
      ),
      basename(
        files
      ),
      ignore.case = TRUE
    )
  ]

  barcode_path <- files[
    grepl(
      paste0(
        "_",
        sample_id,
        "_barcodes\\.tsv\\.gz$"
      ),
      basename(
        files
      ),
      ignore.case = TRUE
    )
  ]

  feature_path <- files[
    grepl(
      paste0(
        "_",
        sample_id,
        "_features\\.tsv\\.gz$"
      ),
      basename(
        files
      ),
      ignore.case = TRUE
    )
  ]

  if (
    length(
      matrix_path
    ) != 1L ||
      length(
        barcode_path
      ) != 1L ||
      length(
        feature_path
      ) != 1L
  ) {
    stop(
      "无法唯一定位",
      sample_id,
      "的matrix/barcodes/features。"
    )
  }

  list(
    matrix = matrix_path,
    barcodes = barcode_path,
    features = feature_path
  )
}

canonicalize_named_vector <- function(x) {
  if (is.null(
    names(
      x
    )
  )) {
    stop(
      "pseudobulk向量缺少基因名。"
    )
  }

  gene <- toupper(
    trimws(
      as.character(
        names(
          x
        )
      )
    )
  )

  value <- as.numeric(
    x
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

  out <- as.numeric(
    summed[, 1]
  )

  names(
    out
  ) <- rownames(
    summed
  )

  out
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

  out <- as.numeric(
    summed[, 1]
  )

  names(
    out
  ) <- rownames(
    summed
  )

  out
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
    nrow = length(
      markers
    ),
    ncol = ncol(
      counts
    ),
    dimnames = list(
      markers,
      colnames(
        counts
      )
    )
  )

  for (
    i in seq_along(
      markers
    )
  ) {
    idx <- which(
      gene_upper ==
        markers[i]
    )

    if (length(
      idx
    )) {
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

merge_named_vector_list <- function(
  vector_list,
  required_names
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

  vector_list <- lapply(
    vector_list,
    canonicalize_named_vector
  )

  common_genes <- Reduce(
    intersect,
    lapply(
      vector_list,
      names
    )
  )

  if (length(
    common_genes
  ) < 10000L) {
    warn_msg(
      "GSE220300样本间共同基因少于10,000：",
      length(
        common_genes
      )
    )
  }

  matrix_result <- do.call(
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
    matrix_result
  ) <- common_genes

  colnames(
    matrix_result
  ) <- names(
    vector_list
  )

  storage.mode(
    matrix_result
  ) <- "double"

  matrix_result
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

  out <- sample_info[
    idx,
    ,
    drop = FALSE
  ]

  if (!identical(
    colnames(
      counts
    ),
    out$key
  )) {
    stop(
      "正常伤口pseudobulk顺序匹配失败。"
    )
  }

  out
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
    n_genes < 20L
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
    library_size <= 0
  )) {
    stop(
      "存在文库总计数≤0。"
    )
  }

  log2(
    sweep(
      counts + 0.5,
      2,
      library_size + 1,
      "/"
    ) * 1e6
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
      ) < 15L) {
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
    new_z = new_z,
    center = center,
    scale = scale_value
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
  distance <- numeric(
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

    distance[i] <- sqrt(
      sum(
        (
          reference_z[i, ] -
            centroid
        ) ^ 2
      )
    )
  }

  distance
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
            ) ^ 2
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
        distance ^ 2
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
        weights["Skin"],
      weight_Wound1 =
        weights["Wound1"],
      weight_Wound7 =
        weights["Wound7"],
      weight_Wound30 =
        weights["Wound30"],
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

aggregate_samples <- function(
  sample_counts,
  aggregation_design
) {
  result <- do.call(
    cbind,
    lapply(
      aggregation_design$aggregate_id,
      function(aggregate_id) {
        sample_ids <- strsplit(
          aggregation_design$component_samples[
            aggregation_design$aggregate_id ==
              aggregate_id
          ][1],
          "\\+"
        )[[1]]

        Matrix::rowSums(
          sample_counts[
            ,
            sample_ids,
            drop = FALSE
          ]
        )
      }
    )
  )

  colnames(
    result
  ) <- aggregation_design$aggregate_id

  result
}

descriptive_two_group <- function(
  value,
  group,
  group_a,
  group_b
) {
  x <- value[
    group ==
      group_a &
      is.finite(
        value
      )
  ]

  y <- value[
    group ==
      group_b &
      is.finite(
        value
      )
  ]

  if (
    length(
      x
    ) < 1L ||
      length(
        y
      ) < 1L
  ) {
    return(
      data.frame(
        n_group_a = length(
          x
        ),
        n_group_b = length(
          y
        ),
        mean_group_a = NA_real_,
        mean_group_b = NA_real_,
        mean_difference = NA_real_,
        median_difference = NA_real_,
        cliffs_delta = NA_real_,
        stringsAsFactors = FALSE
      )
    )
  }

  pair_difference <- outer(
    x,
    y,
    "-"
  )

  data.frame(
    n_group_a = length(
      x
    ),
    n_group_b = length(
      y
    ),
    mean_group_a = mean(
      x
    ),
    mean_group_b = mean(
      y
    ),
    mean_difference = mean(
      x
    ) -
      mean(
        y
      ),
    median_difference = stats::median(
      x
    ) -
      stats::median(
        y
      ),
    cliffs_delta = mean(
      pair_difference >
        0
    ) -
      mean(
        pair_difference <
          0
      ),
    stringsAsFactors = FALSE
  )
}

paired_difference_summary <- function(
  data,
  id_column,
  condition_column,
  value_column,
  condition_high,
  condition_low
) {
  ids <- unique(
    data[[id_column]]
  )

  rows <- list()

  for (
    id in ids
  ) {
    x <- data[
      data[[id_column]] ==
        id,
      ,
      drop = FALSE
    ]

    high_value <- x[[value_column]][
      x[[condition_column]] ==
        condition_high
    ]

    low_value <- x[[value_column]][
      x[[condition_column]] ==
        condition_low
    ]

    if (
      length(
        high_value
      ) ==
      1L &&
        length(
          low_value
        ) ==
        1L
    ) {
      rows[[length(
        rows
      ) + 1L]] <- data.frame(
        pair_id = id,
        condition_high =
          condition_high,
        condition_low =
          condition_low,
        high_value =
          high_value,
        low_value =
          low_value,
        paired_difference =
          high_value -
          low_value,
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(
    rows
  )) {
    return(
      data.frame()
    )
  }

  do.call(
    rbind,
    rows
  )
}

analyse_set <- function(
  set_name,
  vector_list,
  sample_design,
  reference_counts,
  reference_design,
  signatures
) {
  log_msg(
    "固定伤口状态映射：",
    set_name
  )

  sample_counts <- merge_named_vector_list(
    vector_list,
    required_names =
      sample_design$sample_id
  )

  common_genes <- intersect(
    rownames(
      reference_counts
    ),
    rownames(
      sample_counts
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

  # 患者级病灶聚合设计。
  patient_lesion_design <- data.frame(
    aggregate_id = c(
      "A1_lesion",
      "A2_lesion",
      "I1_lesion",
      "I2_lesion",
      "I1_mature_scar",
      "I2_mature_scar",
      "NS1_normal_scar"
    ),
    component_samples = c(
      "AC01+AP01",
      "AC02+AP02",
      "IC01+IP01",
      "IC02+IP02",
      "MS01",
      "MS02",
      "NS01"
    ),
    patient_id = c(
      "A1",
      "A2",
      "I1",
      "I2",
      "I1",
      "I2",
      "NS1"
    ),
    aggregate_state = c(
      "active_lesion",
      "active_lesion",
      "inactive_lesion",
      "inactive_lesion",
      "mature_scar",
      "mature_scar",
      "normal_scar"
    ),
    stringsAsFactors = FALSE
  )

  aggregate_counts <- aggregate_samples(
    sample_counts,
    patient_lesion_design
  )

  reference_expression <- log_cpm(
    reference_counts[
      common_genes,
      ,
      drop = FALSE
    ]
  )

  sample_expression <- log_cpm(
    sample_counts[
      common_genes,
      ,
      drop = FALSE
    ]
  )

  aggregate_expression <- log_cpm(
    aggregate_counts[
      common_genes,
      ,
      drop = FALSE
    ]
  )

  reference_scores <- score_samples_by_rank(
    reference_expression,
    signatures
  )

  sample_scores <- score_samples_by_rank(
    sample_expression,
    signatures
  )

  aggregate_scores <- score_samples_by_rank(
    aggregate_expression,
    signatures
  )

  standardized_sample <-
    standardize_by_reference(
      reference_scores,
      sample_scores
    )

  standardized_aggregate <-
    standardize_by_reference(
      reference_scores,
      aggregate_scores
    )

  reference_z <-
    standardized_sample$reference_z

  sample_z <-
    standardized_sample$new_z

  aggregate_z <-
    standardized_aggregate$new_z

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

  sample_mapping <- map_to_reference(
    sample_z,
    centroids,
    reference_distance_limit
  )

  sample_mapping <- add_fixed_contrasts(
    sample_mapping,
    sample_z
  )

  sample_mapping$analysis_set <-
    set_name

  sample_mapping <- merge(
    sample_mapping,
    sample_design,
    by.x = "sample_key",
    by.y = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )

  sample_mapping <- sample_mapping[
    match(
      rownames(
        sample_z
      ),
      sample_mapping$sample_key
    ),
    ,
    drop = FALSE
  ]

  aggregate_mapping <- map_to_reference(
    aggregate_z,
    centroids,
    reference_distance_limit
  )

  aggregate_mapping <- add_fixed_contrasts(
    aggregate_mapping,
    aggregate_z
  )

  aggregate_mapping$analysis_set <-
    set_name

  aggregate_mapping <- merge(
    aggregate_mapping,
    patient_lesion_design,
    by.x = "sample_key",
    by.y = "aggregate_id",
    all.x = TRUE,
    sort = FALSE
  )

  aggregate_mapping <- aggregate_mapping[
    match(
      rownames(
        aggregate_z
      ),
      aggregate_mapping$sample_key
    ),
    ,
    drop = FALSE
  ]

  coverage <- data.frame(
    analysis_set = set_name,
    stage = names(
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
    sample_counts = sample_counts,
    aggregate_counts =
      aggregate_counts,
    sample_mapping =
      sample_mapping,
    aggregate_mapping =
      aggregate_mapping,
    coverage = coverage,
    calibration = list(
      common_genes =
        common_genes,
      signatures =
        signatures,
      reference_center =
        standardized_sample$center,
      reference_scale =
        standardized_sample$scale,
      stage_centroids =
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
  "第七阶段开始"
)

required_files <- c(
  RAW_TAR,
  REFERENCE_PB_RDS,
  REFERENCE_MODEL_RDS,
  STAGE6A_RDS
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

stage6a_object <- readRDS(
  STAGE6A_RDS
)

if (!"marker_definition" %in%
    names(
      stage6a_object
    )) {
  stop(
    "第六阶段A对象缺少marker_definition。"
  )
}

marker_definition <-
  stage6a_object$marker_definition

required_marker_sets <- c(
  "fibro_core",
  "fibro_support",
  "immune",
  "epithelial",
  "endothelial",
  "neural"
)

marker_sets <- setNames(
  lapply(
    required_marker_sets,
    function(marker_set) {
      unique(
        toupper(
          trimws(
            as.character(
              marker_definition$gene[
                marker_definition$marker_set ==
                  marker_set
              ]
            )
          )
        )
      )
    }
  ),
  required_marker_sets
)

if (any(
  vapply(
    marker_sets,
    length,
    integer(1)
  ) ==
    0L
)) {
  stop(
    "至少一个冻结标志物集合为空。"
  )
}

safe_write_csv(
  marker_definition,
  file.path(
    REPORT_DIR,
    "02_reused_frozen_marker_definition.csv"
  )
)

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
  ) != 12L
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

signatures <-
  get_stable_fibroblast_signatures(
    reference_model
  )

signature_summary <- data.frame(
  stage = names(
    signatures
  ),
  n_stable_genes = vapply(
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
    "03_fixed_signature_summary.csv"
  )
)

# ----------------------------- 解包 --------------------------------------------
outer_dir <- file.path(
  WORK_DIR,
  "GSE220300_RAW"
)

dir.create(
  outer_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

log_msg(
  "解包GSE220300_RAW.tar"
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

detected_ids <- vapply(
  matrix_files,
  extract_sample_id,
  character(1)
)

if (
  length(
    matrix_files
  ) != 11L ||
    !setequal(
      detected_ids,
      sample_design$sample_id
    )
) {
  stop(
    "GSE220300样本结构与预期不一致。实际：",
    paste(
      detected_ids,
      collapse = ", "
    )
  )
}

# ----------------------------- 逐样本QC与成纤维细胞提取 ------------------------
qc_rows <- list()
marker_rows <- list()
feature_rows <- list()
fibro_high_list <- list()
fibro_broad_list <- list()

gene_reference_id <- NULL
gene_reference_symbol <- NULL

pdf(
  file.path(
    FIG_DIR,
    "Figure_GSE220300_QC_and_marker_audit.pdf"
  ),
  width = 10,
  height = 8,
  onefile = TRUE
)

for (
  sample_id in sample_design$sample_id
) {
  log_msg(
    "处理样本：",
    sample_id
  )

  triplet <- locate_triplet(
    outer_dir,
    sample_id
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

  raw_barcodes <- as.character(
    barcodes[[1]]
  )

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
        raw_barcodes
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
    gene_reference_id <- gene_id
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
  ) <- raw_barcodes

  feature_rows[[length(
    feature_rows
  ) + 1L]] <- data.frame(
    sample_id = sample_id,
    n_features = length(
      gene_id
    ),
    identical_gene_ids_to_first =
      identical_id,
    identical_gene_symbols_to_first =
      identical_symbol,
    stringsAsFactors = FALSE
  )

  n_count <- Matrix::colSums(
    counts
  )

  n_feature <- Matrix::colSums(
    counts > 0
  )

  gene_upper <- toupper(
    rownames(
      counts
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
      counts[
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
        counts
      )
    )
  }

  # 采用该队列原研究公开的QC范围：
  # UMI >500，基因>200且<6500，线粒体比例<15%。
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
      n_count > 500 &
      n_feature > 200 &
      n_feature < 6500 &
      pct_mt < 15
  )

  qc_keep[
    is.na(
      qc_keep
    )
  ] <- FALSE

  if (
    sum(
      qc_keep
    ) <
      100L
  ) {
    stop(
      sample_id,
      "QC后细胞少于100：",
      sum(
        qc_keep
      )
    )
  }

  counts_qc <- counts[
    ,
    qc_keep,
    drop = FALSE
  ]

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
      qc_keep
    ]
  )

  marker_results <- list()

  for (
    marker_name in names(
      marker_sets
    )
  ) {
    marker_counts <- marker_count_matrix(
      counts_qc,
      gene_reference_symbol,
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
    marker_results$fibro_core$detected >=
      2 &
      marker_results$fibro_core$score >
      max_exclusion_score
  )

  high_specificity_fibroblast <- (
    marker_results$fibro_core$detected >=
      3 &
      (
        marker_results$fibro_support$detected >=
          1 |
          marker_results$fibro_core$detected >=
          4
      ) &
      marker_results$fibro_core$score >
      max_exclusion_score +
      0.15 &
      marker_results$immune$detected <=
      1 &
      marker_results$epithelial$detected <=
      1 &
      marker_results$endothelial$detected <=
      1 &
      marker_results$neural$detected <=
      1
  )

  broad_fibroblast[
    is.na(
      broad_fibroblast
    )
  ] <- FALSE

  high_specificity_fibroblast[
    is.na(
      high_specificity_fibroblast
    )
  ] <- FALSE

  high_specificity_fibroblast <-
    high_specificity_fibroblast &
    broad_fibroblast

  n_high <- sum(
    high_specificity_fibroblast
  )

  n_broad <- sum(
    broad_fibroblast
  )

  if (
    n_high <
      50L
  ) {
    warn_msg(
      sample_id,
      "高特异性fibroblast-like少于50：",
      n_high
    )
  }

  if (
    n_broad <
      100L
  ) {
    warn_msg(
      sample_id,
      "宽松fibroblast-like少于100：",
      n_broad
    )
  }

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

  fibro_high_list[[sample_id]] <-
    collapse_vector_by_symbol(
      pb_high_raw,
      gene_reference_symbol
    )

  fibro_broad_list[[sample_id]] <-
    collapse_vector_by_symbol(
      pb_broad_raw,
      gene_reference_symbol
    )

  sample_row <- sample_design[
    sample_design$sample_id ==
      sample_id,
    ,
    drop = FALSE
  ]

  qc_rows[[length(
    qc_rows
  ) + 1L]] <- data.frame(
    sample_id = sample_id,
    patient_id =
      sample_row$patient_id,
    tissue_state =
      sample_row$tissue_state,
    region =
      sample_row$region,
    n_input_barcodes = ncol(
      counts
    ),
    n_QC_keep = ncol(
      counts_qc
    ),
    QC_retention_rate = ncol(
      counts_qc
    ) /
      ncol(
        counts
      ),
    n_fibroblast_high =
      n_high,
    n_fibroblast_broad =
      n_broad,
    fibroblast_high_fraction =
      n_high /
      ncol(
        counts_qc
      ),
    fibroblast_broad_fraction =
      n_broad /
      ncol(
        counts_qc
      ),
    median_nCount =
      stats::median(
        n_count_qc
      ),
    median_nFeature =
      stats::median(
        n_feature_qc
      ),
    median_percent_mt =
      stats::median(
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
    subset_index in seq_along(
      subset_names
    )
  ) {
    selected <- subset_indices[[
      subset_index
    ]]

    if (!any(
      selected
    )) {
      next
    }

    marker_rows[[length(
      marker_rows
    ) + 1L]] <- data.frame(
      sample_id =
        sample_id,
      tissue_state =
        sample_row$tissue_state,
      subset =
        subset_names[
          subset_index
        ],
      n_cells =
        sum(
          selected
        ),
      median_fibro_core_detected =
        stats::median(
          marker_results$fibro_core$detected[
            selected
          ]
        ),
      median_fibro_support_detected =
        stats::median(
          marker_results$fibro_support$detected[
            selected
          ]
        ),
      median_fibro_core_score =
        stats::median(
          marker_results$fibro_core$score[
            selected
          ]
        ),
      median_max_exclusion_score =
        stats::median(
          max_exclusion_score[
            selected
          ]
        ),
      stringsAsFactors = FALSE
    )
  }

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
      n_count + 1
    ),
    breaks = 50,
    main = paste0(
      sample_id,
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
    n_feature,
    breaks = 50,
    main = paste0(
      sample_id,
      " nFeature"
    ),
    xlab = "Detected genes"
  )

  abline(
    v = c(
      200,
      6500
    ),
    lty = 2
  )

  hist(
    pct_mt,
    breaks = 50,
    main = paste0(
      sample_id,
      " mitochondrial %"
    ),
    xlab = "percent.mt"
  )

  abline(
    v = 15,
    lty = 2
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
    xlab =
      "Fibroblast core score",
    ylab =
      "Maximum exclusion score",
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

  par(
    op
  )

  rm(
    counts,
    counts_qc,
    marker_results
  )

  gc(
    verbose = FALSE
  )
}

dev.off()

qc_summary <- do.call(
  rbind,
  qc_rows
)

marker_summary <- do.call(
  rbind,
  marker_rows
)

feature_summary <- do.call(
  rbind,
  feature_rows
)

rownames(
  qc_summary
) <- NULL

rownames(
  marker_summary
) <- NULL

rownames(
  feature_summary
) <- NULL

safe_write_csv(
  qc_summary,
  file.path(
    REPORT_DIR,
    "04_sample_QC_fibroblast_summary.csv"
  )
)

safe_write_csv(
  marker_summary,
  file.path(
    REPORT_DIR,
    "05_marker_identity_summary.csv"
  )
)

safe_write_csv(
  feature_summary,
  file.path(
    REPORT_DIR,
    "06_feature_consistency.csv"
  )
)

# ----------------------------- 固定映射 ----------------------------------------
high_result <- analyse_set(
  set_name =
    "high_specificity",
  vector_list =
    fibro_high_list,
  sample_design =
    sample_design,
  reference_counts =
    reference_counts,
  reference_design =
    reference_design,
  signatures =
    signatures
)

broad_result <- analyse_set(
  set_name =
    "broad_sensitivity",
  vector_list =
    fibro_broad_list,
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
  high_result$sample_mapping,
  broad_result$sample_mapping
)

aggregate_scores <- rbind(
  high_result$aggregate_mapping,
  broad_result$aggregate_mapping
)

coverage <- rbind(
  high_result$coverage,
  broad_result$coverage
)

rownames(
  sample_scores
) <- NULL

rownames(
  aggregate_scores
) <- NULL

rownames(
  coverage
) <- NULL

safe_write_csv(
  coverage,
  file.path(
    REPORT_DIR,
    "07_signature_gene_coverage.csv"
  )
)

safe_write_csv(
  sample_scores,
  file.path(
    REPORT_DIR,
    "08_sample_level_wound_state_scores.csv"
  )
)

safe_write_csv(
  aggregate_scores,
  file.path(
    REPORT_DIR,
    "09_patient_lesion_and_scar_scores.csv"
  )
)

# ----------------------------- 患者内与活动度分析 ------------------------------
metrics <- c(
  "late_remodeling_state",
  "ordinal_wound_state_position",
  "z_Wound30",
  "wound_activation",
  "z_Skin",
  "skin_return_state",
  "off_trajectory_ratio",
  "early_state_persistence"
)

region_pair_rows <- list()
inactive_scar_pair_rows <- list()
activity_rows <- list()
keloid_scar_rows <- list()

for (
  analysis_set in c(
    "high_specificity",
    "broad_sensitivity"
  )
) {
  sample_subset <- sample_scores[
    sample_scores$analysis_set ==
      analysis_set,
    ,
    drop = FALSE
  ]

  aggregate_subset <- aggregate_scores[
    aggregate_scores$analysis_set ==
      analysis_set,
    ,
    drop = FALSE
  ]

  for (
    metric in metrics
  ) {
    active_region <- paired_difference_summary(
      data = sample_subset[
        sample_subset$tissue_state ==
          "active_keloid",
        ,
        drop = FALSE
      ],
      id_column =
        "patient_id",
      condition_column =
        "region",
      value_column =
        metric,
      condition_high =
        "periphery",
      condition_low =
        "center"
    )

    if (nrow(
      active_region
    )) {
      active_region$analysis_set <-
        analysis_set
      active_region$activity_group <-
        "active"
      active_region$metric <-
        metric

      region_pair_rows[[length(
        region_pair_rows
      ) + 1L]] <-
        active_region
    }

    inactive_region <- paired_difference_summary(
      data = sample_subset[
        sample_subset$tissue_state ==
          "inactive_keloid",
        ,
        drop = FALSE
      ],
      id_column =
        "patient_id",
      condition_column =
        "region",
      value_column =
        metric,
      condition_high =
        "periphery",
      condition_low =
        "center"
    )

    if (nrow(
      inactive_region
    )) {
      inactive_region$analysis_set <-
        analysis_set
      inactive_region$activity_group <-
        "inactive"
      inactive_region$metric <-
        metric

      region_pair_rows[[length(
        region_pair_rows
      ) + 1L]] <-
        inactive_region
    }

    inactive_pair <- paired_difference_summary(
      data = aggregate_subset[
        aggregate_subset$aggregate_state %in%
          c(
            "inactive_lesion",
            "mature_scar"
          ),
        ,
        drop = FALSE
      ],
      id_column =
        "patient_id",
      condition_column =
        "aggregate_state",
      value_column =
        metric,
      condition_high =
        "inactive_lesion",
      condition_low =
        "mature_scar"
    )

    if (nrow(
      inactive_pair
    )) {
      inactive_pair$analysis_set <-
        analysis_set
      inactive_pair$metric <-
        metric

      inactive_scar_pair_rows[[length(
        inactive_scar_pair_rows
      ) + 1L]] <-
        inactive_pair
    }

    activity_comparison <-
      descriptive_two_group(
        value =
          aggregate_subset[[metric]],
        group =
          aggregate_subset$aggregate_state,
        group_a =
          "active_lesion",
        group_b =
          "inactive_lesion"
      )

    activity_comparison$analysis_set <-
      analysis_set
    activity_comparison$metric <-
      metric
    activity_comparison$group_a <-
      "active_lesion"
    activity_comparison$group_b <-
      "inactive_lesion"

    activity_rows[[length(
      activity_rows
    ) + 1L]] <-
      activity_comparison

    keloid_vs_scar_group <- ifelse(
      aggregate_subset$aggregate_state %in%
        c(
          "active_lesion",
          "inactive_lesion"
        ),
      "keloid_lesion",
      "scar_reference"
    )

    lesion_scar_comparison <-
      descriptive_two_group(
        value =
          aggregate_subset[[metric]],
        group =
          keloid_vs_scar_group,
        group_a =
          "keloid_lesion",
        group_b =
          "scar_reference"
      )

    lesion_scar_comparison$analysis_set <-
      analysis_set
    lesion_scar_comparison$metric <-
      metric
    lesion_scar_comparison$group_a <-
      "keloid_lesion"
    lesion_scar_comparison$group_b <-
      "scar_reference"

    keloid_scar_rows[[length(
      keloid_scar_rows
    ) + 1L]] <-
      lesion_scar_comparison
  }
}

region_pair_results <- do.call(
  rbind,
  region_pair_rows
)

inactive_scar_pair_results <- do.call(
  rbind,
  inactive_scar_pair_rows
)

activity_results <- do.call(
  rbind,
  activity_rows
)

keloid_scar_results <- do.call(
  rbind,
  keloid_scar_rows
)

rownames(
  region_pair_results
) <- NULL

rownames(
  inactive_scar_pair_results
) <- NULL

rownames(
  activity_results
) <- NULL

rownames(
  keloid_scar_results
) <- NULL

safe_write_csv(
  region_pair_results,
  file.path(
    REPORT_DIR,
    "10_center_periphery_paired_differences.csv"
  )
)

safe_write_csv(
  inactive_scar_pair_results,
  file.path(
    REPORT_DIR,
    "11_inactive_lesion_vs_paired_mature_scar.csv"
  )
)

safe_write_csv(
  activity_results,
  file.path(
    REPORT_DIR,
    "12_active_vs_inactive_descriptive_comparisons.csv"
  )
)

safe_write_csv(
  keloid_scar_results,
  file.path(
    REPORT_DIR,
    "13_keloid_lesion_vs_scar_reference_descriptive.csv"
  )
)

# ----------------------------- 跨队列综合表 ------------------------------------
cross_cohort_rows <- list()

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

    x <- x[
      x$cell_type ==
        "Fibroblast" &
        x$metric %in%
        c(
          "remodeling_completion",
          "ordinal_position",
          "z_Wound30",
          "wound_activation",
          "z_Skin"
        ),
      ,
      drop = FALSE
    ]

    if (nrow(
      x
    )) {
      cross_cohort_rows[[
        "GSE163973"
      ]] <- data.frame(
        cohort =
          "GSE163973",
        analysis_set =
          "author_annotated_fibroblast",
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
      x$metric %in%
        c(
          "late_remodeling_state",
          "ordinal_wound_state_position",
          "z_Wound30",
          "wound_activation",
          "z_Skin"
        ),
      ,
      drop = FALSE
    ]

    if (nrow(
      x
    )) {
      cross_cohort_rows[[
        "GSE181316"
      ]] <- data.frame(
        cohort =
          "GSE181316",
        analysis_set =
          x$analysis_set,
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

stage7_cross <- keloid_scar_results[
  keloid_scar_results$metric %in%
    c(
      "late_remodeling_state",
      "ordinal_wound_state_position",
      "z_Wound30",
      "wound_activation",
      "z_Skin"
    ),
  ,
  drop = FALSE
]

if (nrow(
  stage7_cross
)) {
  cross_cohort_rows[[
    "GSE220300"
  ]] <- data.frame(
    cohort =
      "GSE220300",
    analysis_set =
      stage7_cross$analysis_set,
    metric =
      stage7_cross$metric,
    mean_difference =
      stage7_cross$mean_difference,
    hedges_g =
      NA_real_,
    cliffs_delta =
      stage7_cross$cliffs_delta,
    stringsAsFactors = FALSE
  )
}

if (length(
  cross_cohort_rows
)) {
  cross_cohort_summary <- do.call(
    rbind,
    cross_cohort_rows
  )

  rownames(
    cross_cohort_summary
  ) <- NULL

  safe_write_csv(
    cross_cohort_summary,
    file.path(
      REPORT_DIR,
      "14_cross_cohort_direction_summary.csv"
    )
  )
} else {
  cross_cohort_summary <- data.frame()
}

# ----------------------------- 决策闸门 ----------------------------------------
get_lesion_scar_row <- function(
  analysis_set,
  metric
) {
  x <- keloid_scar_results[
    keloid_scar_results$analysis_set ==
      analysis_set &
      keloid_scar_results$metric ==
      metric,
    ,
    drop = FALSE
  ]

  if (nrow(
    x
  ) !=
      1L) {
    stop(
      "无法唯一获得病灶-瘢痕比较：",
      analysis_set,
      " / ",
      metric
    )
  }

  x
}

get_inactive_pair <- function(
  analysis_set,
  metric
) {
  inactive_scar_pair_results[
    inactive_scar_pair_results$analysis_set ==
      analysis_set &
      inactive_scar_pair_results$metric ==
      metric,
    ,
    drop = FALSE
  ]
}

late_high <- get_lesion_scar_row(
  "high_specificity",
  "late_remodeling_state"
)

late_broad <- get_lesion_scar_row(
  "broad_sensitivity",
  "late_remodeling_state"
)

ordinal_high <- get_lesion_scar_row(
  "high_specificity",
  "ordinal_wound_state_position"
)

ordinal_broad <- get_lesion_scar_row(
  "broad_sensitivity",
  "ordinal_wound_state_position"
)

inactive_late_high <- get_inactive_pair(
  "high_specificity",
  "late_remodeling_state"
)

inactive_late_broad <- get_inactive_pair(
  "broad_sensitivity",
  "late_remodeling_state"
)

high_lesion_support <- (
  late_high$mean_difference >
    0 &&
    late_high$cliffs_delta >=
    0.5
)

broad_lesion_support <- (
  late_broad$mean_difference >
    0 &&
    late_broad$cliffs_delta >=
    0.33
)

high_ordinal_support <- (
  ordinal_high$mean_difference >
    0 &&
    ordinal_high$cliffs_delta >=
    0.5
)

broad_ordinal_support <- (
  ordinal_broad$mean_difference >
    0 &&
    ordinal_broad$cliffs_delta >=
    0.33
)

paired_mature_scar_support <- (
  nrow(
    inactive_late_high
  ) ==
    2L &&
    all(
      inactive_late_high$paired_difference >
        0
    )
)

paired_mature_scar_broad_support <- (
  nrow(
    inactive_late_broad
  ) ==
    2L &&
    all(
      inactive_late_broad$paired_difference >
        0
    )
)

fibroblast_count_pass <- all(
  qc_summary$n_fibroblast_high >=
    50 &
    qc_summary$n_fibroblast_broad >=
    100
)

if (
  fibroblast_count_pass &&
    high_lesion_support &&
    broad_lesion_support &&
    high_ordinal_support &&
    broad_ordinal_support &&
    paired_mature_scar_support &&
    paired_mature_scar_broad_support
) {
  project_gate <-
    "PASS_ACTIVITY_REGION_AND_PAIRED_SCAR_SUPPORT"

  interpretation <- paste(
    "GSE220300支持晚期重塑状态在活动性与非活动性瘢痕疙瘩病灶中持续存在，",
    "且两名非活动性患者的病灶均高于其配对成熟瘢痕；",
    "高特异性与宽松成纤维细胞定义方向一致。"
  )
} else if (
  fibroblast_count_pass &&
    high_lesion_support &&
    broad_lesion_support
) {
  project_gate <-
    "PASS_GENERAL_LESION_SUPPORT_ONLY"

  interpretation <- paste(
    "瘢痕疙瘩病灶总体高于瘢痕参照，",
    "但患者内成熟瘢痕或序数位置支持不完整。"
  )
} else if (
  fibroblast_count_pass &&
    (
      high_lesion_support ||
        paired_mature_scar_support
    )
) {
  project_gate <-
    "PARTIAL_ACTIVITY_REGION_SUPPORT"

  interpretation <- paste(
    "GSE220300仅提供部分方向支持，",
    "活动度和区域结果必须保持探索性。"
  )
} else {
  project_gate <-
    "NO_GSE220300_SUPPORT"

  interpretation <- paste(
    "GSE220300未支持晚期重塑状态相对瘢痕参照升高；",
    "不能把活动度和区域证据写入主要结论。"
  )
}

decision_criteria <- data.frame(
  criterion = c(
    "fibroblast_count_pass",
    "high_lesion_support",
    "broad_lesion_support",
    "high_ordinal_support",
    "broad_ordinal_support",
    "paired_mature_scar_support",
    "paired_mature_scar_broad_support"
  ),
  passed = c(
    fibroblast_count_pass,
    high_lesion_support,
    broad_lesion_support,
    high_ordinal_support,
    broad_ordinal_support,
    paired_mature_scar_support,
    paired_mature_scar_broad_support
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  decision_criteria,
  file.path(
    REPORT_DIR,
    "15_support_criteria.csv"
  )
)

# ----------------------------- 图形 --------------------------------------------
# 1. 样本级固定状态热图
for (
  analysis_set in c(
    "high_specificity",
    "broad_sensitivity"
  )
) {
  x <- sample_scores[
    sample_scores$analysis_set ==
      analysis_set,
    ,
    drop = FALSE
  ]

  order_ids <- sample_design$sample_id

  x <- x[
    match(
      order_ids,
      x$sample_key
    ),
    ,
    drop = FALSE
  ]

  heat_matrix <- as.matrix(
    x[
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
    x$sample_key,
    x$tissue_state,
    x$region,
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
      paste0(
        "Figure_",
        analysis_set,
        "_stage_score_heatmap.pdf"
      )
    ),
    width = 8,
    height = 9
  )

  par(
    mar = c(
      5,
      12,
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
    main = paste0(
      "GSE220300 wound-state scores: ",
      analysis_set
    )
  )

  axis(
    1,
    at = seq_len(
      ncol(
        heat_matrix
      )
    ),
    labels = colnames(
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
    labels = rownames(
      heat_matrix
    ),
    las = 2,
    cex.axis = 0.7
  )

  box()

  dev.off()
}

# 2. 中心-周边患者内连线
pdf(
  file.path(
    FIG_DIR,
    "Figure_center_periphery_paired_late_remodeling.pdf"
  ),
  width = 10,
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
    5,
    4,
    3,
    1
  )
)

for (
  analysis_set in c(
    "high_specificity",
    "broad_sensitivity"
  )
) {
  for (
    tissue_state in c(
      "active_keloid",
      "inactive_keloid"
    )
  ) {
    x <- sample_scores[
      sample_scores$analysis_set ==
        analysis_set &
        sample_scores$tissue_state ==
        tissue_state,
      ,
      drop = FALSE
    ]

    y_range <- range(
      x$late_remodeling_state,
      finite = TRUE
    )

    y_pad <- max(
      0.1,
      diff(
        y_range
      ) *
        0.15
    )

    plot(
      NA,
      xlim = c(
        0.8,
        2.2
      ),
      ylim = c(
        y_range[1] -
          y_pad,
        y_range[2] +
          y_pad
      ),
      xaxt = "n",
      xlab = "",
      ylab =
        "Late-remodeling state",
      main = paste0(
        analysis_set,
        "\n",
        tissue_state
      )
    )

    axis(
      1,
      at = c(
        1,
        2
      ),
      labels = c(
        "Center",
        "Periphery"
      )
    )

    for (
      patient in unique(
        x$patient_id
      )
    ) {
      xp <- x[
        x$patient_id ==
          patient,
        ,
        drop = FALSE
      ]

      center_value <-
        xp$late_remodeling_state[
          xp$region ==
            "center"
        ]

      peripheral_value <-
        xp$late_remodeling_state[
          xp$region ==
            "periphery"
        ]

      if (
        length(
          center_value
        ) ==
          1L &&
          length(
            peripheral_value
          ) ==
          1L
      ) {
        lines(
          c(
            1,
            2
          ),
          c(
            center_value,
            peripheral_value
          )
        )

        points(
          c(
            1,
            2
          ),
          c(
            center_value,
            peripheral_value
          ),
          pch = 16
        )

        text(
          2,
          peripheral_value,
          labels = patient,
          pos = 4,
          cex = 0.75
        )
      }
    }
  }
}

par(
  op
)

dev.off()

# 3. 非活动病灶与同患者成熟瘢痕
pdf(
  file.path(
    FIG_DIR,
    "Figure_inactive_lesion_vs_paired_mature_scar.pdf"
  ),
  width = 10,
  height = 6
)

op <- par(
  no.readonly = TRUE
)

par(
  mfrow = c(
    1,
    2
  ),
  mar = c(
    5,
    4,
    3,
    1
  )
)

for (
  analysis_set in c(
    "high_specificity",
    "broad_sensitivity"
  )
) {
  x <- aggregate_scores[
    aggregate_scores$analysis_set ==
      analysis_set &
      aggregate_scores$aggregate_state %in%
      c(
        "inactive_lesion",
        "mature_scar"
      ),
    ,
    drop = FALSE
  ]

  y_range <- range(
    x$late_remodeling_state,
    finite = TRUE
  )

  y_pad <- max(
    0.1,
    diff(
      y_range
    ) *
      0.15
  )

  plot(
    NA,
    xlim = c(
      0.8,
      2.2
    ),
    ylim = c(
      y_range[1] -
        y_pad,
      y_range[2] +
        y_pad
    ),
    xaxt = "n",
    xlab = "",
    ylab =
      "Late-remodeling state",
    main =
      analysis_set
  )

  axis(
    1,
    at = c(
      1,
      2
    ),
    labels = c(
      "Mature scar",
      "Inactive lesion"
    )
  )

  for (
    patient in c(
      "I1",
      "I2"
    )
  ) {
    xp <- x[
      x$patient_id ==
        patient,
      ,
      drop = FALSE
    ]

    scar_value <-
      xp$late_remodeling_state[
        xp$aggregate_state ==
          "mature_scar"
      ]

    lesion_value <-
      xp$late_remodeling_state[
        xp$aggregate_state ==
          "inactive_lesion"
      ]

    if (
      length(
        scar_value
      ) ==
        1L &&
        length(
          lesion_value
        ) ==
        1L
    ) {
      lines(
        c(
          1,
          2
        ),
        c(
          scar_value,
          lesion_value
        )
      )

      points(
        c(
          1,
          2
        ),
        c(
          scar_value,
          lesion_value
        ),
        pch = 16
      )

      text(
        2,
        lesion_value,
        labels =
          patient,
        pos = 4,
        cex = 0.8
      )
    }
  }
}

par(
  op
)

dev.off()

# 4. 活动与非活动患者级探索图
pdf(
  file.path(
    FIG_DIR,
    "Figure_active_vs_inactive_patient_level.pdf"
  ),
  width = 10,
  height = 6
)

op <- par(
  no.readonly = TRUE
)

par(
  mfrow = c(
    1,
    2
  ),
  mar = c(
    6,
    4,
    3,
    1
  )
)

for (
  analysis_set in c(
    "high_specificity",
    "broad_sensitivity"
  )
) {
  x <- aggregate_scores[
    aggregate_scores$analysis_set ==
      analysis_set &
      aggregate_scores$aggregate_state %in%
      c(
        "active_lesion",
        "inactive_lesion"
      ),
    ,
    drop = FALSE
  ]

  group_position <- ifelse(
    x$aggregate_state ==
      "inactive_lesion",
    1,
    2
  )

  y_range <- range(
    x$late_remodeling_state,
    finite = TRUE
  )

  y_pad <- max(
    0.1,
    diff(
      y_range
    ) *
      0.15
  )

  plot(
    jitter(
      group_position,
      amount = 0.05
    ),
    x$late_remodeling_state,
    pch = 16,
    xlim = c(
      0.5,
      2.5
    ),
    ylim = c(
      y_range[1] -
        y_pad,
      y_range[2] +
        y_pad
    ),
    xaxt = "n",
    xlab = "",
    ylab =
      "Late-remodeling state",
    main =
      analysis_set
  )

  axis(
    1,
    at = c(
      1,
      2
    ),
    labels = c(
      "Inactive",
      "Active"
    ),
    las = 2
  )

  text(
    jitter(
      group_position,
      amount = 0.05
    ),
    x$late_remodeling_state,
    labels =
      x$patient_id,
    pos = 3,
    cex = 0.8
  )
}

par(
  op
)

dev.off()

# ----------------------------- 保存对象 ----------------------------------------
saveRDS(
  list(
    project_gate =
      project_gate,
    interpretation =
      interpretation,
    sample_design =
      sample_design,
    qc_summary =
      qc_summary,
    marker_summary =
      marker_summary,
    high_specificity =
      high_result,
    broad_sensitivity =
      broad_result,
    sample_scores =
      sample_scores,
    aggregate_scores =
      aggregate_scores,
    region_pair_results =
      region_pair_results,
    inactive_scar_pair_results =
      inactive_scar_pair_results,
    activity_results =
      activity_results,
    keloid_scar_results =
      keloid_scar_results,
    cross_cohort_summary =
      cross_cohort_summary,
    decision_criteria =
      decision_criteria,
    frozen_endpoints = list(
      primary_supportive =
        "late_remodeling_state",
      key_secondary =
        "ordinal_wound_state_position",
      supporting = c(
        "z_Wound30",
        "wound_activation",
        "z_Skin",
        "skin_return_state"
      ),
      exploratory = c(
        "active_vs_inactive",
        "center_vs_periphery",
        "off_trajectory_ratio",
        "early_state_persistence"
      )
    )
  ),
  file = file.path(
    OBJECT_DIR,
    "GSE220300_activity_region_wound_state_support.rds"
  ),
  compress = "gzip"
)

# ----------------------------- 最终结论 ----------------------------------------
decision_lines <- c(
  "第七阶段：GSE220300活动度、区域与成熟瘢痕支持分析结论",
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
  "高特异性Fibroblast病灶相对瘢痕参照：",
  paste0(
    "- late-remodeling mean difference = ",
    round(
      late_high$mean_difference,
      4
    )
  ),
  paste0(
    "- late-remodeling Cliff delta = ",
    round(
      late_high$cliffs_delta,
      3
    )
  ),
  paste0(
    "- ordinal mean difference = ",
    round(
      ordinal_high$mean_difference,
      4
    )
  ),
  paste0(
    "- ordinal Cliff delta = ",
    round(
      ordinal_high$cliffs_delta,
      3
    )
  ),
  "",
  "患者内非活动病灶相对成熟瘢痕：",
  paste0(
    "- high-specificity两名患者是否均升高：",
    paired_mature_scar_support
  ),
  paste0(
    "- broad两名患者是否均升高：",
    paired_mature_scar_broad_support
  ),
  "",
  "解释边界：",
  "- GSE220300独立患者数很少，本阶段不承担正式假设验证。",
  "- active vs inactive只有2例对2例，只报告效应方向，不报告显著性结论。",
  "- center vs periphery每组只有2个患者，只作患者内描述。",
  "- MS01/MS02与I1/I2患者内比较是支持性证据，不等同于纵向随访。",
  "- high-specificity是主要集合，broad是预设敏感性集合。",
  "- 未在GSE220300中重新筛选签名基因或修改成纤维细胞规则。",
  "- ordinal position不是精确术后日龄。"
)

safe_write_lines(
  decision_lines,
  file.path(
    REPORT_DIR,
    "16_STAGE7_DECISION.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    REPORT_DIR,
    "17_SESSION_INFO.txt"
  )
)

safe_write_lines(
  c(
    paste0(
      "STAGE7_COMPLETED=",
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
      fibroblast_count_pass
    ),
    paste0(
      "HIGH_LESION_SUPPORT=",
      high_lesion_support
    ),
    paste0(
      "BROAD_LESION_SUPPORT=",
      broad_lesion_support
    ),
    paste0(
      "HIGH_ORDINAL_SUPPORT=",
      high_ordinal_support
    ),
    paste0(
      "BROAD_ORDINAL_SUPPORT=",
      broad_ordinal_support
    ),
    paste0(
      "PAIRED_MATURE_SCAR_SUPPORT=",
      paired_mature_scar_support
    ),
    paste0(
      "PAIRED_MATURE_SCAR_BROAD_SUPPORT=",
      paired_mature_scar_broad_support
    )
  ),
  file.path(
    REPORT_DIR,
    "STAGE7_COMPLETED.txt"
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
  "第七阶段完成"
)

cat(
  "\n============================================================\n"
)

cat(
  "第七阶段运行完成。\n"
)

cat(
  "本地分析对象保存在：",
  file.path(
    OBJECT_DIR,
    "GSE220300_activity_region_wound_state_support.rds"
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
