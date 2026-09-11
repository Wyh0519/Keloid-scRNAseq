# =============================================================================
# 项目：瘢痕疙瘩“晚期重塑状态持续与伤口状态终止失败”
# 第六阶段B：GSE181316固定伤口状态独立复制分析
#
# 根目录：
#   C:/Users/33652/Desktop/LZZ/1、yssj
#
# 输入：
#   1. 第六阶段A生成的GSE181316 fibroblast-like pseudobulk对象；
#   2. 第三阶段GSE241132正常伤口pseudobulk；
#   3. 第四阶段固定正常伤口状态模型；
#   4. 第五阶段GSE163973发现结果（可选，用于综合复制表）。
#
# 冻结终点：
#   Primary:
#     late_remodeling_state
#     = D30 score - mean(D1 score, D7 score)
#
#   Key secondary:
#     ordinal_wound_state_position
#
#   Supporting:
#     D30-like score
#     wound activation
#     Skin-like score
#
#   Exploratory:
#     off-trajectory ratio
#     early-state persistence
#
# 关键设计：
#   - high-specificity fibroblast-like为主要集合；
#   - broad fibroblast-like为预设敏感性集合；
#   - keloid_3L和keloid_3R先在原始计数层面合并为患者K3；
#   - 统计单位为患者：3例keloid vs 3例normal scar；
#   - healthy skin仅作描述；
#   - 不在GSE181316中重新选择基因、修改签名或调整细胞规则；
#   - 不将ordinal position解释为精确伤口日龄。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANG = "en")

ROOT_DIR <- "C:/Users/33652/Desktop/LZZ/1、yssj"

GSE181316_RDS <- file.path(
  ROOT_DIR,
  "06A_STAGE6A_GSE181316_CELL_CALLING",
  "objects",
  "GSE181316_fibroblast_like_pseudobulk_list.rds"
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

DISCOVERY_RDS <- file.path(
  ROOT_DIR,
  "05_STAGE5_GSE163973_MAPPING",
  "objects",
  "GSE163973_fixed_wound_state_mapping.rds"
)

STAGE_DIR <- file.path(
  ROOT_DIR,
  "06B_STAGE6B_GSE181316_REPLICATION"
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

PACKAGE_ZIP <- file.path(
  ROOT_DIR,
  "第六阶段B_GSE181316晚期重塑状态独立复制检查包.zip"
)

PACKAGE_TARGZ <- file.path(
  ROOT_DIR,
  "第六阶段B_GSE181316晚期重塑状态独立复制检查包.tar.gz"
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

LOG_FILE <- file.path(
  REPORT_DIR,
  "00_stage6B_log.txt"
)

WARNING_FILE <- file.path(
  REPORT_DIR,
  "00_stage6B_warnings.txt"
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

# ----------------------------- 基础函数 ----------------------------------------
canonicalize_named_vector <- function(x) {
  if (is.null(names(x))) {
    stop(
      "pseudobulk向量缺少基因名。"
    )
  }

  gene <- toupper(
    trimws(
      as.character(
        names(x)
      )
    )
  )

  value <- as.numeric(
    x
  )

  valid <- (
    !is.na(gene) &
      nzchar(gene) &
      is.finite(value)
  )

  if (!any(valid)) {
    stop(
      "pseudobulk向量没有有效基因。"
    )
  }

  summed <- rowsum(
    matrix(
      value[valid],
      ncol = 1
    ),
    group = gene[valid],
    reorder = FALSE
  )

  out <- as.numeric(
    summed[, 1]
  )

  names(out) <- rownames(
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
    !is.na(gene) &
      nzchar(gene)
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
  required_names = NULL
) {
  if (!length(vector_list)) {
    stop(
      "pseudobulk列表为空。"
    )
  }

  vector_list <- lapply(
    vector_list,
    canonicalize_named_vector
  )

  if (!is.null(required_names)) {
    missing_samples <- setdiff(
      required_names,
      names(vector_list)
    )

    if (length(missing_samples)) {
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
  }

  common_genes <- Reduce(
    intersect,
    lapply(
      vector_list,
      names
    )
  )

  if (length(common_genes) < 10000L) {
    warn_msg(
      "GSE181316样本间共同基因少于10,000：",
      length(common_genes)
    )
  }

  mat <- do.call(
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

  rownames(mat) <- common_genes
  colnames(mat) <- names(
    vector_list
  )
  storage.mode(mat) <- "double"

  mat
}

aggregate_samples_to_patients <- function(
  sample_counts,
  sample_design
) {
  required_samples <- sample_design$sample_id

  missing_samples <- setdiff(
    required_samples,
    colnames(
      sample_counts
    )
  )

  if (length(missing_samples)) {
    stop(
      "计数矩阵缺少样本：",
      paste(
        missing_samples,
        collapse = ", "
      )
    )
  }

  patient_order <- c(
    "H7",
    "K1",
    "K2",
    "K3",
    "S1",
    "S2",
    "S3"
  )

  patient_counts <- do.call(
    cbind,
    lapply(
      patient_order,
      function(patient) {
        sample_ids <- sample_design$sample_id[
          sample_design$patient_id == patient
        ]

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
    patient_counts
  ) <- patient_order

  patient_info <- do.call(
    rbind,
    lapply(
      patient_order,
      function(patient) {
        x <- sample_design[
          sample_design$patient_id == patient,
          ,
          drop = FALSE
        ]

        data.frame(
          patient_id = patient,
          group = unique(
            x$group
          )[1],
          component_samples = paste(
            x$sample_id,
            collapse = "+"
          ),
          n_component_samples = nrow(x),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  rownames(
    patient_info
  ) <- NULL

  list(
    counts = patient_counts,
    patient_info = patient_info
  )
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
      "正常伤口pseudobulk列名无法全部匹配sample_info。"
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
      "正常伤口pseudobulk设计顺序匹配失败。"
    )
  }

  out
}

get_stable_fibroblast_signatures <- function(
  model
) {
  stability <- model$signature_stability

  required_cols <- c(
    "cell_type",
    "stage",
    "gene",
    "stable_all_3_folds"
  )

  if (!all(
    required_cols %in%
      names(stability)
  )) {
    stop(
      "第四阶段模型缺少稳定签名字段。"
    )
  }

  x <- stability[
    stability$cell_type == "Fibroblast" &
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
                  x$stage == stage
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

  if (any(n_genes < 20L)) {
    stop(
      "稳定Fibroblast签名不足：",
      paste(
        names(n_genes),
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
      "存在文库总计数≤0的pseudobulk。"
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

  if (!all(
    stages %in%
      names(signatures)
  )) {
    stop(
      "签名缺少四个正常伤口阶段。"
    )
  }

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

    names(ranks) <- rownames(
      expression_matrix
    )

    for (
      stage in stages
    ) {
      genes <- intersect(
        signatures[[stage]],
        names(ranks)
      )

      if (length(genes) < 15L) {
        stop(
          "阶段",
          stage,
          "在共同基因背景中少于15个稳定签名基因：",
          length(genes)
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
      scale_value == 0
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
    reference_z = reference_z,
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
            reference_design$condition == stage,
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
      reference_design$condition == stage &
        reference_design$donor != donor
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
        sum(weights) == 0
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
      nearest_stage = nearest_stage,
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
      mixed_state_entropy = entropy,
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

  if (anyNA(idx)) {
    stop(
      "映射结果无法与z-score匹配。"
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

exact_patient_comparison <- function(
  value,
  group,
  group_a = "keloid",
  group_b = "normal_scar"
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
    group == group_a
  ]

  y <- value[
    group == group_b
  ]

  n1 <- length(x)
  n2 <- length(y)

  if (
    n1 != 3L ||
      n2 != 3L
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
        n_a_above_median_b = NA_integer_,
        complete_separation = NA
      )
    )
  }

  observed_difference <- mean(x) -
    mean(y)

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
      (n1 - 1) *
        stats::var(x) +
        (n2 - 1) *
        stats::var(y)
    ) /
      degrees_freedom
  )

  hedges_g <- NA_real_

  if (
    is.finite(
      pooled_sd
    ) &&
      pooled_sd > 0
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
    pairwise_difference > 0
  ) -
    mean(
      pairwise_difference < 0
    )

  list(
    n_a = n1,
    n_b = n2,
    mean_a = mean(x),
    mean_b = mean(y),
    median_a = stats::median(x),
    median_b = stats::median(y),
    mean_difference =
      observed_difference,
    hodges_lehmann_shift =
      stats::median(
        as.vector(
          pairwise_difference
        )
      ),
    hedges_g = hedges_g,
    cliffs_delta = cliffs_delta,
    exact_two_sided_p =
      exact_two_sided_p,
    exact_one_sided_greater_p =
      exact_one_sided_greater_p,
    n_a_above_median_b =
      sum(
        x >
          stats::median(y)
      ),
    complete_separation =
      min(x) >
      max(y)
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
        ) >= 0.8
    ) ||
      (
        is.finite(
          cliffs_delta
        ) &&
          abs(
            cliffs_delta
          ) >= 0.56
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
        ) >= 0.5
    ) ||
      (
        is.finite(
          cliffs_delta
        ) &&
          abs(
            cliffs_delta
          ) >= 0.33
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
  sample_vector_list,
  sample_design,
  reference_counts,
  reference_design,
  signatures
) {
  log_msg(
    "分析fibroblast集合：",
    set_name
  )

  required_samples <- sample_design$sample_id

  sample_counts <- merge_named_vector_list(
    sample_vector_list,
    required_names = required_samples
  )

  patient_result <- aggregate_samples_to_patients(
    sample_counts,
    sample_design
  )

  patient_counts <- patient_result$counts
  patient_info <- patient_result$patient_info

  common_genes <- Reduce(
    intersect,
    list(
      rownames(
        reference_counts
      ),
      rownames(
        sample_counts
      )
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
    ) < 5000L
  ) {
    stop(
      set_name,
      "与正常伤口参考的共同背景基因少于5000：",
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
    signature_coverage < 15L
  )) {
    stop(
      set_name,
      "至少一个阶段的稳定签名共同基因少于15：",
      paste(
        names(
          signature_coverage
        ),
        signature_coverage,
        sep = "=",
        collapse = ", "
      )
    )
  }

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

  patient_expression <- log_cpm(
    patient_counts[
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

  patient_scores <- score_samples_by_rank(
    patient_expression,
    signatures
  )

  standardized_sample <-
    standardize_by_reference(
      reference_scores,
      sample_scores
    )

  standardized_patient <-
    standardize_by_reference(
      reference_scores,
      patient_scores
    )

  reference_z <-
    standardized_sample$reference_z

  sample_z <-
    standardized_sample$new_z

  patient_z <-
    standardized_patient$new_z

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
    sample_design[
      ,
      c(
        "sample_id",
        "group",
        "patient_id",
        "independent_patient"
      )
    ],
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

  patient_mapping <- map_to_reference(
    patient_z,
    centroids,
    reference_distance_limit
  )

  patient_mapping <- add_fixed_contrasts(
    patient_mapping,
    patient_z
  )

  patient_mapping$analysis_set <-
    set_name

  patient_mapping <- merge(
    patient_mapping,
    patient_info,
    by.x = "sample_key",
    by.y = "patient_id",
    all.x = TRUE,
    sort = FALSE
  )

  patient_mapping <- patient_mapping[
    match(
      rownames(
        patient_z
      ),
      patient_mapping$sample_key
    ),
    ,
    drop = FALSE
  ]

  patient_mapping$patient_id <-
    patient_mapping$sample_key

  comparison_metrics <- c(
    "late_remodeling_state",
    "ordinal_wound_state_position",
    "z_Wound30",
    "wound_activation",
    "z_Skin",
    "skin_return_state",
    "off_trajectory_ratio",
    "early_state_persistence",
    "mixed_state_entropy"
  )

  comparison_rows <- list()

  for (
    metric in comparison_metrics
  ) {
    test_result <- exact_patient_comparison(
      patient_mapping[[metric]],
      patient_mapping$group
    )

    comparison_rows[[length(comparison_rows) + 1L]] <-
      data.frame(
        analysis_set = set_name,
        metric = metric,
        endpoint_role = if (
          metric == "late_remodeling_state"
        ) {
          "primary"
        } else if (
          metric ==
            "ordinal_wound_state_position"
        ) {
          "key_secondary"
        } else if (
          metric %in%
            c(
              "z_Wound30",
              "wound_activation",
              "z_Skin",
              "skin_return_state"
            )
        ) {
          "supporting"
        } else {
          "exploratory"
        },
        n_keloid = test_result$n_a,
        n_normal_scar = test_result$n_b,
        mean_keloid = test_result$mean_a,
        mean_normal_scar = test_result$mean_b,
        median_keloid = test_result$median_a,
        median_normal_scar = test_result$median_b,
        mean_difference =
          test_result$mean_difference,
        hodges_lehmann_shift =
          test_result$hodges_lehmann_shift,
        hedges_g =
          test_result$hedges_g,
        cliffs_delta =
          test_result$cliffs_delta,
        exact_two_sided_p =
          test_result$exact_two_sided_p,
        exact_one_sided_greater_p =
          test_result$exact_one_sided_greater_p,
        n_keloid_above_normal_median =
          test_result$n_a_above_median_b,
        complete_separation =
          test_result$complete_separation,
        effect_strength = effect_strength(
          test_result$hedges_g,
          test_result$cliffs_delta
        ),
        stringsAsFactors = FALSE
      )
  }

  comparison <- do.call(
    rbind,
    comparison_rows
  )

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
    patient_counts = patient_counts,
    patient_info = patient_info,
    sample_mapping = sample_mapping,
    patient_mapping = patient_mapping,
    comparison = comparison,
    coverage = coverage,
    calibration = list(
      common_genes = common_genes,
      signatures = signatures,
      reference_center =
        standardized_sample$center,
      reference_scale =
        standardized_sample$scale,
      centroids = centroids,
      reference_lodo_distance =
        reference_lodo_distance,
      reference_distance_limit =
        reference_distance_limit
    )
  )
}

# ----------------------------- 输入检查 ----------------------------------------
log_msg(
  "第六阶段B开始"
)

required_files <- c(
  GSE181316_RDS,
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

gse181316_object <- readRDS(
  GSE181316_RDS
)

reference_pb_object <- readRDS(
  REFERENCE_PB_RDS
)

reference_model <- readRDS(
  REFERENCE_MODEL_RDS
)

required_181316_names <- c(
  "high_specificity",
  "broad_sensitivity",
  "sample_design",
  "cell_call_summary"
)

missing_181316_names <- setdiff(
  required_181316_names,
  names(
    gse181316_object
  )
)

if (length(
  missing_181316_names
)) {
  stop(
    "第六阶段A对象缺少：",
    paste(
      missing_181316_names,
      collapse = ", "
    )
  )
}

sample_design <-
  gse181316_object$sample_design

required_sample_ids <- c(
  "skin_7",
  "keloid_1",
  "keloid_2",
  "keloid_3L",
  "keloid_3R",
  "scar_1",
  "scar_2",
  "scar_3"
)

if (!setequal(
  sample_design$sample_id,
  required_sample_ids
)) {
  stop(
    "第六阶段A样本设计与预期不一致。"
  )
}

reference_counts_raw <-
  reference_pb_object$counts

reference_gene_symbol <-
  reference_pb_object$gene_symbol

reference_sample_info <-
  reference_pb_object$sample_info

reference_design_all <- get_reference_design(
  reference_counts_raw,
  reference_sample_info
)

reference_counts_all <-
  canonicalize_reference_counts(
    reference_counts_raw,
    reference_gene_symbol
  )

reference_fibroblast_idx <- which(
  reference_design_all$main_cell_type ==
    "Fibroblast"
)

if (length(
  reference_fibroblast_idx
) != 12L) {
  stop(
    "正常伤口Fibroblast pseudobulk不是12个，实际：",
    length(
      reference_fibroblast_idx
    )
  )
}

reference_counts <-
  reference_counts_all[
    ,
    reference_fibroblast_idx,
    drop = FALSE
  ]

reference_design <-
  reference_design_all[
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
    "01_fixed_signature_summary.csv"
  )
)

# ----------------------------- 两套预设分析 ------------------------------------
high_result <- analyse_fibroblast_set(
  set_name = "high_specificity",
  sample_vector_list =
    gse181316_object$high_specificity,
  sample_design = sample_design,
  reference_counts = reference_counts,
  reference_design = reference_design,
  signatures = signatures
)

broad_result <- analyse_fibroblast_set(
  set_name = "broad_sensitivity",
  sample_vector_list =
    gse181316_object$broad_sensitivity,
  sample_design = sample_design,
  reference_counts = reference_counts,
  reference_design = reference_design,
  signatures = signatures
)

sample_scores <- rbind(
  high_result$sample_mapping,
  broad_result$sample_mapping
)

patient_scores <- rbind(
  high_result$patient_mapping,
  broad_result$patient_mapping
)

comparisons <- rbind(
  high_result$comparison,
  broad_result$comparison
)

coverage <- rbind(
  high_result$coverage,
  broad_result$coverage
)

rownames(
  sample_scores
) <- NULL

rownames(
  patient_scores
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
    "02_signature_gene_coverage.csv"
  )
)

safe_write_csv(
  sample_scores,
  file.path(
    REPORT_DIR,
    "03_sample_level_scores_before_K3_merge.csv"
  )
)

safe_write_csv(
  patient_scores,
  file.path(
    REPORT_DIR,
    "04_patient_level_scores_after_K3_merge.csv"
  )
)

safe_write_csv(
  comparisons,
  file.path(
    REPORT_DIR,
    "05_patient_level_group_comparisons.csv"
  )
)

# ----------------------------- K3左右一致性 ------------------------------------
k3_concordance <- sample_scores[
  sample_scores$sample_key %in%
    c(
      "keloid_3L",
      "keloid_3R"
    ),
  c(
    "analysis_set",
    "sample_key",
    "late_remodeling_state",
    "ordinal_wound_state_position",
    "z_Wound30",
    "z_Skin",
    "off_trajectory_ratio"
  ),
  drop = FALSE
]

safe_write_csv(
  k3_concordance,
  file.path(
    REPORT_DIR,
    "06_K3_left_right_concordance.csv"
  )
)

# ----------------------------- 发现-验证综合 -----------------------------------
integrated_replication <- data.frame()

if (file.exists(
  DISCOVERY_RDS
)) {
  discovery <- readRDS(
    DISCOVERY_RDS
  )

  if (
    "comparisons" %in%
      names(discovery)
  ) {
    discovery_comparison <-
      discovery$comparisons

    discovery_map <- data.frame(
      discovery_metric = c(
        "remodeling_completion",
        "ordinal_position",
        "z_Wound30",
        "wound_activation",
        "z_Skin",
        "off_trajectory_ratio",
        "early_state_persistence"
      ),
      validation_metric = c(
        "late_remodeling_state",
        "ordinal_wound_state_position",
        "z_Wound30",
        "wound_activation",
        "z_Skin",
        "off_trajectory_ratio",
        "early_state_persistence"
      ),
      stringsAsFactors = FALSE
    )

    discovery_fib <- discovery_comparison[
      discovery_comparison$cell_type ==
        "Fibroblast",
      ,
      drop = FALSE
    ]

    integrated_rows <- list()

    for (
      i in seq_len(
        nrow(
          discovery_map
        )
      )
    ) {
      discovery_row <- discovery_fib[
        discovery_fib$metric ==
          discovery_map$discovery_metric[i],
        ,
        drop = FALSE
      ]

      validation_rows <- comparisons[
        comparisons$metric ==
          discovery_map$validation_metric[i],
        ,
        drop = FALSE
      ]

      if (
        nrow(
          discovery_row
        ) == 1L &&
          nrow(
            validation_rows
          ) > 0L
      ) {
        for (
          j in seq_len(
            nrow(
              validation_rows
            )
          )
        ) {
          integrated_rows[[length(
            integrated_rows
          ) + 1L]] <-
            data.frame(
              metric =
                discovery_map$validation_metric[i],
              validation_analysis_set =
                validation_rows$analysis_set[j],
              discovery_mean_difference =
                discovery_row$mean_difference,
              discovery_hedges_g =
                discovery_row$hedges_g,
              discovery_cliffs_delta =
                discovery_row$cliffs_delta,
              validation_mean_difference =
                validation_rows$mean_difference[j],
              validation_hedges_g =
                validation_rows$hedges_g[j],
              validation_cliffs_delta =
                validation_rows$cliffs_delta[j],
              same_direction =
                sign(
                  discovery_row$mean_difference
                ) ==
                sign(
                  validation_rows$mean_difference[j]
                ),
              stringsAsFactors = FALSE
            )
        }
      }
    }

    if (length(
      integrated_rows
    )) {
      integrated_replication <- do.call(
        rbind,
        integrated_rows
      )

      rownames(
        integrated_replication
      ) <- NULL

      safe_write_csv(
        integrated_replication,
        file.path(
          REPORT_DIR,
          "07_discovery_validation_replication_table.csv"
        )
      )
    }
  }
} else {
  warn_msg(
    "未找到第五阶段RDS，不生成发现-验证综合表；不影响第六阶段B主要分析。"
  )
}

# ----------------------------- 复制判定 ----------------------------------------
get_comparison_row <- function(
  analysis_set,
  metric
) {
  x <- comparisons[
    comparisons$analysis_set ==
      analysis_set &
      comparisons$metric ==
      metric,
    ,
    drop = FALSE
  ]

  if (nrow(
    x
  ) != 1L) {
    stop(
      "无法唯一获得比较结果：",
      analysis_set,
      " / ",
      metric
    )
  }

  x
}

primary_high <- get_comparison_row(
  "high_specificity",
  "late_remodeling_state"
)

primary_broad <- get_comparison_row(
  "broad_sensitivity",
  "late_remodeling_state"
)

ordinal_high <- get_comparison_row(
  "high_specificity",
  "ordinal_wound_state_position"
)

ordinal_broad <- get_comparison_row(
  "broad_sensitivity",
  "ordinal_wound_state_position"
)

primary_high_positive <- (
  primary_high$mean_difference > 0 &&
    primary_high$cliffs_delta >= 0.56 &&
    primary_high$n_keloid_above_normal_median >= 2L &&
    primary_high$exact_one_sided_greater_p <= 0.10
)

primary_broad_concordant <- (
  primary_broad$mean_difference > 0 &&
    primary_broad$cliffs_delta >= 0.33
)

ordinal_high_positive <- (
  ordinal_high$mean_difference > 0 &&
    ordinal_high$cliffs_delta >= 0.56
)

ordinal_broad_concordant <- (
  ordinal_broad$mean_difference > 0 &&
    ordinal_broad$cliffs_delta >= 0.33
)

k3_high <- k3_concordance[
  k3_concordance$analysis_set ==
    "high_specificity",
  ,
  drop = FALSE
]

normal_sample_median_high <- stats::median(
  sample_scores$late_remodeling_state[
    sample_scores$analysis_set ==
      "high_specificity" &
      sample_scores$group ==
      "normal_scar"
  ]
)

k3_left_right_concordant <- (
  nrow(
    k3_high
  ) == 2L &&
    all(
      k3_high$late_remodeling_state >
        normal_sample_median_high
    )
)

if (
  primary_high_positive &&
    primary_broad_concordant &&
    ordinal_high_positive &&
    ordinal_broad_concordant &&
    k3_left_right_concordant
) {
  project_gate <-
    "PASS_STRONG_INDEPENDENT_REPLICATION"

  interpretation <- paste(
    "GSE181316在高特异性与宽松成纤维细胞集合中均复制了",
    "瘢痕疙瘩晚期重塑状态升高及序数伤口状态后移；",
    "K3左右病灶方向一致。"
  )
} else if (
  primary_high_positive &&
    primary_broad_concordant
) {
  project_gate <-
    "PASS_PRIMARY_REPLICATION"

  interpretation <- paste(
    "主要late-remodeling终点得到复制，",
    "但序数位置或K3左右一致性未达到全部强标准。"
  )
} else if (
  primary_high$mean_difference > 0 &&
    (
      primary_high$cliffs_delta >= 0.33 ||
        primary_high$hedges_g >= 0.5
    )
) {
  project_gate <-
    "PARTIAL_REPLICATION"

  interpretation <- paste(
    "主要终点方向一致且至少为中等效应，",
    "但未达到预设独立复制标准。"
  )
} else {
  project_gate <-
    "FAIL_INDEPENDENT_REPLICATION"

  interpretation <- paste(
    "GSE181316未复制GSE163973中的晚期重塑状态升高；",
    "不得将该机制作为多队列稳定结论。"
  )
}

replication_decision <- data.frame(
  criterion = c(
    "primary_high_positive",
    "primary_broad_concordant",
    "ordinal_high_positive",
    "ordinal_broad_concordant",
    "K3_left_right_concordant"
  ),
  passed = c(
    primary_high_positive,
    primary_broad_concordant,
    ordinal_high_positive,
    ordinal_broad_concordant,
    k3_left_right_concordant
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  replication_decision,
  file.path(
    REPORT_DIR,
    "08_replication_criteria.csv"
  )
)

# ----------------------------- 图形 --------------------------------------------
plot_patient_metric <- function(
  data,
  metric,
  title,
  ylab
) {
  analysis_sets <- c(
    "high_specificity",
    "broad_sensitivity"
  )

  for (
    analysis_set in analysis_sets
  ) {
    x <- data[
      data$analysis_set ==
        analysis_set &
        data$group %in%
        c(
          "normal_scar",
          "keloid"
        ),
      ,
      drop = FALSE
    ]

    x_position <- ifelse(
      x$group ==
        "normal_scar",
      1,
      2
    )

    plot(
      jitter(
        x_position,
        amount = 0.05
      ),
      x[[metric]],
      pch = ifelse(
        x$group ==
          "normal_scar",
        1,
        16
      ),
      xaxt = "n",
      xlim = c(
        0.5,
        2.5
      ),
      xlab = "",
      ylab = ylab,
      main = paste0(
        title,
        "\n",
        analysis_set
      )
    )

    axis(
      1,
      at = c(
        1,
        2
      ),
      labels = c(
        "Normal scar",
        "Keloid"
      ),
      las = 2
    )

    text(
      jitter(
        x_position,
        amount = 0.05
      ),
      x[[metric]],
      labels = x$patient_id,
      pos = 3,
      cex = 0.75
    )

    segments(
      0.85,
      stats::median(
        x[[metric]][
          x$group ==
            "normal_scar"
        ]
      ),
      1.15,
      stats::median(
        x[[metric]][
          x$group ==
            "normal_scar"
        ]
      ),
      lwd = 2
    )

    segments(
      1.85,
      stats::median(
        x[[metric]][
          x$group ==
            "keloid"
        ]
      ),
      2.15,
      stats::median(
        x[[metric]][
          x$group ==
            "keloid"
        ]
      ),
      lwd = 2
    )
  }
}

pdf(
  file.path(
    FIG_DIR,
    "Figure_patient_level_replication_metrics.pdf"
  ),
  width = 11,
  height = 8,
  onefile = TRUE
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
    4,
    1
  )
)

plot_patient_metric(
  patient_scores,
  "late_remodeling_state",
  "Primary endpoint",
  "Late-remodeling state"
)

plot_patient_metric(
  patient_scores,
  "ordinal_wound_state_position",
  "Key secondary endpoint",
  "Ordinal wound-state position"
)

par(op)
dev.off()

# 高特异性四阶段热图
high_patient <- patient_scores[
  patient_scores$analysis_set ==
    "high_specificity",
  ,
  drop = FALSE
]

high_patient <- high_patient[
  match(
    c(
      "H7",
      "S1",
      "S2",
      "S3",
      "K1",
      "K2",
      "K3"
    ),
    high_patient$patient_id
  ),
  ,
  drop = FALSE
]

heat_matrix <- as.matrix(
  high_patient[
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
  high_patient$patient_id,
  high_patient$group,
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
    "Figure_high_specificity_stage_score_heatmap.pdf"
  ),
  width = 7,
  height = 7
)

par(
  mar = c(
    5,
    8,
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
  main = "GSE181316 fixed wound-state scores"
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
  cex.axis = 0.8
)

box()
dev.off()

# K3左右病灶一致性图
pdf(
  file.path(
    FIG_DIR,
    "Figure_K3_left_right_concordance.pdf"
  ),
  width = 8,
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
  x <- sample_scores[
    sample_scores$analysis_set ==
      analysis_set,
    ,
    drop = FALSE
  ]

  normal_values <- x$late_remodeling_state[
    x$group ==
      "normal_scar"
  ]

  k3_values <- x[
    x$sample_key %in%
      c(
        "keloid_3L",
        "keloid_3R"
      ),
    ,
    drop = FALSE
  ]

  plot(
    seq_along(
      normal_values
    ),
    normal_values,
    pch = 1,
    xlim = c(
      0.5,
      5.5
    ),
    xaxt = "n",
    xlab = "",
    ylab = "Late-remodeling state",
    main = analysis_set
  )

  points(
    c(
      4,
      5
    ),
    k3_values$late_remodeling_state,
    pch = 16
  )

  axis(
    1,
    at = 1:5,
    labels = c(
      "S1",
      "S2",
      "S3",
      "K3L",
      "K3R"
    ),
    las = 2
  )

  abline(
    h = stats::median(
      normal_values
    ),
    lty = 2
  )
}

par(op)
dev.off()

# 发现-验证效应方向图
if (
  nrow(
    integrated_replication
  ) > 0L
) {
  late_integrated <- integrated_replication[
    integrated_replication$metric ==
      "late_remodeling_state",
    ,
    drop = FALSE
  ]

  if (
    nrow(
      late_integrated
    ) > 0L
  ) {
    pdf(
      file.path(
        FIG_DIR,
        "Figure_discovery_validation_effects.pdf"
      ),
      width = 8,
      height = 6
    )

    effect_matrix <- rbind(
      Discovery_GSE163973 =
        late_integrated$discovery_hedges_g[1],
      Validation_high =
        late_integrated$validation_hedges_g[
          late_integrated$validation_analysis_set ==
            "high_specificity"
        ][1],
      Validation_broad =
        late_integrated$validation_hedges_g[
          late_integrated$validation_analysis_set ==
            "broad_sensitivity"
        ][1]
    )

    barplot(
      effect_matrix,
      ylab = "Hedges g",
      main = "Late-remodeling state: discovery and validation",
      las = 2
    )

    abline(
      h = 0,
      lty = 2
    )

    dev.off()
  }
}

# ----------------------------- 保存对象 ----------------------------------------
saveRDS(
  list(
    project_gate = project_gate,
    interpretation = interpretation,
    high_specificity = high_result,
    broad_sensitivity = broad_result,
    patient_scores = patient_scores,
    sample_scores = sample_scores,
    comparisons = comparisons,
    replication_criteria =
      replication_decision,
    integrated_replication =
      integrated_replication,
    frozen_endpoints = list(
      primary =
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
        "off_trajectory_ratio",
        "early_state_persistence"
      )
    )
  ),
  file = file.path(
    OBJECT_DIR,
    "GSE181316_fixed_wound_state_replication.rds"
  ),
  compress = "gzip"
)

# ----------------------------- 最终决策 ----------------------------------------
decision_lines <- c(
  "第六阶段B：GSE181316固定伤口状态独立复制结论",
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
  "高特异性Fibroblast主要终点：",
  paste0(
    "- mean difference = ",
    round(
      primary_high$mean_difference,
      4
    )
  ),
  paste0(
    "- Hedges g = ",
    round(
      primary_high$hedges_g,
      3
    )
  ),
  paste0(
    "- Cliff delta = ",
    round(
      primary_high$cliffs_delta,
      3
    )
  ),
  paste0(
    "- exact two-sided P = ",
    round(
      primary_high$exact_two_sided_p,
      3
    )
  ),
  paste0(
    "- exact directional one-sided P = ",
    round(
      primary_high$exact_one_sided_greater_p,
      3
    )
  ),
  paste0(
    "- complete separation = ",
    primary_high$complete_separation
  ),
  "",
  "高特异性Fibroblast关键次要终点：",
  paste0(
    "- ordinal mean difference = ",
    round(
      ordinal_high$mean_difference,
      4
    )
  ),
  paste0(
    "- ordinal Hedges g = ",
    round(
      ordinal_high$hedges_g,
      3
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
  "敏感性与同一患者一致性：",
  paste0(
    "- broad primary direction concordant = ",
    primary_broad_concordant
  ),
  paste0(
    "- broad ordinal direction concordant = ",
    ordinal_broad_concordant
  ),
  paste0(
    "- K3 left/right concordant = ",
    k3_left_right_concordant
  ),
  "",
  "固定解释边界：",
  "- 统计单位为3例keloid与3例normal scar患者。",
  "- K3左右病灶已在计数层面合并后进入患者级比较。",
  "- 双侧exact P最小通常为0.10；单侧P仅因验证方向在运行前已冻结而作为支持。",
  "- high-specificity为主要集合，broad集合只作为敏感性分析。",
  "- ordinal position表示相对正常伤口状态，不是精确术后日龄。",
  "- healthy skin H7只作描述，不进入3 vs 3检验。",
  "- 未使用GSE181316重新筛选签名基因或调整成纤维细胞规则。"
)

safe_write_lines(
  decision_lines,
  file.path(
    REPORT_DIR,
    "09_STAGE6B_DECISION.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    REPORT_DIR,
    "10_SESSION_INFO.txt"
  )
)

safe_write_lines(
  c(
    paste0(
      "STAGE6B_COMPLETED=",
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
      "PRIMARY_HIGH_POSITIVE=",
      primary_high_positive
    ),
    paste0(
      "PRIMARY_BROAD_CONCORDANT=",
      primary_broad_concordant
    ),
    paste0(
      "ORDINAL_HIGH_POSITIVE=",
      ordinal_high_positive
    ),
    paste0(
      "ORDINAL_BROAD_CONCORDANT=",
      ordinal_broad_concordant
    ),
    paste0(
      "K3_LEFT_RIGHT_CONCORDANT=",
      k3_left_right_concordant
    )
  ),
  file.path(
    REPORT_DIR,
    "STAGE6B_COMPLETED.txt"
  )
)

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
    )$size > 0
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
  "第六阶段B完成"
)

cat(
  "\n============================================================\n"
)

cat(
  "第六阶段B运行完成。\n"
)

cat(
  "独立复制对象保存在：",
  file.path(
    OBJECT_DIR,
    "GSE181316_fixed_wound_state_replication.rds"
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
