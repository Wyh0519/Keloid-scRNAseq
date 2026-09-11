# =============================================================================
# 项目：瘢痕疙瘩“晚期重塑状态持续与伤口状态终止失败”
# 第十三阶段：late-remodeling 固定签名稳健性敏感性分析
#
# 目的：
#   1. 不改变第四阶段冻结的 Skin / Wound1 / Wound7 / Wound30 基本定义；
#   2. 对 Wound30 稳定签名做 leave-one-gene-out (LOGO)；
#   3. 对 Wound30 稳定签名做随机 80% 基因子采样；
#   4. 每次都重新计算正常伤口参考分数、参考标准化和疾病队列 late-remodeling score；
#   5. 分别在 GSE163973（发现）和 GSE181316（独立复制）中评估方向、效应量与完全分离；
#   6. 对两个 3-vs-3 核心队列做 20×20=400 个联合标签分配的分层精确整合；
#   7. 判断主结论是否由单个 Wound30 基因或少数基因组合驱动。
#
# 说明：
#   - 本代码直接衔接你现有 Stage3/4/5/6B 的 RDS 对象；
#   - 不重新筛选疾病相关基因；
#   - 不以单队列 exact P=0.10 作为“失败”，因为 3 vs 3 的双侧精确置换最小非零P即0.10；
#   - 重点看：两个队列是否同方向、是否仍完全分离、整合效应是否稳定；
#   - 随机子采样默认 1000 次，可按运行速度改为 500 或 2000。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANG = "en")

# ----------------------------- 用户配置 ----------------------------------------
ROOT_DIR <- "C:/Users/33652/Desktop/LZZ/1、yssj"

N_RANDOM <- 1000L
KEEP_FRACTION <- 0.80
RANDOM_SEED <- 20260908L

# 若只想先快速测试代码，可临时设：
# N_RANDOM <- 100L

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

VALIDATION_RDS <- file.path(
  ROOT_DIR,
  "06B_STAGE6B_GSE181316_REPLICATION",
  "objects",
  "GSE181316_fixed_wound_state_replication.rds"
)

STAGE_DIR <- file.path(
  ROOT_DIR,
  "13_STAGE13_SIGNATURE_ROBUSTNESS"
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
  "第十三阶段_Wound30固定签名稳健性敏感性分析检查包.zip"
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
  "00_stage13_log.txt"
)

WARNING_FILE <- file.path(
  REPORT_DIR,
  "00_stage13_warnings.txt"
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

# ----------------------------- 输入检查 ----------------------------------------
required_files <- c(
  REFERENCE_PB_RDS,
  REFERENCE_MODEL_RDS,
  DISCOVERY_RDS,
  VALIDATION_RDS
)

missing_files <- required_files[
  !file.exists(
    required_files
  )
]

if (length(missing_files)) {
  stop(
    "缺少必要RDS文件：\n",
    paste(
      missing_files,
      collapse = "\n"
    )
  )
}

if (
  !is.finite(
    KEEP_FRACTION
  ) ||
    KEEP_FRACTION <= 0 ||
    KEEP_FRACTION >= 1
) {
  stop(
    "KEEP_FRACTION必须在0和1之间。"
  )
}

if (
  !is.numeric(
    N_RANDOM
  ) ||
    N_RANDOM < 1
) {
  stop(
    "N_RANDOM必须>=1。"
  )
}

log_msg(
  "第十三阶段开始"
)

# ----------------------------- 基础函数 ----------------------------------------
collapse_duplicate_symbols <- function(
  counts,
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
      )
  )

  counts <- counts[
    valid,
    ,
    drop = FALSE
  ]

  gene_symbol <- gene_symbol[
    valid
  ]

  if (!is.matrix(counts)) {
    counts <- as.matrix(
      counts
    )
  }

  collapsed <- rowsum(
    counts,
    group = gene_symbol,
    reorder = FALSE
  )

  storage.mode(
    collapsed
  ) <- "double"

  collapsed
}

log_cpm <- function(
  counts
) {
  counts <- as.matrix(
    counts
  )

  storage.mode(
    counts
  ) <- "double"

  lib <- colSums(
    counts
  )

  if (any(
    !is.finite(
      lib
    ) |
      lib <= 0
  )) {
    stop(
      "发现文库总计数<=0或非有限值。"
    )
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

exclude_background_gene <- function(
  gene
) {
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

get_reference_design <- function(
  counts,
  sample_info
) {
  sample_info <- as.data.frame(
    sample_info,
    stringsAsFactors = FALSE
  )

  required_cols <- c(
    "sample_id",
    "main_cell_type",
    "condition"
  )

  missing_cols <- setdiff(
    required_cols,
    names(
      sample_info
    )
  )

  if (length(
    missing_cols
  )) {
    stop(
      "参考sample_info缺少字段：",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }

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
      "参考pseudobulk列名无法全部匹配sample_info。"
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
      "参考pseudobulk设计顺序匹配失败。"
    )
  }

  out
}

get_stable_fibroblast_signatures <- function(
  model
) {
  st <- model$signature_stability

  required_cols <- c(
    "cell_type",
    "stage",
    "gene",
    "stable_all_3_folds"
  )

  missing_cols <- setdiff(
    required_cols,
    names(
      st
    )
  )

  if (length(
    missing_cols
  )) {
    stop(
      "固定模型signature_stability缺少字段：",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  }

  x <- st[
    st$cell_type == "Fibroblast" &
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

  signatures <- setNames(
    lapply(
      stages,
      function(stage) {
        unique(
          as.character(
            x$gene[
              x$stage == stage
            ]
          )
        )
      }
    ),
    stages
  )

  n_gene <- vapply(
    signatures,
    length,
    integer(1)
  )

  if (any(
    n_gene < 20L
  )) {
    stop(
      "至少一个Fibroblast稳定签名少于20个基因：",
      paste(
        names(
          n_gene
        ),
        n_gene,
        sep = "=",
        collapse = ", "
      )
    )
  }

  signatures
}

make_rank_matrix <- function(
  expr
) {
  expr <- as.matrix(
    expr
  )

  out <- matrix(
    NA_real_,
    nrow = nrow(
      expr
    ),
    ncol = ncol(
      expr
    ),
    dimnames = dimnames(
      expr
    )
  )

  for (
    j in seq_len(
      ncol(
        expr
      )
    )
  ) {
    r <- rank(
      expr[
        ,
        j
      ],
      ties.method = "average",
      na.last = "keep"
    )

    denom <- sum(
      is.finite(
        r
      )
    )

    if (
      !is.finite(
        denom
      ) ||
        denom <= 0
    ) {
      stop(
        "秩评分背景基因为空。"
      )
    }

    out[
      ,
      j
    ] <- r / denom
  }

  out
}

score_from_rank_matrix <- function(
  rank_matrix,
  signatures
) {
  stages <- c(
    "Skin",
    "Wound1",
    "Wound7",
    "Wound30"
  )

  if (!all(
    stages %in% names(
      signatures
    )
  )) {
    stop(
      "签名缺少Skin/Wound1/Wound7/Wound30。"
    )
  }

  out <- matrix(
    NA_real_,
    nrow = ncol(
      rank_matrix
    ),
    ncol = length(
      stages
    ),
    dimnames = list(
      colnames(
        rank_matrix
      ),
      stages
    )
  )

  for (
    stage in stages
  ) {
    genes <- intersect(
      signatures[[stage]],
      rownames(
        rank_matrix
      )
    )

    if (length(
      genes
    ) < 10L) {
      stop(
        "阶段",
        stage,
        "可测签名基因不足10个，实际：",
        length(
          genes
        )
      )
    }

    out[
      ,
      stage
    ] <- colMeans(
      rank_matrix[
        genes,
        ,
        drop = FALSE
      ],
      na.rm = TRUE
    )
  }

  out
}

standardize_by_reference <- function(
  reference_scores,
  new_scores
) {
  mu <- colMeans(
    reference_scores,
    na.rm = TRUE
  )

  sigma <- apply(
    reference_scores,
    2,
    stats::sd,
    na.rm = TRUE
  )

  sigma[
    !is.finite(
      sigma
    ) |
      sigma == 0
  ] <- 1

  reference_z <- sweep(
    sweep(
      reference_scores,
      2,
      mu,
      "-"
    ),
    2,
    sigma,
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
    sigma,
    "/"
  )

  list(
    reference_z = reference_z,
    new_z = new_z,
    center = mu,
    scale = sigma
  )
}

hedges_g <- function(
  x,
  y
) {
  x <- x[
    is.finite(
      x
    )
  ]

  y <- y[
    is.finite(
      y
    )
  ]

  n1 <- length(
    x
  )

  n2 <- length(
    y
  )

  if (
    n1 < 2L ||
      n2 < 2L
  ) {
    return(
      NA_real_
    )
  }

  df <- n1 + n2 - 2L

  pooled_sd <- sqrt(
    (
      (n1 - 1) *
        stats::var(
          x
        ) +
        (n2 - 1) *
          stats::var(
            y
          )
    ) / df
  )

  if (
    !is.finite(
      pooled_sd
    ) ||
      pooled_sd <= 0
  ) {
    return(
      NA_real_
    )
  }

  d <- (
    mean(
      x
    ) -
      mean(
        y
      )
  ) / pooled_sd

  correction <- 1 -
    3 /
      (
        4 * df -
          1
      )

  correction * d
}

cliffs_delta <- function(
  x,
  y
) {
  x <- x[
    is.finite(
      x
    )
  ]

  y <- y[
    is.finite(
      y
    )
  ]

  if (
    !length(
      x
    ) ||
      !length(
        y
      )
  ) {
    return(
      NA_real_
    )
  }

  diff_mat <- outer(
    x,
    y,
    "-"
  )

  mean(
    diff_mat > 0
  ) -
    mean(
      diff_mat < 0
    )
}

exact_single_cohort <- function(
  value,
  group,
  expected_direction = 1
) {
  ok <- (
    is.finite(
      value
    ) &
      group %in% c(
        "keloid",
        "normal_scar"
      )
  )

  value <- value[
    ok
  ]

  group <- as.character(
    group[
      ok
    ]
  )

  x <- value[
    group == "keloid"
  ]

  y <- value[
    group == "normal_scar"
  ]

  if (
    length(
      x
    ) != 3L ||
      length(
        y
      ) != 3L
  ) {
    stop(
      "核心队列必须是3例keloid与3例normal scar；当前为 ",
      length(
        x
      ),
      " vs ",
      length(
        y
      )
    )
  }

  observed_difference <- mean(
    x
  ) -
    mean(
      y
    )

  values <- c(
    x,
    y
  )

  combinations <- utils::combn(
    seq_along(
      values
    ),
    3L
  )

  null_difference <- apply(
    combinations,
    2,
    function(
      keloid_indices
    ) {
      mean(
        values[
          keloid_indices
        ]
      ) -
        mean(
          values[
            -keloid_indices
          ]
        )
    }
  )

  oriented_observed <-
    expected_direction *
    observed_difference

  oriented_null <-
    expected_direction *
    null_difference

  complete_separation_expected <- if (
    expected_direction > 0
  ) {
    min(
      x
    ) >
      max(
        y
      )
  } else {
    max(
      x
    ) <
      min(
        y
      )
  }

  list(
    n_keloid =
      length(
        x
      ),
    n_normal_scar =
      length(
        y
      ),
    mean_keloid =
      mean(
        x
      ),
    mean_normal_scar =
      mean(
        y
      ),
    mean_difference =
      observed_difference,
    hedges_g =
      hedges_g(
        x,
        y
      ),
    cliffs_delta =
      cliffs_delta(
        x,
        y
      ),
    exact_two_sided_p =
      mean(
        abs(
          null_difference
        ) >=
          abs(
            observed_difference
          ) -
          1e-12
      ),
    exact_directional_p =
      mean(
        oriented_null >=
          oriented_observed -
          1e-12
      ),
    complete_separation_expected =
      complete_separation_expected,
    null_distribution =
      null_difference
  )
}

integrate_two_cohorts <- function(
  discovery_result,
  validation_result,
  expected_direction = 1
) {
  observed <- mean(
    c(
      discovery_result$mean_difference,
      validation_result$mean_difference
    )
  )

  null_grid <- expand.grid(
    discovery =
      discovery_result$null_distribution,
    validation =
      validation_result$null_distribution,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  null_statistic <- rowMeans(
    as.matrix(
      null_grid
    )
  )

  oriented_observed <-
    expected_direction *
    observed

  oriented_null <-
    expected_direction *
    null_statistic

  list(
    equal_weight_mean_difference =
      observed,
    mean_hedges_g =
      mean(
        c(
          discovery_result$hedges_g,
          validation_result$hedges_g
        ),
        na.rm = TRUE
      ),
    mean_cliffs_delta =
      mean(
        c(
          discovery_result$cliffs_delta,
          validation_result$cliffs_delta
        ),
        na.rm = TRUE
      ),
    exact_stratified_two_sided_p =
      mean(
        abs(
          null_statistic
        ) >=
          abs(
            observed
          ) -
          1e-12
      ),
    exact_stratified_directional_p =
      mean(
        oriented_null >=
          oriented_observed -
          1e-12
      ),
    n_exact_joint_permutations =
      length(
        null_statistic
      ),
    all_cohorts_expected_direction =
      (
        expected_direction *
          discovery_result$mean_difference >
          0
      ) &&
        (
          expected_direction *
            validation_result$mean_difference >
            0
        ),
    all_cohorts_complete_separation =
      isTRUE(
        discovery_result$complete_separation_expected
      ) &&
        isTRUE(
          validation_result$complete_separation_expected
        )
  )
}

prepare_rank_data <- function(
  reference_counts,
  reference_design,
  disease_counts,
  disease_group,
  disease_id,
  cohort_name
) {
  disease_counts <- as.matrix(
    disease_counts
  )

  if (
    is.null(
      rownames(
        disease_counts
      )
    ) ||
      is.null(
        colnames(
          disease_counts
        )
      )
  ) {
    stop(
      cohort_name,
      " disease_counts缺少行名或列名。"
    )
  }

  idx <- match(
    colnames(
      disease_counts
    ),
    disease_id
  )

  if (anyNA(
    idx
  )) {
    stop(
      cohort_name,
      " 疾病矩阵列名无法全部匹配疾病sample/patient ID。"
    )
  }

  disease_group <- as.character(
    disease_group[
      idx
    ]
  )

  disease_id <- as.character(
    disease_id[
      idx
    ]
  )

  keep_sample <- disease_group %in%
    c(
      "keloid",
      "normal_scar"
    )

  disease_counts <- disease_counts[
    ,
    keep_sample,
    drop = FALSE
  ]

  disease_group <- disease_group[
    keep_sample
  ]

  disease_id <- disease_id[
    keep_sample
  ]

  if (
    sum(
      disease_group == "keloid"
    ) != 3L ||
      sum(
        disease_group == "normal_scar"
      ) != 3L
  ) {
    stop(
      cohort_name,
      " 核心样本数不是3 vs 3。"
    )
  }

  common_genes <- intersect(
    rownames(
      reference_counts
    ),
    rownames(
      disease_counts
    )
  )

  common_genes <- common_genes[
    !exclude_background_gene(
      common_genes
    )
  ]

  if (length(
    common_genes
  ) < 5000L) {
    stop(
      cohort_name,
      " 参考与疾病共同背景基因少于5000，实际：",
      length(
        common_genes
      )
    )
  }

  reference_expr <- log_cpm(
    reference_counts[
      common_genes,
      ,
      drop = FALSE
    ]
  )

  disease_expr <- log_cpm(
    disease_counts[
      common_genes,
      ,
      drop = FALSE
    ]
  )

  list(
    cohort =
      cohort_name,
    common_genes =
      common_genes,
    reference_rank =
      make_rank_matrix(
        reference_expr
      ),
    disease_rank =
      make_rank_matrix(
        disease_expr
      ),
    disease_group =
      disease_group,
    disease_id =
      disease_id,
    reference_design =
      reference_design
  )
}

score_late_remodeling <- function(
  prepared,
  signatures
) {
  reference_scores <- score_from_rank_matrix(
    prepared$reference_rank,
    signatures
  )

  disease_scores <- score_from_rank_matrix(
    prepared$disease_rank,
    signatures
  )

  standardized <- standardize_by_reference(
    reference_scores,
    disease_scores
  )

  z <- standardized$new_z

  late <- z[
    ,
    "Wound30"
  ] -
    (
      z[
        ,
        "Wound1"
      ] +
        z[
          ,
          "Wound7"
        ]
    ) / 2

  data.frame(
    cohort =
      prepared$cohort,
    patient_id =
      prepared$disease_id,
    group =
      prepared$disease_group,
    z_Skin =
      z[
        ,
        "Skin"
      ],
    z_Wound1 =
      z[
        ,
        "Wound1"
      ],
    z_Wound7 =
      z[
        ,
        "Wound7"
      ],
    z_Wound30 =
      z[
        ,
        "Wound30"
      ],
    late_remodeling_state =
      late,
    stringsAsFactors = FALSE
  )
}

evaluate_signatures <- function(
  discovery_prepared,
  validation_prepared,
  signatures
) {
  discovery_scores <- score_late_remodeling(
    discovery_prepared,
    signatures
  )

  validation_scores <- score_late_remodeling(
    validation_prepared,
    signatures
  )

  discovery_test <- exact_single_cohort(
    discovery_scores$late_remodeling_state,
    discovery_scores$group,
    expected_direction = 1
  )

  validation_test <- exact_single_cohort(
    validation_scores$late_remodeling_state,
    validation_scores$group,
    expected_direction = 1
  )

  integrated <- integrate_two_cohorts(
    discovery_test,
    validation_test,
    expected_direction = 1
  )

  list(
    discovery_scores =
      discovery_scores,
    validation_scores =
      validation_scores,
    discovery_test =
      discovery_test,
    validation_test =
      validation_test,
    integrated =
      integrated
  )
}

flatten_result <- function(
  id,
  result,
  baseline_integrated_difference,
  omitted_gene = NA_character_,
  n_w30_genes = NA_integer_
) {
  data.frame(
    iteration_id =
      id,
    omitted_gene =
      omitted_gene,
    n_Wound30_genes_used =
      n_w30_genes,

    discovery_mean_difference =
      result$discovery_test$mean_difference,
    discovery_hedges_g =
      result$discovery_test$hedges_g,
    discovery_cliffs_delta =
      result$discovery_test$cliffs_delta,
    discovery_exact_two_sided_p =
      result$discovery_test$exact_two_sided_p,
    discovery_complete_separation =
      result$discovery_test$complete_separation_expected,

    validation_mean_difference =
      result$validation_test$mean_difference,
    validation_hedges_g =
      result$validation_test$hedges_g,
    validation_cliffs_delta =
      result$validation_test$cliffs_delta,
    validation_exact_two_sided_p =
      result$validation_test$exact_two_sided_p,
    validation_complete_separation =
      result$validation_test$complete_separation_expected,

    integrated_mean_difference =
      result$integrated$equal_weight_mean_difference,
    integrated_mean_hedges_g =
      result$integrated$mean_hedges_g,
    integrated_mean_cliffs_delta =
      result$integrated$mean_cliffs_delta,
    integrated_exact_two_sided_p =
      result$integrated$exact_stratified_two_sided_p,
    integrated_exact_directional_p =
      result$integrated$exact_stratified_directional_p,
    both_cohorts_expected_direction =
      result$integrated$all_cohorts_expected_direction,
    both_cohorts_complete_separation =
      result$integrated$all_cohorts_complete_separation,

    integrated_effect_fraction_of_baseline =
      if (
        is.finite(
          baseline_integrated_difference
        ) &&
          baseline_integrated_difference != 0
      ) {
        result$integrated$equal_weight_mean_difference /
          baseline_integrated_difference
      } else {
        NA_real_
      },

    stringsAsFactors = FALSE
  )
}

# ----------------------------- 读取对象 ----------------------------------------
reference_pb <- readRDS(
  REFERENCE_PB_RDS
)

reference_model <- readRDS(
  REFERENCE_MODEL_RDS
)

discovery_object <- readRDS(
  DISCOVERY_RDS
)

validation_object <- readRDS(
  VALIDATION_RDS
)

# ----------------------------- 正常伤口参考：Fibroblast -------------------------
if (
  is.null(
    reference_pb$counts
  ) ||
    is.null(
      reference_pb$sample_info
    ) ||
    is.null(
      reference_pb$gene_symbol
    )
) {
  stop(
    "Stage3参考RDS结构不完整。"
  )
}

reference_design_all <- get_reference_design(
  reference_pb$counts,
  reference_pb$sample_info
)

reference_counts_all <- collapse_duplicate_symbols(
  reference_pb$counts,
  reference_pb$gene_symbol
)

fib_ref_idx <- which(
  reference_design_all$main_cell_type ==
    "Fibroblast"
)

if (!length(
  fib_ref_idx
)) {
  stop(
    "正常伤口参考中没有Fibroblast。"
  )
}

reference_counts_fib <- reference_counts_all[
  ,
  fib_ref_idx,
  drop = FALSE
]

reference_design_fib <- reference_design_all[
  fib_ref_idx,
  ,
  drop = FALSE
]

# ----------------------------- 固定签名 ----------------------------------------
signatures_full <- get_stable_fibroblast_signatures(
  reference_model
)

signature_summary <- data.frame(
  stage =
    names(
      signatures_full
    ),
  n_stable_genes =
    vapply(
      signatures_full,
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

log_msg(
  "固定Fibroblast Wound30基因数：",
  length(
    signatures_full$Wound30
  )
)

# ----------------------------- GSE163973发现队列 -------------------------------
if (
  is.null(
    discovery_object$disease_pseudobulk$counts
  ) ||
    is.null(
      discovery_object$disease_pseudobulk$sample_info
    )
) {
  stop(
    "Stage5 RDS缺少disease_pseudobulk。"
  )
}

discovery_counts_all <-
  discovery_object$disease_pseudobulk$counts

discovery_info_all <-
  as.data.frame(
    discovery_object$disease_pseudobulk$sample_info,
    stringsAsFactors = FALSE
  )

required_discovery_cols <- c(
  "sample_key",
  "main_cell_type",
  "group"
)

missing_discovery_cols <- setdiff(
  required_discovery_cols,
  names(
    discovery_info_all
  )
)

if (length(
  missing_discovery_cols
)) {
  stop(
    "Stage5 sample_info缺少字段：",
    paste(
      missing_discovery_cols,
      collapse = ", "
    )
  )
}

discovery_info_fib <- discovery_info_all[
  discovery_info_all$main_cell_type ==
    "Fibroblast",
  ,
  drop = FALSE
]

discovery_col_idx <- match(
  discovery_info_fib$sample_key,
  colnames(
    discovery_counts_all
  )
)

if (anyNA(
  discovery_col_idx
)) {
  stop(
    "Stage5 Fibroblast sample_key无法匹配counts列。"
  )
}

discovery_counts_fib <- discovery_counts_all[
  ,
  discovery_col_idx,
  drop = FALSE
]

# ----------------------------- GSE181316复制队列 -------------------------------
if (
  is.null(
    validation_object$high_specificity$patient_counts
  ) ||
    is.null(
      validation_object$high_specificity$patient_info
    )
) {
  stop(
    "Stage6B RDS缺少high_specificity patient_counts/patient_info。"
  )
}

validation_counts_fib <-
  validation_object$high_specificity$patient_counts

validation_info_fib <-
  as.data.frame(
    validation_object$high_specificity$patient_info,
    stringsAsFactors = FALSE
  )

if (!"group" %in%
  names(
    validation_info_fib
  )) {
  stop(
    "Stage6B patient_info缺少group字段。"
  )
}

id_candidates <- c(
  "patient_id",
  "sample_key",
  "sample_id"
)

id_col <- id_candidates[
  id_candidates %in%
    names(
      validation_info_fib
    )
][1]

if (
  length(
    id_col
  ) == 0 ||
    is.na(
      id_col
    )
) {
  stop(
    "Stage6B patient_info中找不到patient_id/sample_key/sample_id。"
  )
}

validation_patient_id <- as.character(
  validation_info_fib[[id_col]]
)

# ----------------------------- 准备秩矩阵 --------------------------------------
log_msg(
  "预计算两个核心队列的共同基因背景与样本内秩矩阵"
)

discovery_prepared <- prepare_rank_data(
  reference_counts =
    reference_counts_fib,
  reference_design =
    reference_design_fib,
  disease_counts =
    discovery_counts_fib,
  disease_group =
    discovery_info_fib$group,
  disease_id =
    discovery_info_fib$sample_key,
  cohort_name =
    "GSE163973"
)

validation_prepared <- prepare_rank_data(
  reference_counts =
    reference_counts_fib,
  reference_design =
    reference_design_fib,
  disease_counts =
    validation_counts_fib,
  disease_group =
    validation_info_fib$group,
  disease_id =
    validation_patient_id,
  cohort_name =
    "GSE181316"
)

coverage_table <- data.frame(
  cohort = c(
    "GSE163973",
    "GSE181316"
  ),
  n_common_background_genes = c(
    length(
      discovery_prepared$common_genes
    ),
    length(
      validation_prepared$common_genes
    )
  ),
  n_Wound30_stable_genes_total =
    length(
      signatures_full$Wound30
    ),
  n_Wound30_genes_measurable = c(
    sum(
      signatures_full$Wound30 %in%
        discovery_prepared$common_genes
    ),
    sum(
      signatures_full$Wound30 %in%
        validation_prepared$common_genes
    )
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  coverage_table,
  file.path(
    REPORT_DIR,
    "02_core_cohort_Wound30_gene_coverage.csv"
  )
)

# ----------------------------- 先复现原始主结果 ---------------------------------
log_msg(
  "复现完整固定签名主结果"
)

baseline <- evaluate_signatures(
  discovery_prepared,
  validation_prepared,
  signatures_full
)

baseline_table <- flatten_result(
  id =
    "FULL_SIGNATURE",
  result =
    baseline,
  baseline_integrated_difference =
    baseline$integrated$equal_weight_mean_difference,
  omitted_gene =
    NA_character_,
  n_w30_genes =
    length(
      signatures_full$Wound30
    )
)

safe_write_csv(
  baseline_table,
  file.path(
    REPORT_DIR,
    "03_baseline_full_signature_reproduction.csv"
  )
)

baseline_patient_scores <- rbind(
  baseline$discovery_scores,
  baseline$validation_scores
)

safe_write_csv(
  baseline_patient_scores,
  file.path(
    REPORT_DIR,
    "04_baseline_patient_level_scores.csv"
  )
)

# 与现有Stage11结果做数值审计；允许微小浮点误差。
KNOWN_DISCOVERY_DIFF <- 1.154000
KNOWN_VALIDATION_DIFF <- 0.582882
KNOWN_INTEGRATED_DIFF <- 0.868441

reproduction_audit <- data.frame(
  quantity = c(
    "GSE163973_mean_difference",
    "GSE181316_mean_difference",
    "integrated_equal_weight_mean_difference"
  ),
  expected_from_existing_pipeline = c(
    KNOWN_DISCOVERY_DIFF,
    KNOWN_VALIDATION_DIFF,
    KNOWN_INTEGRATED_DIFF
  ),
  stage13_recomputed = c(
    baseline$discovery_test$mean_difference,
    baseline$validation_test$mean_difference,
    baseline$integrated$equal_weight_mean_difference
  ),
  absolute_difference = abs(
    c(
      baseline$discovery_test$mean_difference,
      baseline$validation_test$mean_difference,
      baseline$integrated$equal_weight_mean_difference
    ) -
      c(
        KNOWN_DISCOVERY_DIFF,
        KNOWN_VALIDATION_DIFF,
        KNOWN_INTEGRATED_DIFF
      )
  ),
  stringsAsFactors = FALSE
)

reproduction_audit$within_tolerance <-
  reproduction_audit$absolute_difference <
  0.005

safe_write_csv(
  reproduction_audit,
  file.path(
    REPORT_DIR,
    "05_baseline_reproduction_audit.csv"
  )
)

if (!all(
  reproduction_audit$within_tolerance
)) {
  warn_msg(
    "Stage13完整签名复现结果与现有Stage5/6B/11结果偏差>0.005；",
    "请先检查RDS版本、high_specificity集合或上游对象是否被替换。"
  )
}

# ----------------------------- LOGO：逐个删除Wound30基因 ------------------------
w30_genes <- signatures_full$Wound30

if (length(
  w30_genes
) < 20L) {
  stop(
    "Wound30稳定签名基因过少，停止LOGO。"
  )
}

logo_rows <- vector(
  "list",
  length(
    w30_genes
  )
)

log_msg(
  "开始leave-one-Wound30-gene-out，共",
  length(
    w30_genes
  ),
  "次"
)

for (
  i in seq_along(
    w30_genes
  )
) {
  gene <- w30_genes[
    i
  ]

  signatures_i <- signatures_full

  signatures_i$Wound30 <- setdiff(
    signatures_i$Wound30,
    gene
  )

  result_i <- evaluate_signatures(
    discovery_prepared,
    validation_prepared,
    signatures_i
  )

  row_i <- flatten_result(
    id =
      paste0(
        "LOGO_",
        i
      ),
    result =
      result_i,
    baseline_integrated_difference =
      baseline$integrated$equal_weight_mean_difference,
    omitted_gene =
      gene,
    n_w30_genes =
      length(
        signatures_i$Wound30
      )
  )

  row_i$gene_measurable_in_discovery <-
    gene %in%
    discovery_prepared$common_genes

  row_i$gene_measurable_in_validation <-
    gene %in%
    validation_prepared$common_genes

  logo_rows[[i]] <- row_i
}

logo_results <- do.call(
  rbind,
  logo_rows
)

logo_results <- logo_results[
  order(
    logo_results$integrated_effect_fraction_of_baseline
  ),
  ,
  drop = FALSE
]

safe_write_csv(
  logo_results,
  file.path(
    REPORT_DIR,
    "06_leave_one_Wound30_gene_out_results.csv"
  )
)

# ----------------------------- 随机80% Wound30子采样 ----------------------------
set.seed(
  RANDOM_SEED
)

n_keep <- max(
  10L,
  floor(
    length(
      w30_genes
    ) *
      KEEP_FRACTION
  )
)

if (
  n_keep >=
    length(
      w30_genes
    )
) {
  stop(
    "随机子采样保留基因数必须小于完整Wound30签名基因数。"
  )
}

random_rows <- vector(
  "list",
  N_RANDOM
)

random_gene_sets <- vector(
  "list",
  N_RANDOM
)

log_msg(
  "开始随机Wound30子采样：",
  N_RANDOM,
  "次；每次保留",
  n_keep,
  "/",
  length(
    w30_genes
  ),
  "个基因"
)

for (
  b in seq_len(
    N_RANDOM
  )
) {
  genes_b <- sample(
    w30_genes,
    size = n_keep,
    replace = FALSE
  )

  signatures_b <- signatures_full

  signatures_b$Wound30 <- genes_b

  result_b <- evaluate_signatures(
    discovery_prepared,
    validation_prepared,
    signatures_b
  )

  random_rows[[b]] <- flatten_result(
    id =
      paste0(
        "RANDOM_",
        b
      ),
    result =
      result_b,
    baseline_integrated_difference =
      baseline$integrated$equal_weight_mean_difference,
    omitted_gene =
      NA_character_,
    n_w30_genes =
      length(
        genes_b
      )
  )

  random_gene_sets[[b]] <- data.frame(
    iteration_id =
      paste0(
        "RANDOM_",
        b
      ),
    gene =
      genes_b,
    stringsAsFactors = FALSE
  )

  if (
    b %% 100L ==
      0L
  ) {
    log_msg(
      "随机子采样完成：",
      b,
      "/",
      N_RANDOM
    )
  }
}

random_results <- do.call(
  rbind,
  random_rows
)

random_gene_table <- do.call(
  rbind,
  random_gene_sets
)

safe_write_csv(
  random_results,
  file.path(
    REPORT_DIR,
    "07_random_80pct_Wound30_subsampling_results.csv"
  )
)

safe_write_csv(
  random_gene_table,
  file.path(
    REPORT_DIR,
    "08_random_80pct_Wound30_gene_sets.csv"
  )
)

# ----------------------------- 汇总稳健性 --------------------------------------
logo_summary <- data.frame(
  analysis =
    "leave_one_Wound30_gene_out",
  n_iterations =
    nrow(
      logo_results
    ),
  fraction_both_cohorts_expected_direction =
    mean(
      logo_results$both_cohorts_expected_direction
    ),
  fraction_both_cohorts_complete_separation =
    mean(
      logo_results$both_cohorts_complete_separation
    ),
  median_integrated_effect_fraction_of_baseline =
    stats::median(
      logo_results$integrated_effect_fraction_of_baseline,
      na.rm = TRUE
    ),
  minimum_integrated_effect_fraction_of_baseline =
    min(
      logo_results$integrated_effect_fraction_of_baseline,
      na.rm = TRUE
    ),
  maximum_integrated_effect_fraction_of_baseline =
    max(
      logo_results$integrated_effect_fraction_of_baseline,
      na.rm = TRUE
    ),
  fraction_integrated_effect_ge_80pct_baseline =
    mean(
      logo_results$integrated_effect_fraction_of_baseline >=
        0.80
    ),
  fraction_integrated_exact_p_le_0_01 =
    mean(
      logo_results$integrated_exact_two_sided_p <=
        0.01
    ),
  stringsAsFactors = FALSE
)

random_summary <- data.frame(
  analysis =
    paste0(
      "random_",
      round(
        KEEP_FRACTION *
          100
      ),
      "pct_Wound30_subsampling"
    ),
  n_iterations =
    nrow(
      random_results
    ),
  fraction_both_cohorts_expected_direction =
    mean(
      random_results$both_cohorts_expected_direction
    ),
  fraction_both_cohorts_complete_separation =
    mean(
      random_results$both_cohorts_complete_separation
    ),
  median_integrated_effect_fraction_of_baseline =
    stats::median(
      random_results$integrated_effect_fraction_of_baseline,
      na.rm = TRUE
    ),
  minimum_integrated_effect_fraction_of_baseline =
    min(
      random_results$integrated_effect_fraction_of_baseline,
      na.rm = TRUE
    ),
  maximum_integrated_effect_fraction_of_baseline =
    max(
      random_results$integrated_effect_fraction_of_baseline,
      na.rm = TRUE
    ),
  fraction_integrated_effect_ge_80pct_baseline =
    mean(
      random_results$integrated_effect_fraction_of_baseline >=
        0.80
    ),
  fraction_integrated_exact_p_le_0_01 =
    mean(
      random_results$integrated_exact_two_sided_p <=
        0.01
    ),
  stringsAsFactors = FALSE
)

robustness_summary <- rbind(
  logo_summary,
  random_summary
)

safe_write_csv(
  robustness_summary,
  file.path(
    REPORT_DIR,
    "09_signature_robustness_summary.csv"
  )
)

# ----------------------------- 单基因影响排序 -----------------------------------
gene_influence <- logo_results[
  ,
  c(
    "omitted_gene",
    "gene_measurable_in_discovery",
    "gene_measurable_in_validation",
    "discovery_mean_difference",
    "validation_mean_difference",
    "integrated_mean_difference",
    "integrated_effect_fraction_of_baseline",
    "both_cohorts_expected_direction",
    "both_cohorts_complete_separation",
    "integrated_exact_two_sided_p"
  )
]

gene_influence$absolute_loss_fraction <-
  1 -
  gene_influence$integrated_effect_fraction_of_baseline

gene_influence <- gene_influence[
  order(
    -gene_influence$absolute_loss_fraction
  ),
  ,
  drop = FALSE
]

safe_write_csv(
  gene_influence,
  file.path(
    REPORT_DIR,
    "10_Wound30_single_gene_influence_ranking.csv"
  )
)

# ----------------------------- 图1：LOGO整合效应 --------------------------------
pdf(
  file.path(
    FIG_DIR,
    "Figure_LOGO_integrated_effect_stability.pdf"
  ),
  width = 10,
  height = 6
)

op <- par(
  mar = c(
    9,
    5,
    3,
    1
  )
)

plot(
  seq_len(
    nrow(
      logo_results
    )
  ),
  logo_results$integrated_mean_difference,
  pch = 16,
  xaxt = "n",
  xlab = "",
  ylab = "Equal-weight integrated mean difference",
  main = "Leave-one-Wound30-gene-out robustness"
)

axis(
  1,
  at = seq_len(
    nrow(
      logo_results
    )
  ),
  labels =
    logo_results$omitted_gene,
  las = 2,
  cex.axis = 0.65
)

abline(
  h =
    baseline$integrated$equal_weight_mean_difference,
  lty = 2
)

abline(
  h = 0,
  lty = 3
)

par(
  op
)

dev.off()

# ----------------------------- 图2：随机子采样分布 -------------------------------
pdf(
  file.path(
    FIG_DIR,
    "Figure_random_subsampling_integrated_effect_distribution.pdf"
  ),
  width = 7,
  height = 6
)

hist(
  random_results$integrated_mean_difference,
  breaks = 30,
  xlab = "Equal-weight integrated mean difference",
  main = paste0(
    "Random ",
    round(
      KEEP_FRACTION *
        100
    ),
    "% Wound30-gene subsampling"
  )
)

abline(
  v =
    baseline$integrated$equal_weight_mean_difference,
  lty = 2
)

abline(
  v = 0,
  lty = 3
)

dev.off()

# ----------------------------- 图3：相对baseline效应 -----------------------------
pdf(
  file.path(
    FIG_DIR,
    "Figure_random_subsampling_relative_effect.pdf"
  ),
  width = 7,
  height = 6
)

boxplot(
  random_results$integrated_effect_fraction_of_baseline,
  ylab = "Integrated effect / full-signature effect",
  main = "Relative stability of late-remodeling effect"
)

abline(
  h = 1,
  lty = 2
)

abline(
  h = 0.8,
  lty = 3
)

dev.off()

# ----------------------------- 自动判定 ----------------------------------------
# 这是稳健性“描述性闸门”，不作为新的确认性P值门槛。
# 推荐判据：
#   A. LOGO所有/绝大多数迭代两个队列保持正方向；
#   B. LOGO最小整合效应仍>=完整签名的80%；
#   C. 随机80%子采样>=95%迭代两个队列同方向；
#   D. 随机80%子采样>=90%迭代保留至少80%的baseline整合效应。
#
# 注意：这些是本次敏感性分析的解释规则，不应倒装成原研究的“预注册阈值”。

logo_direction_pass <-
  logo_summary$fraction_both_cohorts_expected_direction >=
  0.95

logo_effect_pass <-
  logo_summary$minimum_integrated_effect_fraction_of_baseline >=
  0.80

random_direction_pass <-
  random_summary$fraction_both_cohorts_expected_direction >=
  0.95

random_effect_pass <-
  random_summary$fraction_integrated_effect_ge_80pct_baseline >=
  0.90

overall_robustness_pass <- all(
  logo_direction_pass,
  logo_effect_pass,
  random_direction_pass,
  random_effect_pass
)

decision_table <- data.frame(
  criterion = c(
    "LOGO >=95% iterations both cohorts expected direction",
    "LOGO minimum integrated effect >=80% baseline",
    "Random subsampling >=95% iterations both cohorts expected direction",
    "Random subsampling >=90% iterations integrated effect >=80% baseline"
  ),
  pass = c(
    logo_direction_pass,
    logo_effect_pass,
    random_direction_pass,
    random_effect_pass
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  decision_table,
  file.path(
    REPORT_DIR,
    "11_STAGE13_robustness_criteria.csv"
  )
)

most_influential <- head(
  gene_influence,
  10L
)

decision_lines <- c(
  "第十三阶段：Wound30固定签名稳健性敏感性分析结论",
  "============================================================",
  paste0(
    "运行时间：",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  "",
  "一、完整固定签名复现",
  paste0(
    "GSE163973 mean difference = ",
    sprintf(
      "%.6f",
      baseline$discovery_test$mean_difference
    )
  ),
  paste0(
    "GSE181316 mean difference = ",
    sprintf(
      "%.6f",
      baseline$validation_test$mean_difference
    )
  ),
  paste0(
    "Integrated equal-weight mean difference = ",
    sprintf(
      "%.6f",
      baseline$integrated$equal_weight_mean_difference
    )
  ),
  paste0(
    "Integrated exact two-sided P = ",
    format(
      baseline$integrated$exact_stratified_two_sided_p,
      scientific = FALSE
    )
  ),
  "",
  "二、Leave-one-Wound30-gene-out",
  paste0(
    "两个队列均保持预期正方向比例 = ",
    sprintf(
      "%.1f%%",
      100 *
        logo_summary$fraction_both_cohorts_expected_direction
    )
  ),
  paste0(
    "两个队列均保持完全分离比例 = ",
    sprintf(
      "%.1f%%",
      100 *
        logo_summary$fraction_both_cohorts_complete_separation
    )
  ),
  paste0(
    "最小整合效应 / baseline = ",
    sprintf(
      "%.3f",
      logo_summary$minimum_integrated_effect_fraction_of_baseline
    )
  ),
  "",
  "三、随机子采样",
  paste0(
    "随机次数 = ",
    N_RANDOM,
    "；每次保留 ",
    n_keep,
    "/",
    length(
      w30_genes
    ),
    " 个Wound30基因"
  ),
  paste0(
    "两个队列均保持预期正方向比例 = ",
    sprintf(
      "%.1f%%",
      100 *
        random_summary$fraction_both_cohorts_expected_direction
    )
  ),
  paste0(
    "两个队列均保持完全分离比例 = ",
    sprintf(
      "%.1f%%",
      100 *
        random_summary$fraction_both_cohorts_complete_separation
    )
  ),
  paste0(
    "整合效应>=baseline 80%的比例 = ",
    sprintf(
      "%.1f%%",
      100 *
        random_summary$fraction_integrated_effect_ge_80pct_baseline
    )
  ),
  "",
  "四、总体判定",
  if (
    overall_robustness_pass
  ) {
    paste(
      "ROBUST:",
      "晚期重塑主信号对单个Wound30基因删除及随机基因子采样均较稳定；",
      "没有证据提示主结论由某一个单独Wound30基因驱动。"
    )
  } else {
    paste(
      "SENSITIVITY_WARNING:",
      "至少一个预设稳健性描述规则未通过；",
      "应查看单基因影响排名和随机子采样分布后再决定如何写入正文。"
    )
  },
  "",
  "五、解释边界",
  paste(
    "本分析检验的是固定Wound30转录程序的内部稳健性，",
    "不能替代独立实验验证，也不能证明POSTN、ADAM12或其他单个基因具有因果作用。"
  )
)

safe_write_lines(
  decision_lines,
  file.path(
    REPORT_DIR,
    "12_STAGE13_DECISION.txt"
  )
)

# ----------------------------- 保存RDS -----------------------------------------
saveRDS(
  list(
    configuration = list(
      N_RANDOM =
        N_RANDOM,
      KEEP_FRACTION =
        KEEP_FRACTION,
      RANDOM_SEED =
        RANDOM_SEED
    ),
    signatures_full =
      signatures_full,
    coverage =
      coverage_table,
    baseline =
      baseline,
    baseline_table =
      baseline_table,
    reproduction_audit =
      reproduction_audit,
    logo_results =
      logo_results,
    random_results =
      random_results,
    robustness_summary =
      robustness_summary,
    gene_influence =
      gene_influence,
    decision_table =
      decision_table,
    overall_robustness_pass =
      overall_robustness_pass
  ),
  file = file.path(
    OBJECT_DIR,
    "Keloid_WoundAge_STAGE13_signature_robustness.rds"
  ),
  compress = "gzip"
)

# ----------------------------- Session info ------------------------------------
capture.output(
  sessionInfo(),
  file =
    file.path(
      REPORT_DIR,
      "13_SESSION_INFO.txt"
    )
)

safe_write_lines(
  c(
    "STAGE13_COMPLETED",
    paste0(
      "time=",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    paste0(
      "overall_robustness_pass=",
      overall_robustness_pass
    )
  ),
  file.path(
    REPORT_DIR,
    "STAGE13_COMPLETED.txt"
  )
)

# ----------------------------- 打包检查包 --------------------------------------
if (file.exists(
  PACKAGE_ZIP
)) {
  unlink(
    PACKAGE_ZIP,
    force = TRUE
  )
}

old_wd <- getwd()

on.exit(
  setwd(
    old_wd
  ),
  add = TRUE
)

setwd(
  STAGE_DIR
)

zip_files <- list.files(
  STAGE_DIR,
  recursive = TRUE,
  full.names = FALSE
)

zip_files <- zip_files[
  !grepl(
    "^objects/",
    gsub(
      "\\\\",
      "/",
      zip_files
    )
  )
]

tryCatch(
  {
    utils::zip(
      zipfile =
        PACKAGE_ZIP,
      files =
        zip_files
    )
  },
  error = function(e) {
    warn_msg(
      "ZIP打包失败：",
      conditionMessage(
        e
      )
    )
  }
)

setwd(
  old_wd
)

log_msg(
  "第十三阶段完成"
)

cat(
  paste(
    decision_lines,
    collapse = "\n"
  ),
  "\n"
)
