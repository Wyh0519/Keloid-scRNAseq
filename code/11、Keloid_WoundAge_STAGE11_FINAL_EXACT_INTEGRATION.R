# =============================================================================
# 项目：瘢痕疙瘩“晚期重塑状态持续与伤口状态终止失败”
# 第十一阶段：核心队列精确整合、证据分级与最终分析锁定
#
# 根目录：
#   C:/Users/33652/Desktop/LZZ/1、yssj
#
# 核心统计对象：
#   GSE163973：发现队列，3例keloid vs 3例normal scar
#   GSE181316：预设终点独立复制队列，3例keloid vs 3例normal scar
#
# 整合原则：
#   1. 只整合两个患者级核心队列；
#   2. 每个队列内部保留3 vs 3标签结构；
#   3. 使用分层精确置换：每队列C(6,3)=20种，共400种联合置换；
#   4. 统计量为两个队列组间均值差的等权平均；
#   5. 不把GSE220300、GSE181297、空间spots或单个细胞纳入主要整合；
#   6. 不把该分析称为传统随机效应Meta分析；
#   7. 主要终点和方向均已在独立复制前冻结。
#
# 证据角色：
#   Core:
#     GSE163973 + GSE181316
#   Supportive keloid:
#     GSE220300
#     GSE181297 Visium（表达层面支持）
#   Normal-reference orthogonal validation:
#     GSE241124 spatial transcriptomics
#   Disease control:
#     GSE265972 venous ulcer
#   Mixed/negative evidence:
#     GSE181297 scRNA high-specificity
#
# 本阶段输出：
#   - 两个核心队列患者级数据；
#   - 各队列效应量；
#   - 400种分层精确置换整合；
#   - 支持性证据矩阵；
#   - 可写与不可写结论边界；
#   - 最终分析完成闸门；
#   - 审计图和检查包。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANG = "en")

ROOT_DIR <- "C:/Users/33652/Desktop/LZZ/1、yssj"

STAGE_DIR <- file.path(
  ROOT_DIR,
  "11_STAGE11_FINAL_INTEGRATION"
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
  "第十一阶段_核心队列精确整合与最终证据锁定检查包.zip"
)

PACKAGE_TARGZ <- file.path(
  ROOT_DIR,
  "第十一阶段_核心队列精确整合与最终证据锁定检查包.tar.gz"
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
  "00_stage11_log.txt"
)

WARNING_FILE <- file.path(
  REPORT_DIR,
  "00_stage11_warnings.txt"
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

# ----------------------------- 文件定位 ----------------------------------------
find_report_file <- function(
  target_basename,
  preferred_fragment
) {
  all_files <- list.files(
    ROOT_DIR,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  )

  candidates <- all_files[
    tolower(
      basename(
        all_files
      )
    ) ==
      tolower(
        target_basename
      )
  ]

  if (!length(
    candidates
  )) {
    stop(
      "未找到必要报告文件：",
      target_basename
    )
  }

  preferred <- candidates[
    grepl(
      preferred_fragment,
      candidates,
      fixed = TRUE
    )
  ]

  if (length(
    preferred
  )) {
    candidates <- preferred
  }

  if (length(
    candidates
  ) >
      1L) {
    file_info <- file.info(
      candidates
    )

    selected <- candidates[
      order(
        file_info$mtime,
        decreasing = TRUE
      )[1]
    ]

    warn_msg(
      "发现多个同名文件：",
      target_basename,
      "；选择最新文件：",
      selected
    )
  } else {
    selected <- candidates[1]
  }

  selected
}

input_registry <- data.frame(
  stage = c(
    "Stage5",
    "Stage6B",
    "Stage7",
    "Stage8_overall",
    "Stage8_donor",
    "Stage9",
    "Stage10"
  ),
  target_basename = c(
    "09_GSE163973_patient_wound_state_scores.csv",
    "04_patient_level_scores_after_K3_merge.csv",
    "13_keloid_lesion_vs_scar_reference_descriptive.csv",
    "10_spatial_validation_overall_metrics.csv",
    "11_spatial_validation_donor_metrics.csv",
    "11_GSE265972_group_comparisons.csv",
    "11_cross_modal_direction_summary.csv"
  ),
  preferred_fragment = c(
    "05_STAGE5_GSE163973_MAPPING",
    "06B_STAGE6B_GSE181316_REPLICATION",
    "07_STAGE7_GSE220300_ACTIVITY_REGION",
    "08_STAGE8_GSE241124_SPATIAL_REFERENCE",
    "08_STAGE8_GSE241124_SPATIAL_REFERENCE",
    "09_STAGE9_GSE265972_REVERSE_CONTROL",
    "10_STAGE10_GSE181297_CROSS_MODAL_SUPPORT"
  ),
  stringsAsFactors = FALSE
)

input_registry$resolved_path <- vapply(
  seq_len(
    nrow(
      input_registry
    )
  ),
  function(i) {
    find_report_file(
      input_registry$target_basename[i],
      input_registry$preferred_fragment[i]
    )
  },
  character(1)
)

safe_write_csv(
  input_registry,
  file.path(
    REPORT_DIR,
    "01_input_report_registry.csv"
  )
)

get_input_path <- function(stage) {
  input_registry$resolved_path[
    input_registry$stage ==
      stage
  ][1]
}

# ----------------------------- 通用统计函数 ------------------------------------
normalize_group <- function(x) {
  y <- tolower(
    trimws(
      as.character(
        x
      )
    )
  )

  y[
    y %in%
      c(
        "kl",
        "keloid",
        "keloid_lesion"
      )
  ] <- "keloid"

  y[
    y %in%
      c(
        "ns",
        "normal_scar",
        "scar_reference"
      )
  ] <- "normal_scar"

  y
}

cliffs_delta <- function(
  x,
  y
) {
  pairwise_difference <- outer(
    x,
    y,
    "-"
  )

  mean(
    pairwise_difference >
      0
  ) -
    mean(
      pairwise_difference <
        0
    )
}

hedges_g <- function(
  x,
  y
) {
  n1 <- length(
    x
  )

  n2 <- length(
    y
  )

  df <- n1 +
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
      df
  )

  if (
    !is.finite(
      pooled_sd
    ) ||
      pooled_sd <=
      0
  ) {
    return(
      NA_real_
    )
  }

  cohen_d <- (
    mean(
      x
    ) -
      mean(
        y
      )
  ) /
    pooled_sd

  correction <- 1 -
    3 /
    (
      4 *
        df -
        1
    )

  correction *
    cohen_d
}

exact_single_cohort <- function(
  value,
  group,
  expected_direction = 1
) {
  x <- value[
    group ==
      "keloid"
  ]

  y <- value[
    group ==
      "normal_scar"
  ]

  if (
    length(
      x
    ) !=
      3L ||
      length(
        y
      ) !=
      3L
  ) {
    stop(
      "核心队列必须是3例keloid与3例normal scar。"
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
    function(keloid_indices) {
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
      if (
        expected_direction >
          0
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
      },
    null_distribution =
      null_difference
  )
}

stratified_exact_integration <- function(
  cohort_list,
  metric,
  expected_direction
) {
  cohort_results <- lapply(
    cohort_list,
    function(data) {
      exact_single_cohort(
        data[[metric]],
        data$group,
        expected_direction
      )
    }
  )

  observed_differences <- vapply(
    cohort_results,
    function(x) {
      x$mean_difference
    },
    numeric(1)
  )

  observed_statistic <- mean(
    observed_differences
  )

  null_grid <- expand.grid(
    lapply(
      cohort_results,
      function(x) {
        x$null_distribution
      }
    ),
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
    observed_statistic

  oriented_null <-
    expected_direction *
    null_statistic

  cohort_effects <- do.call(
    rbind,
    lapply(
      seq_along(
        cohort_results
      ),
      function(i) {
        x <- cohort_results[[i]]

        data.frame(
          cohort =
            names(
              cohort_list
            )[i],
          metric =
            metric,
          expected_direction =
            expected_direction,
          n_keloid =
            x$n_keloid,
          n_normal_scar =
            x$n_normal_scar,
          mean_keloid =
            x$mean_keloid,
          mean_normal_scar =
            x$mean_normal_scar,
          mean_difference =
            x$mean_difference,
          hedges_g =
            x$hedges_g,
          cliffs_delta =
            x$cliffs_delta,
          exact_two_sided_p =
            x$exact_two_sided_p,
          exact_directional_p =
            x$exact_directional_p,
          complete_separation_expected =
            x$complete_separation_expected,
          stringsAsFactors = FALSE
        )
      }
    )
  )

  integrated <- data.frame(
    metric =
      metric,
    expected_direction =
      expected_direction,
    n_cohorts =
      length(
        cohort_list
      ),
    total_keloid =
      sum(
        cohort_effects$n_keloid
      ),
    total_normal_scar =
      sum(
        cohort_effects$n_normal_scar
      ),
    equal_weight_mean_difference =
      observed_statistic,
    mean_hedges_g =
      mean(
        cohort_effects$hedges_g,
        na.rm = TRUE
      ),
    mean_cliffs_delta =
      mean(
        cohort_effects$cliffs_delta,
        na.rm = TRUE
      ),
    exact_stratified_two_sided_p =
      mean(
        abs(
          null_statistic
        ) >=
          abs(
            observed_statistic
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
      all(
        expected_direction *
          cohort_effects$mean_difference >
          0
      ),
    all_cohorts_complete_separation =
      all(
        cohort_effects$complete_separation_expected
      ),
    stringsAsFactors = FALSE
  )

  list(
    cohort_effects =
      cohort_effects,
    integrated =
      integrated,
    null_distribution =
      null_statistic
  )
}

# ----------------------------- 读取核心队列 ------------------------------------
log_msg(
  "第十一阶段开始"
)

stage5 <- utils::read.csv(
  get_input_path(
    "Stage5"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stage6 <- utils::read.csv(
  get_input_path(
    "Stage6B"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_stage5_columns <- c(
  "sample_id",
  "group",
  "cell_type",
  "ordinal_position",
  "remodeling_completion",
  "z_Wound30",
  "wound_activation",
  "z_Skin",
  "early_state_persistence",
  "off_trajectory_ratio"
)

missing_stage5 <- setdiff(
  required_stage5_columns,
  names(
    stage5
  )
)

if (length(
  missing_stage5
)) {
  stop(
    "Stage5患者表缺少字段：",
    paste(
      missing_stage5,
      collapse = ", "
    )
  )
}

stage5_core <- stage5[
  stage5$cell_type ==
    "Fibroblast",
  ,
  drop = FALSE
]

stage5_core$cohort <-
  "GSE163973"

stage5_core$patient_id <-
  stage5_core$sample_id

stage5_core$group <-
  normalize_group(
    stage5_core$group
  )

stage5_core$late_remodeling_state <-
  stage5_core$remodeling_completion

stage5_core$ordinal_wound_state_position <-
  stage5_core$ordinal_position

required_stage6_columns <- c(
  "patient_id",
  "group",
  "analysis_set",
  "ordinal_wound_state_position",
  "late_remodeling_state",
  "z_Wound30",
  "wound_activation",
  "z_Skin",
  "early_state_persistence",
  "off_trajectory_ratio"
)

missing_stage6 <- setdiff(
  required_stage6_columns,
  names(
    stage6
  )
)

if (length(
  missing_stage6
)) {
  stop(
    "Stage6B患者表缺少字段：",
    paste(
      missing_stage6,
      collapse = ", "
    )
  )
}

stage6_core <- stage6[
  stage6$analysis_set ==
    "high_specificity" &
    stage6$group %in%
    c(
      "keloid",
      "normal_scar"
    ),
  ,
  drop = FALSE
]

stage6_core$cohort <-
  "GSE181316"

stage6_core$group <-
  normalize_group(
    stage6_core$group
  )

core_metrics <- data.frame(
  metric = c(
    "late_remodeling_state",
    "ordinal_wound_state_position",
    "z_Wound30",
    "wound_activation",
    "z_Skin",
    "early_state_persistence",
    "off_trajectory_ratio"
  ),
  endpoint_role = c(
    "primary",
    "key_secondary",
    "supporting",
    "supporting",
    "supporting",
    "exploratory",
    "exploratory"
  ),
  expected_direction = c(
    1,
    1,
    1,
    1,
    -1,
    1,
    1
  ),
  stringsAsFactors = FALSE
)

required_core_columns <- c(
  "cohort",
  "patient_id",
  "group",
  core_metrics$metric
)

stage5_core_export <- stage5_core[
  ,
  required_core_columns,
  drop = FALSE
]

stage6_core_export <- stage6_core[
  ,
  required_core_columns,
  drop = FALSE
]

core_patient_data <- rbind(
  stage5_core_export,
  stage6_core_export
)

rownames(
  core_patient_data
) <- NULL

if (
  nrow(
    stage5_core_export
  ) !=
    6L ||
    nrow(
      stage6_core_export
    ) !=
    6L ||
    !all(
      table(
        stage5_core_export$group
      ) ==
        3L
    ) ||
    !all(
      table(
        stage6_core_export$group
      ) ==
        3L
    )
) {
  stop(
    "核心患者结构不是两个独立3 vs 3队列。"
  )
}

safe_write_csv(
  core_patient_data,
  file.path(
    REPORT_DIR,
    "02_core_patient_level_data.csv"
  )
)

# ----------------------------- 精确分层整合 ------------------------------------
cohort_list <- list(
  GSE163973 =
    stage5_core_export,
  GSE181316 =
    stage6_core_export
)

integration_results <- lapply(
  seq_len(
    nrow(
      core_metrics
    )
  ),
  function(i) {
    stratified_exact_integration(
      cohort_list =
        cohort_list,
      metric =
        core_metrics$metric[i],
      expected_direction =
        core_metrics$expected_direction[i]
    )
  }
)

cohort_effects <- do.call(
  rbind,
  lapply(
    integration_results,
    function(x) {
      x$cohort_effects
    }
  )
)

integrated_effects <- do.call(
  rbind,
  lapply(
    integration_results,
    function(x) {
      x$integrated
    }
  )
)

integrated_effects <- merge(
  integrated_effects,
  core_metrics[
    ,
    c(
      "metric",
      "endpoint_role"
    )
  ],
  by = "metric",
  all.x = TRUE,
  sort = FALSE
)

integrated_effects <- integrated_effects[
  match(
    core_metrics$metric,
    integrated_effects$metric
  ),
  ,
  drop = FALSE
]

rownames(
  cohort_effects
) <- NULL

rownames(
  integrated_effects
) <- NULL

safe_write_csv(
  cohort_effects,
  file.path(
    REPORT_DIR,
    "03_core_cohort_specific_effects.csv"
  )
)

safe_write_csv(
  integrated_effects,
  file.path(
    REPORT_DIR,
    "04_stratified_exact_core_integration.csv"
  )
)

# 保存每个主要指标的400个联合置换统计量，供审计。
null_distribution_rows <- do.call(
  rbind,
  lapply(
    seq_along(
      integration_results
    ),
    function(i) {
      data.frame(
        metric =
          core_metrics$metric[i],
        permutation_id =
          seq_along(
            integration_results[[i]]$null_distribution
          ),
        null_equal_weight_mean_difference =
          integration_results[[i]]$null_distribution,
        stringsAsFactors = FALSE
      )
    }
  )
)

safe_write_csv(
  null_distribution_rows,
  file.path(
    REPORT_DIR,
    "05_exact_joint_permutation_distributions.csv"
  )
)

# ----------------------------- 支持性证据 --------------------------------------
stage7 <- utils::read.csv(
  get_input_path(
    "Stage7"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stage8_overall <- utils::read.csv(
  get_input_path(
    "Stage8_overall"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stage8_donor <- utils::read.csv(
  get_input_path(
    "Stage8_donor"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stage9 <- utils::read.csv(
  get_input_path(
    "Stage9"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stage10 <- utils::read.csv(
  get_input_path(
    "Stage10"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

support_rows <- list()

# GSE220300支持队列
for (
  set_name in c(
    "high_specificity",
    "broad_sensitivity"
  )
) {
  x <- stage7[
    stage7$analysis_set ==
      set_name &
    stage7$metric %in%
      c(
        "late_remodeling_state",
        "ordinal_wound_state_position",
        "z_Wound30",
        "wound_activation",
        "z_Skin",
        "early_state_persistence",
        "off_trajectory_ratio"
      ),
    ,
    drop = FALSE
  ]

  if (nrow(
    x
  )) {
    support_rows[[length(
      support_rows
    ) + 1L]] <- data.frame(
      evidence_domain =
        "supportive_keloid_cohort",
      dataset =
        "GSE220300",
      modality =
        "scRNA",
      analysis_set =
        set_name,
      metric =
        x$metric,
      estimate =
        x$mean_difference,
      secondary_estimate =
        x$cliffs_delta,
      direction =
        sign(
          x$mean_difference
        ),
      evidence_role =
        "supportive_not_formal_replication",
      stringsAsFactors = FALSE
    )
  }
}

# GSE241124正常空间验证
for (
  i in seq_len(
    nrow(
      stage8_overall
    )
  )
) {
  support_rows[[length(
    support_rows
  ) + 1L]] <- data.frame(
    evidence_domain =
      "normal_reference_validation",
    dataset =
      "GSE241124",
    modality =
      "Visium",
    analysis_set =
      stage8_overall$analysis_set[i],
    metric =
      c(
        "exact_accuracy",
        "adjacent_accuracy",
        "pooled_spearman"
      ),
    estimate = c(
      stage8_overall$exact_accuracy[i],
      stage8_overall$adjacent_accuracy[i],
      stage8_overall$pooled_spearman[i]
    ),
    secondary_estimate =
      NA_real_,
    direction =
      1,
    evidence_role =
      "orthogonal_reference_validation",
    stringsAsFactors = FALSE
  )
}

# GSE265972疾病对照
x9 <- stage9[
  stage9$analysis_set ==
    "FB1_FB2_combined_primary" &
  stage9$metric %in%
    c(
      "late_remodeling_state",
      "early_state_persistence",
      "late_vs_early_balance",
      "ordinal_wound_state_position",
      "off_trajectory_ratio"
    ),
  ,
  drop = FALSE
]

support_rows[[length(
  support_rows
) + 1L]] <- data.frame(
  evidence_domain =
    "disease_control",
  dataset =
    "GSE265972",
  modality =
    "scRNA",
  analysis_set =
    "FB1_FB2_combined_primary",
  metric =
    x9$metric,
  estimate =
    x9$mean_difference,
  secondary_estimate =
    x9$cliffs_delta,
  direction =
    sign(
      x9$mean_difference
    ),
  evidence_role =
    "failure_mode_contrast",
  stringsAsFactors = FALSE
)

# GSE181297混合跨模态支持
for (
  i in seq_len(
    nrow(
      stage10
    )
  )
) {
  support_rows[[length(
    support_rows
  ) + 1L]] <- rbind(
    data.frame(
      evidence_domain =
        "cross_modal_support",
      dataset =
        "GSE181297",
      modality =
        stage10$modality[i],
      analysis_set =
        stage10$analysis_set[i],
      metric =
        "late_remodeling_state",
      estimate =
        stage10$late_remodeling_mean_difference[i],
      secondary_estimate =
        stage10$late_remodeling_cliffs_delta[i],
      direction =
        sign(
          stage10$late_remodeling_mean_difference[i]
        ),
      evidence_role =
        if (
          stage10$modality[i] ==
            "Visium"
        ) {
          "section_level_expression_support"
        } else {
          "mixed_small_scRNA_evidence"
        },
      stringsAsFactors = FALSE
    ),
    data.frame(
      evidence_domain =
        "cross_modal_support",
      dataset =
        "GSE181297",
      modality =
        stage10$modality[i],
      analysis_set =
        stage10$analysis_set[i],
      metric =
        "ordinal_wound_state_position",
      estimate =
        stage10$ordinal_mean_difference[i],
      secondary_estimate =
        stage10$ordinal_cliffs_delta[i],
      direction =
        sign(
          stage10$ordinal_mean_difference[i]
        ),
      evidence_role =
        if (
          stage10$modality[i] ==
            "Visium"
        ) {
          "section_level_expression_support"
        } else {
          "mixed_small_scRNA_evidence"
        },
      stringsAsFactors = FALSE
    )
  )
}

support_evidence <- do.call(
  rbind,
  support_rows
)

rownames(
  support_evidence
) <- NULL

safe_write_csv(
  support_evidence,
  file.path(
    REPORT_DIR,
    "06_supportive_evidence_registry.csv"
  )
)

# ----------------------------- 证据等级表 --------------------------------------
primary_late <- integrated_effects[
  integrated_effects$metric ==
    "late_remodeling_state",
  ,
  drop = FALSE
]

primary_ordinal <- integrated_effects[
  integrated_effects$metric ==
    "ordinal_wound_state_position",
  ,
  drop = FALSE
]

support_d30 <- integrated_effects[
  integrated_effects$metric ==
    "z_Wound30",
  ,
  drop = FALSE
]

support_activation <- integrated_effects[
  integrated_effects$metric ==
    "wound_activation",
  ,
  drop = FALSE
]

support_skin <- integrated_effects[
  integrated_effects$metric ==
    "z_Skin",
  ,
  drop = FALSE
]

early_exploratory <- integrated_effects[
  integrated_effects$metric ==
    "early_state_persistence",
  ,
  drop = FALSE
]

offtrajectory_exploratory <- integrated_effects[
  integrated_effects$metric ==
    "off_trajectory_ratio",
  ,
  drop = FALSE
]

stage7_late_high <- stage7[
  stage7$analysis_set ==
    "high_specificity" &
    stage7$metric ==
    "late_remodeling_state",
  ,
  drop = FALSE
]

stage7_late_broad <- stage7[
  stage7$analysis_set ==
    "broad_sensitivity" &
    stage7$metric ==
    "late_remodeling_state",
  ,
  drop = FALSE
]

stage8_high <- stage8_overall[
  stage8_overall$analysis_set ==
    "high_specificity_spots",
  ,
  drop = FALSE
]

stage8_broad <- stage8_overall[
  stage8_overall$analysis_set ==
    "broad_top_quartile_spots",
  ,
  drop = FALSE
]

stage9_early <- stage9[
  stage9$analysis_set ==
    "FB1_FB2_combined_primary" &
    stage9$metric ==
    "early_state_persistence",
  ,
  drop = FALSE
]

stage9_balance <- stage9[
  stage9$analysis_set ==
    "FB1_FB2_combined_primary" &
    stage9$metric ==
    "late_vs_early_balance",
  ,
  drop = FALSE
]

stage10_visium_high <- stage10[
  stage10$analysis_set ==
    "Visium_high_specificity_spots",
  ,
  drop = FALSE
]

stage10_scrna_high <- stage10[
  stage10$analysis_set ==
    "scRNA_high_specificity",
  ,
  drop = FALSE
]

evidence_hierarchy <- data.frame(
  evidence_item = c(
    "Core late-remodeling integration",
    "Core ordinal-position integration",
    "Core D30-like supporting signal",
    "Core wound-activation supporting signal",
    "Core loss of Skin-like state",
    "GSE220300 high-specificity support",
    "GSE220300 broad support",
    "GSE241124 high spatial ordinal validation",
    "GSE241124 broad spatial ordinal validation",
    "GSE265972 early-state disease-control signal",
    "GSE265972 late-vs-early balance",
    "GSE181297 Visium high support",
    "GSE181297 scRNA high result",
    "Early-state persistence across core cohorts",
    "Off-trajectory across core cohorts"
  ),
  evidence_level = c(
    "core_confirmatory",
    "core_confirmatory",
    "core_supporting",
    "core_supporting",
    "core_supporting",
    "supportive",
    "supportive",
    "orthogonal_validation",
    "orthogonal_validation",
    "disease_control",
    "disease_control",
    "supportive",
    "mixed_negative",
    "exploratory_inconsistent",
    "exploratory_inconsistent"
  ),
  estimate = c(
    primary_late$equal_weight_mean_difference,
    primary_ordinal$equal_weight_mean_difference,
    support_d30$equal_weight_mean_difference,
    support_activation$equal_weight_mean_difference,
    support_skin$equal_weight_mean_difference,
    stage7_late_high$mean_difference,
    stage7_late_broad$mean_difference,
    stage8_high$pooled_spearman,
    stage8_broad$pooled_spearman,
    stage9_early$mean_difference,
    stage9_balance$mean_difference,
    stage10_visium_high$late_remodeling_mean_difference,
    stage10_scrna_high$late_remodeling_mean_difference,
    early_exploratory$equal_weight_mean_difference,
    offtrajectory_exploratory$equal_weight_mean_difference
  ),
  inferential_value = c(
    primary_late$exact_stratified_two_sided_p,
    primary_ordinal$exact_stratified_two_sided_p,
    support_d30$exact_stratified_two_sided_p,
    support_activation$exact_stratified_two_sided_p,
    support_skin$exact_stratified_two_sided_p,
    NA,
    NA,
    stage8_high$exact_accuracy,
    stage8_broad$exact_accuracy,
    stage9_early$exact_two_sided_p,
    stage9_balance$exact_two_sided_p,
    NA,
    NA,
    early_exploratory$exact_stratified_two_sided_p,
    offtrajectory_exploratory$exact_stratified_two_sided_p
  ),
  conclusion = c(
    "confirmed",
    "confirmed",
    "confirmed",
    "confirmed",
    "confirmed",
    "supportive_positive",
    "supportive_positive",
    "moderate_to_strong_validation",
    "moderate_to_strong_validation",
    "distinct_failure_mode_support",
    "early_dominance_support",
    "section_level_positive",
    "not_replicated_in_small_scRNA_subset",
    "not_a_stable_core_feature",
    "not_a_stable_core_feature"
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  evidence_hierarchy,
  file.path(
    REPORT_DIR,
    "07_final_evidence_hierarchy.csv"
  )
)

# ----------------------------- 结论边界 ----------------------------------------
claim_boundary <- data.frame(
  claim = c(
    "Keloid fibroblasts show a persistent late-remodeling state relative to mature/normal scars.",
    "The late-remodeling signal was independently replicated in two patient-level scRNA cohorts.",
    "Keloid fibroblasts shift toward a D30-like state and away from a Skin-like state.",
    "The fixed normal-wound signature shows ordinal validity in independent spatial transcriptomic sections.",
    "GSE220300 provides additional lesion-level support.",
    "GSE181297 Visium provides section-level expression support.",
    "Venous ulcers show a different balance, with stronger early-state persistence.",
    "The score measures exact molecular wound age in days.",
    "All keloids universally exhibit the same state.",
    "Clinical activity is proportional to the late-remodeling score.",
    "Keloid centers or peripheries consistently have higher late-remodeling activity.",
    "GSE181297 proves precise spatial localization of the state.",
    "The analysis proves a causal molecular mechanism.",
    "The analysis predicts postoperative recurrence.",
    "Early-state persistence is a stable defining feature of keloids.",
    "Off-trajectory distance is a stable defining feature of keloids."
  ),
  status = c(
    rep(
      "allowed",
      7
    ),
    rep(
      "not_allowed",
      9
    )
  ),
  manuscript_use = c(
    "Title/Abstract/Results/Discussion",
    "Abstract/Results/Discussion",
    "Abstract/Results/Discussion",
    "Results/Discussion",
    "Results/Discussion",
    "Supplementary or supporting Results",
    "Results/Discussion",
    "Do not use",
    "Do not use",
    "Do not use",
    "Do not use",
    "Do not use",
    "Do not use",
    "Do not use",
    "Exploratory only",
    "Exploratory only"
  ),
  reason = c(
    "Confirmed by exact stratified integration of GSE163973 and GSE181316.",
    "Both core cohorts showed the pre-specified direction with complete patient-level separation.",
    "D30, wound activation and Skin-like supporting metrics were concordant.",
    "GSE241124 recovered the ordinal process with 75-81% exact and 87.5% adjacent accuracy.",
    "Both high-specificity and broad analyses were directionally positive.",
    "Both keloid sections exceeded both adjacent-normal sections, but coordinates were unavailable.",
    "GSE265972 showed a stronger early-state signal and negative late-vs-early balance.",
    "Ordinal position is relative, not calibrated chronological age.",
    "Available cohorts are small and anatomically heterogeneous.",
    "GSE220300 active versus inactive differences were unstable.",
    "GSE220300 center-periphery directions were inconsistent.",
    "No tissue_positions or scalefactors were available.",
    "All data are observational and computational.",
    "No usable recurrence-expression cohort was available.",
    "The two core cohorts had opposite early-state directions.",
    "The two core cohorts had opposite off-trajectory directions."
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  claim_boundary,
  file.path(
    REPORT_DIR,
    "08_manuscript_claim_boundary_matrix.csv"
  )
)

# ----------------------------- 最终闸门 ----------------------------------------
core_confirmed <- (
  primary_late$exact_stratified_two_sided_p <=
    0.01 &&
    primary_late$all_cohorts_expected_direction &&
    primary_late$all_cohorts_complete_separation &&
    primary_ordinal$exact_stratified_two_sided_p <=
    0.01 &&
    support_d30$exact_stratified_two_sided_p <=
    0.01 &&
    support_activation$exact_stratified_two_sided_p <=
    0.01 &&
    support_skin$exact_stratified_two_sided_p <=
    0.01
)

supportive_keloid_positive <- (
  nrow(
    stage7_late_high
  ) ==
    1L &&
    nrow(
      stage7_late_broad
    ) ==
    1L &&
    stage7_late_high$mean_difference >
    0 &&
    stage7_late_broad$mean_difference >
    0
)

normal_spatial_valid <- (
  nrow(
    stage8_high
  ) ==
    1L &&
    nrow(
      stage8_broad
    ) ==
    1L &&
    stage8_high$exact_accuracy >=
    0.75 &&
    stage8_high$adjacent_accuracy >=
    0.875 &&
    stage8_high$pooled_spearman >=
    0.65 &&
    stage8_broad$exact_accuracy >=
    0.70 &&
    stage8_broad$adjacent_accuracy >=
    0.875
)

disease_control_distinct <- (
  nrow(
    stage9_early
  ) ==
    1L &&
    nrow(
      stage9_balance
    ) ==
    1L &&
    stage9_early$mean_difference >
    0 &&
    stage9_early$exact_two_sided_p <
    0.05 &&
    stage9_balance$mean_difference <
    0
)

stage10_correctly_bounded <- (
  nrow(
    stage10_visium_high
  ) ==
    1L &&
    nrow(
      stage10_scrna_high
    ) ==
    1L &&
    stage10_visium_high$late_remodeling_mean_difference >
    0 &&
    stage10_scrna_high$late_remodeling_mean_difference <
    0
)

if (
  core_confirmed &&
    supportive_keloid_positive &&
    normal_spatial_valid &&
    disease_control_distinct &&
    stage10_correctly_bounded
) {
  project_gate <-
    "ANALYSIS_COMPLETE_CORE_EVIDENCE_STRONG"

  interpretation <- paste(
    "两个核心患者级单细胞队列通过分层精确置换整合，",
    "确认瘢痕疙瘩成纤维细胞具有持续的D30-like晚期重塑状态并丧失Skin-like稳态；",
    "第三瘢痕疙瘩队列、正常伤口空间转录组和静脉溃疡疾病对照提供了互补支持。",
    "GSE181297仅保留为Visium阳性、scRNA高特异性未复制的混合补充证据。"
  )
} else if (
  core_confirmed
) {
  project_gate <-
    "ANALYSIS_COMPLETE_CORE_EVIDENCE_CONFIRMED"

  interpretation <- paste(
    "两个核心患者级队列已确认主要结论，",
    "但至少一个支持性、空间或疾病对照模块未达到全部预设标准。"
  )
} else {
  project_gate <-
    "CORE_EVIDENCE_NOT_CONFIRMED"

  interpretation <- paste(
    "两个核心队列的分层精确整合未达到预设确认标准；",
    "不得进入主要论文结论锁定。"
  )
}

final_criteria <- data.frame(
  criterion = c(
    "core_confirmed",
    "supportive_keloid_positive",
    "normal_spatial_valid",
    "disease_control_distinct",
    "stage10_correctly_bounded"
  ),
  passed = c(
    core_confirmed,
    supportive_keloid_positive,
    normal_spatial_valid,
    disease_control_distinct,
    stage10_correctly_bounded
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  final_criteria,
  file.path(
    REPORT_DIR,
    "09_final_analysis_criteria.csv"
  )
)

# ----------------------------- 审计图 ------------------------------------------
# 1. 两个核心队列患者级主要终点
pdf(
  file.path(
    FIG_DIR,
    "Figure_core_patient_level_primary_endpoints.pdf"
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

for (
  cohort in c(
    "GSE163973",
    "GSE181316"
  )
) {
  data <- core_patient_data[
    core_patient_data$cohort ==
      cohort,
    ,
    drop = FALSE
  ]

  group_position <- ifelse(
    data$group ==
      "normal_scar",
    1,
    2
  )

  plot(
    jitter(
      group_position,
      amount = 0.05
    ),
    data$late_remodeling_state,
    pch = ifelse(
      data$group ==
        "keloid",
      16,
      1
    ),
    xlim = c(
      0.5,
      2.5
    ),
    xaxt = "n",
    xlab = "",
    ylab =
      "Late-remodeling state",
    main =
      cohort
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
      group_position,
      amount = 0.05
    ),
    data$late_remodeling_state,
    labels =
      data$patient_id,
    pos = 3,
    cex = 0.75
  )
}

for (
  cohort in c(
    "GSE163973",
    "GSE181316"
  )
) {
  data <- core_patient_data[
    core_patient_data$cohort ==
      cohort,
    ,
    drop = FALSE
  ]

  group_position <- ifelse(
    data$group ==
      "normal_scar",
    1,
    2
  )

  plot(
    jitter(
      group_position,
      amount = 0.05
    ),
    data$ordinal_wound_state_position,
    pch = ifelse(
      data$group ==
        "keloid",
      16,
      1
    ),
    xlim = c(
      0.5,
      2.5
    ),
    xaxt = "n",
    xlab = "",
    ylab =
      "Ordinal wound-state position",
    main =
      cohort
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
      group_position,
      amount = 0.05
    ),
    data$ordinal_wound_state_position,
    labels =
      data$patient_id,
    pos = 3,
    cex = 0.75
  )
}

par(
  op
)

dev.off()

# 2. 核心队列效应方向图
plot_metrics <- c(
  "late_remodeling_state",
  "ordinal_wound_state_position",
  "z_Wound30",
  "wound_activation",
  "z_Skin"
)

effect_plot_data <- cohort_effects[
  cohort_effects$metric %in%
    plot_metrics,
  ,
  drop = FALSE
]

effect_matrix <- matrix(
  NA_real_,
  nrow = length(
    plot_metrics
  ),
  ncol = 2,
  dimnames = list(
    plot_metrics,
    c(
      "GSE163973",
      "GSE181316"
    )
  )
)

for (
  metric in plot_metrics
) {
  for (
    cohort in colnames(
      effect_matrix
    )
  ) {
    x <- effect_plot_data[
      effect_plot_data$metric ==
        metric &
      effect_plot_data$cohort ==
        cohort,
      ,
      drop = FALSE
    ]

    effect_matrix[
      metric,
      cohort
    ] <- x$hedges_g
  }
}

pdf(
  file.path(
    FIG_DIR,
    "Figure_core_cohort_effect_sizes.pdf"
  ),
  width = 10,
  height = 7
)

matplot(
  x = seq_len(
    nrow(
      effect_matrix
    )
  ),
  y = effect_matrix,
  type = "b",
  pch = c(
    16,
    1
  ),
  lty = 1,
  xaxt = "n",
  xlab = "",
  ylab =
    "Hedges g: keloid minus normal scar",
  main =
    "Core cohort effect sizes"
)

axis(
  1,
  at = seq_len(
    nrow(
      effect_matrix
    )
  ),
  labels = rownames(
    effect_matrix
  ),
  las = 2,
  cex.axis = 0.8
)

abline(
  h = 0,
  lty = 2
)

legend(
  "topright",
  legend = colnames(
    effect_matrix
  ),
  pch = c(
    16,
    1
  ),
  lty = 1,
  bty = "n"
)

dev.off()

# 3. 早期—晚期修复失败模式
keloid_late <- primary_late$equal_weight_mean_difference
keloid_early <- early_exploratory$equal_weight_mean_difference

vu_late <- stage9[
  stage9$analysis_set ==
    "FB1_FB2_combined_primary" &
    stage9$metric ==
    "late_remodeling_state",
  "mean_difference"
]

vu_early <- stage9[
  stage9$analysis_set ==
    "FB1_FB2_combined_primary" &
    stage9$metric ==
    "early_state_persistence",
  "mean_difference"
]

failure_mode_data <- data.frame(
  disease = c(
    "Keloid_core_integrated",
    "Venous_ulcer"
  ),
  early_state_effect = c(
    keloid_early,
    vu_early
  ),
  late_remodeling_effect = c(
    keloid_late,
    vu_late
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  failure_mode_data,
  file.path(
    REPORT_DIR,
    "10_failure_mode_contrast_source_data.csv"
  )
)

pdf(
  file.path(
    FIG_DIR,
    "Figure_keloid_vs_venous_ulcer_failure_modes.pdf"
  ),
  width = 8,
  height = 7
)

plot(
  failure_mode_data$early_state_effect,
  failure_mode_data$late_remodeling_effect,
  pch = c(
    16,
    1
  ),
  cex = 1.4,
  xlab =
    "Early-state persistence effect",
  ylab =
    "Late-remodeling effect",
  main =
    "Contrasting chronic wound failure modes"
)

text(
  failure_mode_data$early_state_effect,
  failure_mode_data$late_remodeling_effect,
  labels =
    failure_mode_data$disease,
  pos = 3,
  cex = 0.8
)

abline(
  h = 0,
  v = 0,
  lty = 2
)

abline(
  a = 0,
  b = 1,
  lty = 3
)

dev.off()

# ----------------------------- 保存对象 ----------------------------------------
saveRDS(
  list(
    project_gate =
      project_gate,
    interpretation =
      interpretation,
    input_registry =
      input_registry,
    core_patient_data =
      core_patient_data,
    cohort_effects =
      cohort_effects,
    integrated_effects =
      integrated_effects,
    support_evidence =
      support_evidence,
    evidence_hierarchy =
      evidence_hierarchy,
    claim_boundary =
      claim_boundary,
    final_criteria =
      final_criteria,
    failure_mode_data =
      failure_mode_data
  ),
  file = file.path(
    OBJECT_DIR,
    "Keloid_WoundAge_final_integrated_evidence.rds"
  ),
  compress = "gzip"
)

# ----------------------------- 最终结论 ----------------------------------------
decision_lines <- c(
  "第十一阶段：核心队列精确整合与最终证据锁定",
  "============================================================",
  paste0(
    "运行时间：",
    format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  paste0(
    "最终闸门：",
    project_gate
  ),
  "",
  "总体解释：",
  interpretation,
  "",
  "核心主要终点：late-remodeling state",
  paste0(
    "- GSE163973 mean difference = ",
    round(
      cohort_effects$mean_difference[
        cohort_effects$cohort ==
          "GSE163973" &
        cohort_effects$metric ==
          "late_remodeling_state"
      ],
      4
    )
  ),
  paste0(
    "- GSE181316 mean difference = ",
    round(
      cohort_effects$mean_difference[
        cohort_effects$cohort ==
          "GSE181316" &
        cohort_effects$metric ==
          "late_remodeling_state"
      ],
      4
    )
  ),
  paste0(
    "- Equal-weight integrated mean difference = ",
    round(
      primary_late$equal_weight_mean_difference,
      4
    )
  ),
  paste0(
    "- Exact stratified two-sided P = ",
    format(
      primary_late$exact_stratified_two_sided_p,
      digits = 4
    )
  ),
  paste0(
    "- Exact stratified directional P = ",
    format(
      primary_late$exact_stratified_directional_p,
      digits = 4
    )
  ),
  paste0(
    "- Joint permutations = ",
    primary_late$n_exact_joint_permutations
  ),
  paste0(
    "- Both cohorts complete separation = ",
    primary_late$all_cohorts_complete_separation
  ),
  "",
  "核心关键次要终点：ordinal wound-state position",
  paste0(
    "- Integrated mean difference = ",
    round(
      primary_ordinal$equal_weight_mean_difference,
      4
    )
  ),
  paste0(
    "- Exact stratified two-sided P = ",
    format(
      primary_ordinal$exact_stratified_two_sided_p,
      digits = 4
    )
  ),
  "",
  "核心支持指标：",
  paste0(
    "- D30-like exact P = ",
    format(
      support_d30$exact_stratified_two_sided_p,
      digits = 4
    )
  ),
  paste0(
    "- Wound activation exact P = ",
    format(
      support_activation$exact_stratified_two_sided_p,
      digits = 4
    )
  ),
  paste0(
    "- Skin-like reduction exact P = ",
    format(
      support_skin$exact_stratified_two_sided_p,
      digits = 4
    )
  ),
  "",
  "不稳定探索指标：",
  paste0(
    "- Early-state persistence integrated exact P = ",
    format(
      early_exploratory$exact_stratified_two_sided_p,
      digits = 4
    )
  ),
  paste0(
    "- Off-trajectory integrated exact P = ",
    format(
      offtrajectory_exploratory$exact_stratified_two_sided_p,
      digits = 4
    )
  ),
  "",
  "最终证据定位：",
  "- GSE163973与GSE181316构成主要患者级确认性证据。",
  "- GSE220300是第三队列支持，不作为正式独立复制。",
  "- GSE241124验证固定正常伤口签名在独立空间模态中的序数效度。",
  "- GSE265972用于比较不同慢性伤口修复失败模式。",
  "- GSE181297 Visium为组织切片表达支持；其scRNA高特异性结果未复制。",
  "",
  "必须保留的结论边界：",
  "- 不称为精确分子伤口日龄。",
  "- 不声称所有瘢痕疙瘩完全一致。",
  "- 不声称临床活动度与评分成正比。",
  "- 不声称存在稳定中心—周边梯度。",
  "- 不声称GSE181297具有精确空间定位证据。",
  "- 不声称因果机制或术后复发预测。",
  "- early-state persistence和off-trajectory不能作为核心结论。"
)

safe_write_lines(
  decision_lines,
  file.path(
    REPORT_DIR,
    "11_STAGE11_FINAL_DECISION.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    REPORT_DIR,
    "12_SESSION_INFO.txt"
  )
)

safe_write_lines(
  c(
    paste0(
      "STAGE11_COMPLETED=",
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
      "CORE_CONFIRMED=",
      core_confirmed
    ),
    paste0(
      "SUPPORTIVE_KELOID_POSITIVE=",
      supportive_keloid_positive
    ),
    paste0(
      "NORMAL_SPATIAL_VALID=",
      normal_spatial_valid
    ),
    paste0(
      "DISEASE_CONTROL_DISTINCT=",
      disease_control_distinct
    ),
    paste0(
      "STAGE10_CORRECTLY_BOUNDED=",
      stage10_correctly_bounded
    )
  ),
  file.path(
    REPORT_DIR,
    "STAGE11_COMPLETED.txt"
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
  "第十一阶段完成"
)

cat(
  "\n============================================================\n"
)

cat(
  "第十一阶段运行完成。\n"
)

cat(
  "最终整合对象保存在：",
  file.path(
    OBJECT_DIR,
    "Keloid_WoundAge_final_integrated_evidence.rds"
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
