# =============================================================================
# 项目：瘢痕疙瘩“晚期重塑状态持续与伤口状态终止失败”
# 第十二阶段：跨核心队列一致性基因、通路与成纤维细胞状态解释（V1.2最终修正版）（V1.1修正版）
#
# 根目录：
#   C:/Users/33652/Desktop/LZZ/1、yssj
#
# 本阶段定位：
#   - 这是最后一个生物学解释分析；
#   - 不再加入新的公共数据集；
#   - 不改变第四阶段冻结的Skin、Wound1、Wound7、Wound30签名；
#   - 不重新优化主要终点；
#   - GSE163973为发现队列；
#   - GSE181316为预设终点独立复制队列；
#   - GSE220300只作第三队列方向支持；
#   - 通路和已知成纤维细胞状态比较均为机制解释，不是因果证明。
#
# 主要任务：
#   1. 提取两个核心患者级Fibroblast pseudobulk；
#   2. 分别做TMM标准化和log2-CPM；
#   3. 对共同可测基因进行两队列分层精确置换（20×20=400种）；
#   4. 解析四个冻结伤口阶段签名的基因贡献；
#   5. 固定定义：
#        persistent late-remodeling genes
#        = Wound30稳定签名中，两个核心队列均为keloid上调的基因；
#        loss-of-homeostasis genes
#        = Skin稳定签名中，两个核心队列均为keloid下调的基因；
#      不以P值作为选基因条件；
#   6. 进行固定阶段签名GSEA、全基因Hallmark/GO BP/Reactome GSEA；
#   7. 对两个共识模块进行ORA；
#   8. 与透明预定义的成纤维细胞状态标志基因面板比较；
#   9. 生成自动检查包。
#
# 统计边界：
#   - 每个核心队列均为3例keloid vs 3例normal scar；
#   - 基因层面P值为探索性，并进行BH校正；
#   - 共识基因依据“冻结签名成员资格 + 两队列方向一致”定义，
#     不依据显著性筛选；
#   - GSEA排序使用两个队列各自logFC的秩正态化后等权平均，
#     避免不同平台绝对尺度直接合并；
#   - 已知成纤维细胞面板为探索性文献解释工具，不能作为新分型模型。
#   - V1.2修正含NA逻辑索引导致的共识基因空白行与计数偏差；
#   - V1.2的富集审计图仅显示FDR<0.05结果。
#   - V1.1修复calculate_cohort_gene_effects中两处矩阵下标缺失右方括号的问题。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANG = "en")
set.seed(20260618)

ROOT_DIR <- "C:/Users/33652/Desktop/LZZ/1、yssj"

STAGE_DIR <- file.path(
  ROOT_DIR,
  "12_STAGE12_GENE_PROGRAM_INTERPRETATION"
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
  "第十二阶段_跨队列一致性基因与通路解释检查包.zip"
)

PACKAGE_TARGZ <- file.path(
  ROOT_DIR,
  "第十二阶段_跨队列一致性基因与通路解释检查包.tar.gz"
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
  "00_stage12_log.txt"
)

WARNING_FILE <- file.path(
  REPORT_DIR,
  "00_stage12_warnings.txt"
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
install_cran_if_missing <- function(pkg) {
  if (!requireNamespace(
    pkg,
    quietly = TRUE
  )) {
    install.packages(
      pkg,
      repos = "https://cloud.r-project.org"
    )
  }
}

install_bioc_if_missing <- function(pkg) {
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
    pkg,
    quietly = TRUE
  )) {
    BiocManager::install(
      pkg,
      ask = FALSE,
      update = FALSE
    )
  }
}

install_cran_if_missing(
  "Matrix"
)

install_cran_if_missing(
  "msigdbr"
)

install_bioc_if_missing(
  "edgeR"
)

install_bioc_if_missing(
  "clusterProfiler"
)

install_bioc_if_missing(
  "org.Hs.eg.db"
)

core_packages <- c(
  "Matrix",
  "edgeR"
)

missing_core_packages <- core_packages[
  !vapply(
    core_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(
  missing_core_packages
)) {
  stop(
    "无法加载核心分析包：",
    paste(
      missing_core_packages,
      collapse = ", "
    )
  )
}

pathway_packages_available <- all(
  vapply(
    c(
      "msigdbr",
      "clusterProfiler"
    ),
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
)

# ----------------------------- 文件定位 ----------------------------------------
find_file <- function(
  target_basename,
  preferred_fragment = NULL,
  required = TRUE
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

  if (
    !is.null(
      preferred_fragment
    ) &&
      length(
        candidates
      )
  ) {
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
  }

  if (!length(
    candidates
  )) {
    if (required) {
      stop(
        "未找到必要文件：",
        target_basename
      )
    }

    return(
      NA_character_
    )
  }

  if (length(
    candidates
  ) >
      1L) {
    info <- file.info(
      candidates
    )

    selected <- candidates[
      order(
        info$mtime,
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
  object_role = c(
    "frozen_reference_model",
    "discovery_core_cohort",
    "validation_core_cohort",
    "supportive_third_cohort",
    "final_integration_optional"
  ),
  target_basename = c(
    "GSE241132_fixed_wound_state_reference.rds",
    "GSE163973_fixed_wound_state_mapping.rds",
    "GSE181316_fixed_wound_state_replication.rds",
    "GSE220300_activity_region_wound_state_support.rds",
    "Keloid_WoundAge_final_integrated_evidence.rds"
  ),
  preferred_fragment = c(
    "04_STAGE4_WOUND_STATE_MODEL",
    "05_STAGE5_GSE163973_MAPPING",
    "06B_STAGE6B_GSE181316_REPLICATION",
    "07_STAGE7_GSE220300_ACTIVITY_REGION",
    "11_STAGE11_FINAL_INTEGRATION"
  ),
  required = c(
    TRUE,
    TRUE,
    TRUE,
    FALSE,
    FALSE
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
    find_file(
      input_registry$target_basename[i],
      input_registry$preferred_fragment[i],
      required =
        input_registry$required[i]
    )
  },
  character(1)
)

safe_write_csv(
  input_registry,
  file.path(
    REPORT_DIR,
    "01_input_object_registry.csv"
  )
)

get_input_path <- function(role) {
  input_registry$resolved_path[
    input_registry$object_role ==
      role
  ][1]
}

# ----------------------------- 基础数据函数 ------------------------------------
canonicalize_count_matrix <- function(counts) {
  if (is.null(
    rownames(
      counts
    )
  )) {
    stop(
      "计数矩阵缺少基因名。"
    )
  }

  gene <- toupper(
    trimws(
      as.character(
        rownames(
          counts
        )
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

  if (inherits(
    counts,
    "sparseMatrix"
  )) {
    counts <- as.matrix(
      counts
    )
  }

  storage.mode(
    counts
  ) <- "double"

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
        "keloid_lesion",
        "active_lesion",
        "inactive_lesion"
      )
  ] <- "keloid"

  y[
    y %in%
      c(
        "ns",
        "normal_scar",
        "scar_reference",
        "mature_scar"
      )
  ] <- "normal_scar"

  y
}

prepare_core_cohort <- function(
  counts,
  sample_info,
  sample_id_column,
  group_column,
  cohort_name
) {
  counts <- canonicalize_count_matrix(
    counts
  )

  sample_ids <- as.character(
    sample_info[[sample_id_column]]
  )

  groups <- normalize_group(
    sample_info[[group_column]]
  )

  keep_samples <- groups %in%
    c(
      "keloid",
      "normal_scar"
    )

  sample_info <- sample_info[
    keep_samples,
    ,
    drop = FALSE
  ]

  sample_ids <- sample_ids[
    keep_samples
  ]

  groups <- groups[
    keep_samples
  ]

  count_index <- match(
    sample_ids,
    colnames(
      counts
    )
  )

  if (anyNA(
    count_index
  )) {
    stop(
      cohort_name,
      "样本信息无法匹配计数矩阵列：",
      paste(
        sample_ids[
          is.na(
            count_index
          )
        ],
        collapse = ", "
      )
    )
  }

  counts <- counts[
    ,
    count_index,
    drop = FALSE
  ]

  colnames(
    counts
  ) <- sample_ids

  sample_info_core <- data.frame(
    cohort =
      cohort_name,
    sample_id =
      sample_ids,
    group =
      groups,
    stringsAsFactors = FALSE
  )

  if (
    sum(
      groups ==
        "keloid"
    ) !=
      3L ||
      sum(
        groups ==
          "normal_scar"
      ) !=
      3L
  ) {
    stop(
      cohort_name,
      "不是3例keloid与3例normal scar。实际：",
      paste(
        names(
          table(
            groups
          )
        ),
        as.integer(
          table(
            groups
          )
        ),
        sep = "=",
        collapse = ", "
      )
    )
  }

  dge <- edgeR::DGEList(
    counts =
      counts,
    group =
      groups
  )

  cpm_raw <- edgeR::cpm(
    dge,
    log = FALSE
  )

  keep_genes <- rowSums(
    cpm_raw >=
      1
  ) >=
    2

  if (
    sum(
      keep_genes
    ) <
      5000L
  ) {
    warn_msg(
      cohort_name,
      "CPM过滤后基因少于5000：",
      sum(
        keep_genes
      )
    )
  }

  dge <- dge[
    keep_genes,
    ,
    keep.lib.sizes = FALSE
  ]

  dge <- edgeR::calcNormFactors(
    dge,
    method = "TMM"
  )

  log_cpm <- edgeR::cpm(
    dge,
    log = TRUE,
    prior.count = 2
  )

  list(
    counts =
      dge$counts,
    log_cpm =
      log_cpm,
    sample_info =
      sample_info_core,
    norm_factors =
      dge$samples
  )
}

prepare_supportive_cohort <- function(
  counts,
  mapping,
  cohort_name
) {
  counts <- canonicalize_count_matrix(
    counts
  )

  required_columns <- c(
    "sample_key",
    "aggregate_state"
  )

  if (!all(
    required_columns %in%
      names(
        mapping
      )
  )) {
    stop(
      cohort_name,
      "支持队列mapping缺少字段。"
    )
  }

  sample_ids <- as.character(
    mapping$sample_key
  )

  group <- ifelse(
    mapping$aggregate_state %in%
      c(
        "active_lesion",
        "inactive_lesion"
      ),
    "keloid",
    "normal_scar"
  )

  index <- match(
    sample_ids,
    colnames(
      counts
    )
  )

  if (anyNA(
    index
  )) {
    stop(
      cohort_name,
      "mapping无法匹配aggregate_counts。"
    )
  }

  counts <- counts[
    ,
    index,
    drop = FALSE
  ]

  colnames(
    counts
  ) <- sample_ids

  dge <- edgeR::DGEList(
    counts =
      counts,
    group =
      group
  )

  cpm_raw <- edgeR::cpm(
    dge,
    log = FALSE
  )

  keep_genes <- rowSums(
    cpm_raw >=
      1
  ) >=
    2

  dge <- dge[
    keep_genes,
    ,
    keep.lib.sizes = FALSE
  ]

  dge <- edgeR::calcNormFactors(
    dge,
    method = "TMM"
  )

  log_cpm <- edgeR::cpm(
    dge,
    log = TRUE,
    prior.count = 2
  )

  list(
    log_cpm =
      log_cpm,
    sample_info = data.frame(
      cohort =
        cohort_name,
      sample_id =
        sample_ids,
      group =
        group,
      stringsAsFactors = FALSE
    )
  )
}

get_stable_signatures <- function(
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
      "稳定Fibroblast签名基因数异常：",
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

# ----------------------------- 基因效应函数 ------------------------------------
build_contrast_matrix_3v3 <- function() {
  combinations <- utils::combn(
    6L,
    3L
  )

  contrast <- matrix(
    -1 / 3,
    nrow = 6L,
    ncol = ncol(
      combinations
    )
  )

  for (
    j in seq_len(
      ncol(
        combinations
      )
    )
  ) {
    contrast[
      combinations[, j],
      j
    ] <- 1 / 3
  }

  contrast
}

vectorized_cliffs_delta <- function(
  expression,
  group
) {
  keloid_index <- which(
    group ==
      "keloid"
  )

  control_index <- which(
    group ==
      "normal_scar"
  )

  score <- rep(
    0,
    nrow(
      expression
    )
  )

  for (
    i in keloid_index
  ) {
    for (
      j in control_index
    ) {
      score <- score +
        sign(
          expression[, i] -
            expression[, j]
        )
    }
  }

  score /
    (
      length(
        keloid_index
      ) *
        length(
          control_index
        )
    )
}

vectorized_hedges_g <- function(
  expression,
  group
) {
  x <- expression[
    ,
    group ==
      "keloid",
    drop = FALSE
  ]

  y <- expression[
    ,
    group ==
      "normal_scar",
    drop = FALSE
  ]

  n1 <- ncol(
    x
  )

  n2 <- ncol(
    y
  )

  mean_difference <- rowMeans(
    x
  ) -
    rowMeans(
      y
    )

  var_x <- apply(
    x,
    1,
    stats::var
  )

  var_y <- apply(
    y,
    1,
    stats::var
  )

  degrees_freedom <- n1 +
    n2 -
    2

  pooled_sd <- sqrt(
    (
      (n1 -
        1) *
        var_x +
        (n2 -
          1) *
        var_y
    ) /
      degrees_freedom
  )

  cohen_d <- mean_difference /
    pooled_sd

  correction <- 1 -
    3 /
    (
      4 *
        degrees_freedom -
        1
    )

  result <- correction *
    cohen_d

  result[
    !is.finite(
      result
    )
  ] <- NA_real_

  result
}

calculate_cohort_gene_effects <- function(
  expression,
  sample_info,
  cohort_name
) {
  group <- sample_info$group

  if (
    ncol(
      expression
    ) !=
      6L ||
      sum(
        group ==
          "keloid"
      ) !=
      3L ||
      sum(
        group ==
          "normal_scar"
      ) !=
      3L
  ) {
    stop(
      cohort_name,
      "不符合3 vs 3基因效应分析结构。"
    )
  }

  observed <- rowMeans(
    expression[
      ,
      group ==
        "keloid",
      drop = FALSE
    ]
  ) -
    rowMeans(
      expression[
        ,
        group ==
          "normal_scar",
      drop = FALSE
    ]
  )

  contrast <- build_contrast_matrix_3v3()

  null_matrix <- expression %*%
    contrast

  exact_two_sided_p <- rowMeans(
    abs(
      null_matrix
    ) >=
      abs(
        observed
      ) -
      1e-12
  )

  exact_positive_p <- rowMeans(
    null_matrix >=
      observed -
      1e-12
  )

  exact_negative_p <- rowMeans(
    null_matrix <=
      observed +
      1e-12
  )

  data.frame(
    gene =
      rownames(
        expression
      ),
    cohort =
      cohort_name,
    log2_cpm_mean_keloid =
      rowMeans(
        expression[
          ,
          group ==
            "keloid",
          drop = FALSE
        ]
      ),
    log2_cpm_mean_normal_scar =
      rowMeans(
        expression[
          ,
          group ==
            "normal_scar",
          drop = FALSE
        ]
      ),
    mean_difference =
      observed,
    hedges_g =
      vectorized_hedges_g(
        expression,
        group
      ),
    cliffs_delta =
      vectorized_cliffs_delta(
        expression,
        group
      ),
    exact_two_sided_p =
      exact_two_sided_p,
    exact_positive_p =
      exact_positive_p,
    exact_negative_p =
      exact_negative_p,
    stringsAsFactors = FALSE
  )
}

calculate_stratified_gene_integration <- function(
  expression_1,
  sample_info_1,
  expression_2,
  sample_info_2,
  cohort_name_1,
  cohort_name_2
) {
  common_genes <- intersect(
    rownames(
      expression_1
    ),
    rownames(
      expression_2
    )
  )

  if (
    length(
      common_genes
    ) <
      5000L
  ) {
    warn_msg(
      "两个核心队列共同可测基因少于5000：",
      length(
        common_genes
      )
    )
  }

  expression_1 <- expression_1[
    common_genes,
    ,
    drop = FALSE
  ]

  expression_2 <- expression_2[
    common_genes,
    ,
    drop = FALSE
  ]

  effect_1 <- calculate_cohort_gene_effects(
    expression_1,
    sample_info_1,
    cohort_name_1
  )

  effect_2 <- calculate_cohort_gene_effects(
    expression_2,
    sample_info_2,
    cohort_name_2
  )

  effect_1 <- effect_1[
    match(
      common_genes,
      effect_1$gene
    ),
    ,
    drop = FALSE
  ]

  effect_2 <- effect_2[
    match(
      common_genes,
      effect_2$gene
    ),
    ,
    drop = FALSE
  ]

  contrast <- build_contrast_matrix_3v3()

  null_1 <- expression_1 %*%
    contrast

  null_2 <- expression_2 %*%
    contrast

  observed_1 <- effect_1$mean_difference
  observed_2 <- effect_2$mean_difference

  observed_integrated <- (
    observed_1 +
      observed_2
  ) /
    2

  count_two_sided <- numeric(
    length(
      common_genes
    )
  )

  count_positive <- numeric(
    length(
      common_genes
    )
  )

  count_negative <- numeric(
    length(
      common_genes
    )
  )

  for (
    j in seq_len(
      ncol(
        null_2
      )
    )
  ) {
    joint_null <- (
      null_1 +
        null_2[, j]
    ) /
      2

    count_two_sided <- count_two_sided +
      rowSums(
        abs(
          joint_null
        ) >=
          abs(
            observed_integrated
          ) -
          1e-12
      )

    count_positive <- count_positive +
      rowSums(
        joint_null >=
          observed_integrated -
          1e-12
      )

    count_negative <- count_negative +
      rowSums(
        joint_null <=
          observed_integrated +
          1e-12
      )
  }

  n_joint <- ncol(
    null_1
  ) *
    ncol(
      null_2
    )

  rank_normal_score <- function(x) {
    stats::qnorm(
      (
        rank(
          x,
          ties.method = "average"
        ) -
          0.5
      ) /
        length(
          x
        )
    )
  }

  integrated_rank_score <- (
    rank_normal_score(
      observed_1
    ) +
      rank_normal_score(
        observed_2
      )
  ) /
    2

  result <- data.frame(
    gene =
      common_genes,
    mean_difference_GSE163973 =
      observed_1,
    hedges_g_GSE163973 =
      effect_1$hedges_g,
    cliffs_delta_GSE163973 =
      effect_1$cliffs_delta,
    exact_p_GSE163973 =
      effect_1$exact_two_sided_p,
    mean_difference_GSE181316 =
      observed_2,
    hedges_g_GSE181316 =
      effect_2$hedges_g,
    cliffs_delta_GSE181316 =
      effect_2$cliffs_delta,
    exact_p_GSE181316 =
      effect_2$exact_two_sided_p,
    equal_weight_mean_difference =
      observed_integrated,
    integrated_rank_score =
      integrated_rank_score,
    exact_stratified_two_sided_p =
      count_two_sided /
      n_joint,
    exact_stratified_positive_p =
      count_positive /
      n_joint,
    exact_stratified_negative_p =
      count_negative /
      n_joint,
    concordance = ifelse(
      observed_1 >
        0 &
        observed_2 >
        0,
      "positive_both",
      ifelse(
        observed_1 <
          0 &
          observed_2 <
          0,
        "negative_both",
        ifelse(
          observed_1 ==
            0 &
            observed_2 ==
            0,
          "zero_both",
          "discordant"
        )
      )
    ),
    stringsAsFactors = FALSE
  )

  result$BH_FDR_stratified <-
    stats::p.adjust(
      result$exact_stratified_two_sided_p,
      method = "BH"
    )

  result
}

# ----------------------------- 模块评分与精确整合 ------------------------------
score_module <- function(
  expression,
  genes
) {
  genes_present <- intersect(
    genes,
    rownames(
      expression
    )
  )

  if (
    length(
      genes_present
    ) <
      2L
  ) {
    return(
      rep(
        NA_real_,
        ncol(
          expression
        )
      )
    )
  }

  gene_matrix <- expression[
    genes_present,
    ,
    drop = FALSE
  ]

  gene_z <- t(
    scale(
      t(
        gene_matrix
      )
    )
  )

  gene_z[
    !is.finite(
      gene_z
    )
  ] <- 0

  colMeans(
    gene_z
  )
}

exact_3v3_score <- function(
  score,
  group
) {
  if (
    sum(
      group ==
        "keloid"
    ) !=
      3L ||
      sum(
        group ==
          "normal_scar"
      ) !=
      3L
  ) {
    stop(
      "模块评分精确检验不是3 vs 3结构。"
    )
  }

  observed <- mean(
    score[
      group ==
        "keloid"
    ]
  ) -
    mean(
      score[
        group ==
          "normal_scar"
    ]
  )

  contrast <- build_contrast_matrix_3v3()

  null <- as.numeric(
    matrix(
      score,
      nrow = 1
    ) %*%
      contrast
  )

  list(
    observed =
      observed,
    null =
      null,
    exact_two_sided_p =
      mean(
        abs(
          null
        ) >=
          abs(
            observed
          ) -
          1e-12
      ),
    exact_positive_p =
      mean(
        null >=
          observed -
          1e-12
      ),
    cliffs_delta =
      {
        x <- score[
          group ==
            "keloid"
        ]

        y <- score[
          group ==
            "normal_scar"
        ]

        mean(
          outer(
            x,
            y,
            "-"
          ) >
            0
        ) -
          mean(
            outer(
              x,
              y,
              "-"
            ) <
              0
          )
      }
  )
}

integrate_module_scores <- function(
  score_1,
  group_1,
  score_2,
  group_2
) {
  test_1 <- exact_3v3_score(
    score_1,
    group_1
  )

  test_2 <- exact_3v3_score(
    score_2,
    group_2
  )

  observed <- (
    test_1$observed +
      test_2$observed
  ) /
    2

  joint_null <- as.vector(
    outer(
      test_1$null,
      test_2$null,
      "+"
    ) /
      2
  )

  list(
    effect_GSE163973 =
      test_1$observed,
    cliff_GSE163973 =
      test_1$cliffs_delta,
    effect_GSE181316 =
      test_2$observed,
    cliff_GSE181316 =
      test_2$cliffs_delta,
    integrated_effect =
      observed,
    exact_stratified_two_sided_p =
      mean(
        abs(
          joint_null
        ) >=
          abs(
            observed
          ) -
          1e-12
      ),
    exact_stratified_positive_p =
      mean(
        joint_null >=
          observed -
          1e-12
      )
  )
}

# ----------------------------- GSEA/ORA函数 ------------------------------------
get_msigdb_term2gene <- function(
  collection,
  subcollection = NULL
) {
  if (!pathway_packages_available) {
    stop(
      "msigdbr或clusterProfiler不可用。"
    )
  }

  fml <- names(
    formals(
      msigdbr::msigdbr
    )
  )

  args <- list(
    species = "Homo sapiens"
  )

  if ("collection" %in%
      fml) {
    args$collection <- collection
  } else {
    args$category <- collection
  }

  if (!is.null(
    subcollection
  )) {
    if ("subcollection" %in%
        fml) {
      args$subcollection <-
        subcollection
    } else if ("subcategory" %in%
               fml) {
      args$subcategory <-
        subcollection
    }
  }

  msig <- tryCatch(
    do.call(
      msigdbr::msigdbr,
      args
    ),
    error = function(e) {
      warn_msg(
        "msigdbr读取失败：",
        collection,
        "/",
        ifelse(
          is.null(
            subcollection
          ),
          "ALL",
          subcollection
        ),
        "|",
        conditionMessage(e)
      )

      NULL
    }
  )

  if (is.null(
    msig
  )) {
    return(
      data.frame()
    )
  }

  gene_column <- if (
    "gene_symbol" %in%
      names(
        msig
      )
  ) {
    "gene_symbol"
  } else if (
    "human_gene_symbol" %in%
      names(
        msig
      )
  ) {
    "human_gene_symbol"
  } else {
    stop(
      "msigdbr结果缺少gene symbol字段。"
    )
  }

  set_column <- if (
    "gs_name" %in%
      names(
        msig
      )
  ) {
    "gs_name"
  } else {
    stop(
      "msigdbr结果缺少gs_name字段。"
    )
  }

  result <- data.frame(
    term =
      as.character(
        msig[[set_column]]
      ),
    gene =
      toupper(
        as.character(
          msig[[gene_column]]
        )
      ),
    stringsAsFactors = FALSE
  )

  result <- result[
    !is.na(
      result$term
    ) &
      nzchar(
        result$term
      ) &
      !is.na(
        result$gene
      ) &
      nzchar(
        result$gene
      ),
    ,
    drop = FALSE
  ]

  unique(
    result
  )
}

run_gsea <- function(
  gene_score,
  term2gene,
  analysis_name
) {
  if (
    !pathway_packages_available ||
      !nrow(
        term2gene
      )
  ) {
    return(
      data.frame()
    )
  }

  gene_score <- gene_score[
    is.finite(
      gene_score
    )
  ]

  gene_score <- sort(
    gene_score,
    decreasing = TRUE
  )

  if (
    anyDuplicated(
      names(
        gene_score
      )
    )
  ) {
    gene_score <- gene_score[
      !duplicated(
        names(
          gene_score
        )
      )
    ]
  }

  # 给并列值增加极小确定性偏移，不改变实质排序。
  gene_score <- gene_score +
    seq_along(
      gene_score
    ) *
    1e-12

  gsea_object <- tryCatch(
    clusterProfiler::GSEA(
      geneList =
        gene_score,
      TERM2GENE =
        term2gene,
      minGSSize =
        10,
      maxGSSize =
        500,
      pvalueCutoff =
        1,
      pAdjustMethod =
        "BH",
      verbose =
        FALSE,
      by =
        "fgsea"
    ),
    error = function(e) {
      warn_msg(
        analysis_name,
        "GSEA失败：",
        conditionMessage(e)
      )

      NULL
    }
  )

  if (is.null(
    gsea_object
  )) {
    return(
      data.frame()
    )
  }

  result <- as.data.frame(
    gsea_object
  )

  if (nrow(
    result
  )) {
    result$analysis <-
      analysis_name
  }

  result
}

run_ora <- function(
  genes,
  universe,
  term2gene,
  analysis_name
) {
  genes <- intersect(
    unique(
      toupper(
        genes
      )
    ),
    universe
  )

  if (
    !pathway_packages_available ||
      !nrow(
        term2gene
      ) ||
      length(
        genes
      ) <
      3L
  ) {
    return(
      data.frame()
    )
  }

  ora_object <- tryCatch(
    clusterProfiler::enricher(
      gene =
        genes,
      universe =
        universe,
      TERM2GENE =
        term2gene,
      minGSSize =
        5,
      maxGSSize =
        500,
      pvalueCutoff =
        1,
      pAdjustMethod =
        "BH",
      qvalueCutoff =
        1
    ),
    error = function(e) {
      warn_msg(
        analysis_name,
        "ORA失败：",
        conditionMessage(e)
      )

      NULL
    }
  )

  if (is.null(
    ora_object
  )) {
    return(
      data.frame()
    )
  }

  result <- as.data.frame(
    ora_object
  )

  if (nrow(
    result
  )) {
    result$analysis <-
      analysis_name
  }

  result
}

safe_msigdb_collection <- function(
  collection,
  subcollection = NULL,
  label
) {
  tryCatch(
    get_msigdb_term2gene(
      collection,
      subcollection
    ),
    error = function(e) {
      warn_msg(
        label,
        "基因集读取失败：",
        conditionMessage(e)
      )

      data.frame()
    }
  )
}

# ----------------------------- 图形函数 ----------------------------------------
plot_top_enrichment <- function(
  result,
  file_path,
  title,
  top_n = 15L,
  fdr_cutoff = 0.05
) {
  create_empty_enrichment_page <- function(
    message
  ) {
    pdf(
      file_path,
      width = 10,
      height = 6
    )

    plot.new()

    title(
      main =
        title
    )

    text(
      0.5,
      0.5,
      labels =
        message,
      cex = 1.1
    )

    dev.off()
  }

  if (
    is.null(
      result
    ) ||
      !nrow(
        result
      )
  ) {
    create_empty_enrichment_page(
      "No enrichment result was generated."
    )

    return(
      invisible(
        FALSE
      )
    )
  }

  if (!"p.adjust" %in%
      names(
        result
      )) {
    create_empty_enrichment_page(
      "Adjusted P values were unavailable."
    )

    return(
      invisible(
        FALSE
      )
    )
  }

  result <- result[
    is.finite(
      result$p.adjust
    ) &
      result$p.adjust <
      fdr_cutoff,
    ,
    drop = FALSE
  ]

  if (!nrow(
    result
  )) {
    create_empty_enrichment_page(
      paste0(
        "No FDR-significant term (FDR < ",
        fdr_cutoff,
        ")."
      )
    )

    return(
      invisible(
        FALSE
      )
    )
  }

  if (
    "NES" %in%
      names(
        result
      )
  ) {
    result <- result[
      is.finite(
        result$NES
      ),
      ,
      drop = FALSE
    ]

    result <- result[
      order(
        result$p.adjust,
        -abs(
          result$NES
        )
      ),
      ,
      drop = FALSE
    ]

    result <- head(
      result,
      top_n
    )

    pdf(
      file_path,
      width = 10,
      height = max(
        6,
        0.35 *
          nrow(
            result
          ) +
          2
      )
    )

    par(
      mar = c(
        5,
        14,
        3,
        2
      )
    )

    barplot(
      rev(
        result$NES
      ),
      horiz = TRUE,
      names.arg =
        rev(
          result$Description
        ),
      las = 1,
      xlab =
        "Normalized enrichment score",
      main =
        title
    )

    abline(
      v = 0,
      lty = 2
    )

    dev.off()
  } else {
    result <- result[
      order(
        result$p.adjust,
        result$pvalue
      ),
      ,
      drop = FALSE
    ]

    result <- head(
      result,
      top_n
    )

    score <- -log10(
      pmax(
        result$p.adjust,
        1e-300
      )
    )

    pdf(
      file_path,
      width = 10,
      height = max(
        6,
        0.35 *
          nrow(
            result
          ) +
          2
      )
    )

    par(
      mar = c(
        5,
        14,
        3,
        2
      )
    )

    barplot(
      rev(
        score
      ),
      horiz = TRUE,
      names.arg =
        rev(
          result$Description
        ),
      las = 1,
      xlab =
        "-log10 adjusted P",
      main =
        title
    )

    dev.off()
  }

  invisible(
    TRUE
  )
}

# ----------------------------- 读取对象 ----------------------------------------
log_msg(
  "第十二阶段开始"
)

reference_model <- readRDS(
  get_input_path(
    "frozen_reference_model"
  )
)

stage5_object <- readRDS(
  get_input_path(
    "discovery_core_cohort"
  )
)

stage6_object <- readRDS(
  get_input_path(
    "validation_core_cohort"
  )
)

stage7_path <- get_input_path(
  "supportive_third_cohort"
)

stage7_object <- if (
  is.na(
    stage7_path
  )
) {
  NULL
} else {
  readRDS(
    stage7_path
  )
}

signatures <- get_stable_signatures(
  reference_model
)

signature_summary <- data.frame(
  stage =
    names(
      signatures
    ),
  n_frozen_stable_genes =
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
    "02_frozen_signature_summary.csv"
  )
)

# ----------------------------- 提取核心pseudobulk ------------------------------
if (
  is.null(
    stage5_object$disease_pseudobulk$counts
  ) ||
    is.null(
      stage5_object$disease_pseudobulk$sample_info
    )
) {
  stop(
    "Stage5对象缺少disease_pseudobulk。"
  )
}

stage5_counts_all <-
  stage5_object$disease_pseudobulk$counts

stage5_info_all <-
  stage5_object$disease_pseudobulk$sample_info

stage5_fibro_index <- which(
  stage5_info_all$main_cell_type ==
    "Fibroblast"
)

if (
  length(
    stage5_fibro_index
  ) !=
    6L
) {
  stop(
    "Stage5 Fibroblast pseudobulk不是6个：",
    length(
      stage5_fibro_index
    )
  )
}

stage5_fibro_info <- stage5_info_all[
  stage5_fibro_index,
  ,
  drop = FALSE
]

stage5_fibro_counts <- stage5_counts_all[
  ,
  match(
    stage5_fibro_info$sample_key,
    colnames(
      stage5_counts_all
    )
  ),
  drop = FALSE
]

stage5_fibro_info$analysis_sample_id <-
  stage5_fibro_info$sample_key

if (
  is.null(
    stage6_object$high_specificity$patient_counts
  ) ||
    is.null(
      stage6_object$high_specificity$patient_info
    )
) {
  stop(
    "Stage6B对象缺少high_specificity patient counts/info。"
  )
}

stage6_patient_counts <-
  stage6_object$high_specificity$patient_counts

stage6_patient_info <-
  stage6_object$high_specificity$patient_info

stage5_core <- prepare_core_cohort(
  counts =
    stage5_fibro_counts,
  sample_info =
    stage5_fibro_info,
  sample_id_column =
    "analysis_sample_id",
  group_column =
    "group",
  cohort_name =
    "GSE163973"
)

stage6_core <- prepare_core_cohort(
  counts =
    stage6_patient_counts,
  sample_info =
    stage6_patient_info,
  sample_id_column =
    "patient_id",
  group_column =
    "group",
  cohort_name =
    "GSE181316"
)

core_sample_structure <- rbind(
  stage5_core$sample_info,
  stage6_core$sample_info
)

safe_write_csv(
  core_sample_structure,
  file.path(
    REPORT_DIR,
    "03_core_cohort_sample_structure.csv"
  )
)

normalization_summary <- rbind(
  data.frame(
    cohort =
      "GSE163973",
    sample_id =
      rownames(
        stage5_core$norm_factors
      ),
    library_size =
      stage5_core$norm_factors$lib.size,
    normalization_factor =
      stage5_core$norm_factors$norm.factors,
    stringsAsFactors = FALSE
  ),
  data.frame(
    cohort =
      "GSE181316",
    sample_id =
      rownames(
        stage6_core$norm_factors
      ),
    library_size =
      stage6_core$norm_factors$lib.size,
    normalization_factor =
      stage6_core$norm_factors$norm.factors,
    stringsAsFactors = FALSE
  )
)

safe_write_csv(
  normalization_summary,
  file.path(
    REPORT_DIR,
    "04_TMM_normalization_summary.csv"
  )
)

# ----------------------------- 全基因核心整合 ----------------------------------
log_msg(
  "计算两个核心队列全基因分层精确效应"
)

all_gene_effects <- calculate_stratified_gene_integration(
  expression_1 =
    stage5_core$log_cpm,
  sample_info_1 =
    stage5_core$sample_info,
  expression_2 =
    stage6_core$log_cpm,
  sample_info_2 =
    stage6_core$sample_info,
  cohort_name_1 =
    "GSE163973",
  cohort_name_2 =
    "GSE181316"
)

safe_write_csv(
  all_gene_effects,
  file.path(
    REPORT_DIR,
    "05_all_gene_core_effects_and_exact_integration.csv"
  )
)

common_gene_universe <- all_gene_effects$gene

# ----------------------------- 第三队列基因方向 --------------------------------
stage7_gene_effects <- data.frame()

if (!is.null(
  stage7_object
)) {
  if (
    !is.null(
      stage7_object$high_specificity$aggregate_counts
    ) &&
      !is.null(
        stage7_object$high_specificity$aggregate_mapping
      )
  ) {
    stage7_support <- prepare_supportive_cohort(
      counts =
        stage7_object$high_specificity$aggregate_counts,
      mapping =
        stage7_object$high_specificity$aggregate_mapping,
      cohort_name =
        "GSE220300"
    )

    expression_7 <- stage7_support$log_cpm
    group_7 <- stage7_support$sample_info$group

    stage7_gene_effects <- data.frame(
      gene =
        rownames(
          expression_7
        ),
      mean_difference_GSE220300 =
        rowMeans(
          expression_7[
            ,
            group_7 ==
              "keloid",
            drop = FALSE
          ]
        ) -
        rowMeans(
          expression_7[
            ,
            group_7 ==
              "normal_scar",
            drop = FALSE
          ]
        ),
      cliffs_delta_GSE220300 =
        vectorized_cliffs_delta(
          expression_7,
          group_7
        ),
      stringsAsFactors = FALSE
    )

    safe_write_csv(
      stage7_gene_effects,
      file.path(
        REPORT_DIR,
        "06_GSE220300_supportive_gene_effects.csv"
      )
    )
  } else {
    warn_msg(
      "Stage7对象缺少high_specificity aggregate counts/mapping；跳过第三队列基因方向。"
    )
  }
}

# ----------------------------- 固定签名基因贡献 --------------------------------
signature_gene_table <- do.call(
  rbind,
  lapply(
    names(
      signatures
    ),
    function(stage) {
      data.frame(
        stage =
          stage,
        gene =
          signatures[[stage]],
        expected_direction = if (
          stage ==
            "Skin"
        ) {
          "negative_in_keloid"
        } else {
          "positive_in_keloid"
        },
        stringsAsFactors = FALSE
      )
    }
  )
)

signature_gene_table <- merge(
  signature_gene_table,
  all_gene_effects,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

if (nrow(
  stage7_gene_effects
)) {
  signature_gene_table <- merge(
    signature_gene_table,
    stage7_gene_effects,
    by = "gene",
    all.x = TRUE,
    sort = FALSE
  )
}

signature_gene_table$core_direction_match <- with(
  signature_gene_table,
  ifelse(
    stage ==
      "Skin",
    mean_difference_GSE163973 <
      0 &
      mean_difference_GSE181316 <
      0,
    mean_difference_GSE163973 >
      0 &
      mean_difference_GSE181316 >
      0
  )
)

signature_gene_table$stage7_direction_match <- if (
  "mean_difference_GSE220300" %in%
    names(
      signature_gene_table
    )
) {
  with(
    signature_gene_table,
    ifelse(
      is.na(
        mean_difference_GSE220300
      ),
      NA,
      ifelse(
        stage ==
          "Skin",
        mean_difference_GSE220300 <
          0,
        mean_difference_GSE220300 >
          0
      )
    )
  )
} else {
  NA
}

signature_gene_table$minimum_oriented_core_effect <- with(
  signature_gene_table,
  ifelse(
    stage ==
      "Skin",
    pmin(
      -mean_difference_GSE163973,
      -mean_difference_GSE181316
    ),
    pmin(
      mean_difference_GSE163973,
      mean_difference_GSE181316
    )
  )
)

signature_gene_table <- signature_gene_table[
  order(
    match(
      signature_gene_table$stage,
      c(
        "Skin",
        "Wound1",
        "Wound7",
        "Wound30"
      )
    ),
    -signature_gene_table$minimum_oriented_core_effect
  ),
  ,
  drop = FALSE
]

safe_write_csv(
  signature_gene_table,
  file.path(
    REPORT_DIR,
    "07_fixed_signature_gene_contribution_table.csv"
  )
)

# 共识模块不以P值筛选。
late_consensus_index <- which(
  signature_gene_table$stage ==
    "Wound30" &
    !is.na(
      signature_gene_table$core_direction_match
    ) &
    signature_gene_table$core_direction_match &
    !is.na(
      signature_gene_table$gene
    ) &
    nzchar(
      signature_gene_table$gene
    )
)

skin_loss_consensus_index <- which(
  signature_gene_table$stage ==
    "Skin" &
    !is.na(
      signature_gene_table$core_direction_match
    ) &
    signature_gene_table$core_direction_match &
    !is.na(
      signature_gene_table$gene
    ) &
    nzchar(
      signature_gene_table$gene
    )
)

consensus_late_table <- signature_gene_table[
  late_consensus_index,
  ,
  drop = FALSE
]

consensus_skin_loss_table <- signature_gene_table[
  skin_loss_consensus_index,
  ,
  drop = FALSE
]

consensus_late_genes <- unique(
  as.character(
    consensus_late_table$gene
  )
)

consensus_skin_loss_genes <- unique(
  as.character(
    consensus_skin_loss_table$gene
  )
)

if (
  anyNA(
    consensus_late_genes
  ) ||
    any(
      !nzchar(
        consensus_late_genes
      )
    ) ||
    anyNA(
      consensus_skin_loss_genes
    ) ||
    any(
      !nzchar(
        consensus_skin_loss_genes
      )
    )
) {
  stop(
    "共识基因列表仍含NA或空字符串。"
  )
}

if (
  nrow(
    consensus_late_table
  ) !=
    length(
      consensus_late_genes
    ) ||
    nrow(
      consensus_skin_loss_table
    ) !=
    length(
      consensus_skin_loss_genes
    )
) {
  stop(
    "共识基因表存在重复基因；请检查冻结签名和合并过程。"
  )
}

consensus_late_table <- consensus_late_table[
  order(
    -consensus_late_table$minimum_oriented_core_effect
  ),
  ,
  drop = FALSE
]

consensus_skin_loss_table <- consensus_skin_loss_table[
  order(
    -consensus_skin_loss_table$minimum_oriented_core_effect
  ),
  ,
  drop = FALSE
]

safe_write_csv(
  consensus_late_table,
  file.path(
    REPORT_DIR,
    "08_consensus_persistent_late_remodeling_genes.csv"
  )
)

safe_write_csv(
  consensus_skin_loss_table,
  file.path(
    REPORT_DIR,
    "09_consensus_loss_of_skin_homeostasis_genes.csv"
  )
)

# ----------------------------- 阶段竞争假设 ------------------------------------
stage_competition_rows <- list()

for (
  stage in c(
    "Skin",
    "Wound1",
    "Wound7",
    "Wound30"
  )
) {
  x <- signature_gene_table[
    signature_gene_table$stage ==
      stage &
      !is.na(
        signature_gene_table$equal_weight_mean_difference
      ),
    ,
    drop = FALSE
  ]

  expected_match <- if (
    stage ==
      "Skin"
  ) {
    x$mean_difference_GSE163973 <
      0 &
      x$mean_difference_GSE181316 <
      0
  } else {
    x$mean_difference_GSE163973 >
      0 &
      x$mean_difference_GSE181316 >
      0
  }

  oriented_integrated <- if (
    stage ==
      "Skin"
  ) {
    -x$equal_weight_mean_difference
  } else {
    x$equal_weight_mean_difference
  }

  stage_competition_rows[[length(
    stage_competition_rows
  ) + 1L]] <- data.frame(
    stage =
      stage,
    n_frozen_genes =
      length(
        signatures[[stage]]
      ),
    n_measured_in_both_core_cohorts =
      nrow(
        x
      ),
    n_expected_direction_in_both =
      sum(
        expected_match,
        na.rm = TRUE
      ),
    proportion_expected_direction_in_both =
      mean(
        expected_match,
        na.rm = TRUE
      ),
    median_oriented_integrated_effect =
      stats::median(
        oriented_integrated,
        na.rm = TRUE
      ),
    mean_oriented_integrated_effect =
      mean(
        oriented_integrated,
        na.rm = TRUE
      ),
    n_stratified_FDR_below_0_10 =
      sum(
        x$BH_FDR_stratified <
          0.10,
        na.rm = TRUE
      ),
    stringsAsFactors = FALSE
  )
}

stage_competition_summary <- do.call(
  rbind,
  stage_competition_rows
)

safe_write_csv(
  stage_competition_summary,
  file.path(
    REPORT_DIR,
    "10_stage_signature_competition_summary.csv"
  )
)

# ----------------------------- 预定义状态面板 ----------------------------------
# 以下面板用于探索性解释，不作为诊断分类器。
# 文献锚点：
# - mesenchymal / papillary：Deng et al., Nat Commun 2021；
# - mechanoresponsive：Cheng et al., Commun Biol 2024；
# - pro-inflammatory、reticular、myofibroblast为透明的经典功能标志物面板。
known_state_panels <- list(
  Keloid_mesenchymal_literature = c(
    "COL11A1",
    "POSTN",
    "COMP"
  ),
  Secretory_papillary_literature = c(
    "COL13A1",
    "COL18A1",
    "COL23A1"
  ),
  Mechanoresponsive_keloid_literature = c(
    "YAP1",
    "WWTR1",
    "PIEZO1",
    "RHOA",
    "ROCK2",
    "FN1",
    "ITGA1",
    "ITGB1",
    "CD44"
  ),
  Pro_inflammatory_fibroblast_canonical = c(
    "IL6",
    "CXCL1",
    "CXCL2",
    "CXCL3",
    "CXCL8",
    "CCL2",
    "PTGS2",
    "ICAM1",
    "NFKBIA",
    "TNFAIP3"
  ),
  Secretory_reticular_canonical = c(
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "DPT",
    "COL14A1",
    "C7",
    "FBLN1"
  ),
  Contractile_myofibroblast_canonical = c(
    "ACTA2",
    "TAGLN",
    "MYL9",
    "TPM2",
    "CNN1",
    "CALD1",
    "COL1A1",
    "COL3A1"
  )
)

panel_definition <- do.call(
  rbind,
  lapply(
    names(
      known_state_panels
    ),
    function(panel) {
      data.frame(
        panel =
          panel,
        gene =
          unique(
            toupper(
              known_state_panels[[panel]]
            )
          ),
        panel_role =
          "exploratory_literature_or_canonical_state_interpretation",
        stringsAsFactors = FALSE
      )
    }
  )
)

safe_write_csv(
  panel_definition,
  file.path(
    REPORT_DIR,
    "11_predefined_fibroblast_state_panels.csv"
  )
)

hypergeometric_overlap <- function(
  query,
  target,
  universe
) {
  query <- intersect(
    unique(
      query
    ),
    universe
  )

  target <- intersect(
    unique(
      target
    ),
    universe
  )

  overlap <- intersect(
    query,
    target
  )

  m <- length(
    target
  )

  n <- length(
    setdiff(
      universe,
      target
    )
  )

  k <- length(
    query
  )

  q <- length(
    overlap
  )

  p <- if (
    k >
      0 &&
      m >
      0
  ) {
    stats::phyper(
      q -
        1,
      m,
      n,
      k,
      lower.tail = FALSE
    )
  } else {
    NA_real_
  }

  union_size <- length(
    union(
      query,
      target
    )
  )

  data.frame(
    query_size =
      length(
        query
      ),
    target_size =
      length(
        target
      ),
    overlap_size =
      q,
    jaccard =
      if (
        union_size >
          0
      ) {
        q /
          union_size
      } else {
        NA_real_
      },
    hypergeometric_p =
      p,
    overlap_genes =
      paste(
        sort(
          overlap
        ),
        collapse = ";"
      ),
    stringsAsFactors = FALSE
  )
}

panel_overlap_rows <- list()

for (
  panel in names(
    known_state_panels
  )
) {
  for (
    module_name in c(
      "persistent_late_remodeling",
      "loss_of_skin_homeostasis"
    )
  ) {
    module_genes <- if (
      module_name ==
        "persistent_late_remodeling"
    ) {
      consensus_late_genes
    } else {
      consensus_skin_loss_genes
    }

    overlap <- hypergeometric_overlap(
      query =
        known_state_panels[[panel]],
      target =
        module_genes,
      universe =
        common_gene_universe
    )

    overlap$panel <-
      panel

    overlap$target_module <-
      module_name

    panel_overlap_rows[[length(
      panel_overlap_rows
    ) + 1L]] <-
      overlap
  }
}

panel_overlap <- do.call(
  rbind,
  panel_overlap_rows
)

panel_overlap$BH_FDR <-
  stats::p.adjust(
    panel_overlap$hypergeometric_p,
    method = "BH"
  )

safe_write_csv(
  panel_overlap,
  file.path(
    REPORT_DIR,
    "12_fibroblast_state_panel_overlap.csv"
  )
)

# 患者级模块评分。
panel_score_rows <- list()
panel_sample_score_rows <- list()

for (
  panel in names(
    known_state_panels
  )
) {
  score_5 <- score_module(
    stage5_core$log_cpm,
    known_state_panels[[panel]]
  )

  score_6 <- score_module(
    stage6_core$log_cpm,
    known_state_panels[[panel]]
  )

  coverage_5 <- length(
    intersect(
      known_state_panels[[panel]],
      rownames(
        stage5_core$log_cpm
      )
    )
  )

  coverage_6 <- length(
    intersect(
      known_state_panels[[panel]],
      rownames(
        stage6_core$log_cpm
      )
    )
  )

  if (
    all(
      is.finite(
        score_5
      )
    ) &&
      all(
        is.finite(
          score_6
        )
      ) &&
      coverage_5 >=
      2L &&
      coverage_6 >=
      2L
  ) {
    integrated <- integrate_module_scores(
      score_1 =
        score_5,
      group_1 =
        stage5_core$sample_info$group,
      score_2 =
        score_6,
      group_2 =
        stage6_core$sample_info$group
    )

    panel_score_rows[[length(
      panel_score_rows
    ) + 1L]] <- data.frame(
      panel =
        panel,
      genes_defined =
        length(
          known_state_panels[[panel]]
        ),
      genes_present_GSE163973 =
        coverage_5,
      genes_present_GSE181316 =
        coverage_6,
      effect_GSE163973 =
        integrated$effect_GSE163973,
      cliffs_delta_GSE163973 =
        integrated$cliff_GSE163973,
      effect_GSE181316 =
        integrated$effect_GSE181316,
      cliffs_delta_GSE181316 =
        integrated$cliff_GSE181316,
      integrated_effect =
        integrated$integrated_effect,
      exact_stratified_two_sided_p =
        integrated$exact_stratified_two_sided_p,
      exact_stratified_positive_p =
        integrated$exact_stratified_positive_p,
      stringsAsFactors = FALSE
    )

    panel_sample_score_rows[[length(
      panel_sample_score_rows
    ) + 1L]] <- rbind(
      data.frame(
        panel =
          panel,
        cohort =
          "GSE163973",
        sample_id =
          stage5_core$sample_info$sample_id,
        group =
          stage5_core$sample_info$group,
        module_score =
          score_5,
        stringsAsFactors = FALSE
      ),
      data.frame(
        panel =
          panel,
        cohort =
          "GSE181316",
        sample_id =
          stage6_core$sample_info$sample_id,
        group =
          stage6_core$sample_info$group,
        module_score =
          score_6,
        stringsAsFactors = FALSE
      )
    )
  } else {
    warn_msg(
      "状态面板覆盖不足，跳过模块评分：",
      panel
    )
  }
}

panel_score_results <- if (length(
  panel_score_rows
)) {
  do.call(
    rbind,
    panel_score_rows
  )
} else {
  data.frame()
}

panel_sample_scores <- if (length(
  panel_sample_score_rows
)) {
  do.call(
    rbind,
    panel_sample_score_rows
  )
} else {
  data.frame()
}

if (nrow(
  panel_score_results
)) {
  panel_score_results$BH_FDR <-
    stats::p.adjust(
      panel_score_results$exact_stratified_two_sided_p,
      method = "BH"
    )
}

safe_write_csv(
  panel_score_results,
  file.path(
    REPORT_DIR,
    "13_fibroblast_state_panel_patient_effects.csv"
  )
)

safe_write_csv(
  panel_sample_scores,
  file.path(
    REPORT_DIR,
    "14_fibroblast_state_panel_sample_scores.csv"
  )
)

# ----------------------------- 固定阶段签名GSEA --------------------------------
gene_score <- all_gene_effects$integrated_rank_score
names(
  gene_score
) <- all_gene_effects$gene

stage_term2gene <- do.call(
  rbind,
  lapply(
    names(
      signatures
    ),
    function(stage) {
      data.frame(
        term =
          paste0(
            "FROZEN_",
            toupper(
              stage
            ),
            "_FIBROBLAST_SIGNATURE"
          ),
        gene =
          signatures[[stage]],
        stringsAsFactors = FALSE
      )
    }
  )
)

stage_gsea <- run_gsea(
  gene_score =
    gene_score,
  term2gene =
    stage_term2gene,
  analysis_name =
    "Frozen_wound_stage_signature_GSEA"
)

safe_write_csv(
  stage_gsea,
  file.path(
    REPORT_DIR,
    "15_frozen_stage_signature_GSEA.csv"
  )
)

# ----------------------------- MSigDB通路分析 ----------------------------------
pathway_availability <- data.frame(
  component = c(
    "msigdbr",
    "clusterProfiler",
    "pathway_analysis_ready"
  ),
  available = c(
    requireNamespace(
      "msigdbr",
      quietly = TRUE
    ),
    requireNamespace(
      "clusterProfiler",
      quietly = TRUE
    ),
    pathway_packages_available
  ),
  stringsAsFactors = FALSE
)

hallmark_term2gene <- data.frame()
gobp_term2gene <- data.frame()
reactome_term2gene <- data.frame()

hallmark_gsea <- data.frame()
gobp_gsea <- data.frame()
reactome_gsea <- data.frame()

late_hallmark_ora <- data.frame()
late_gobp_ora <- data.frame()
late_reactome_ora <- data.frame()

skin_hallmark_ora <- data.frame()
skin_gobp_ora <- data.frame()
skin_reactome_ora <- data.frame()

if (pathway_packages_available) {
  log_msg(
    "读取MSigDB Hallmark、GO BP与Reactome基因集"
  )

  hallmark_term2gene <- safe_msigdb_collection(
    collection =
      "H",
    subcollection =
      NULL,
    label =
      "Hallmark"
  )

  gobp_term2gene <- safe_msigdb_collection(
    collection =
      "C5",
    subcollection =
      "GO:BP",
    label =
      "GO_BP"
  )

  reactome_term2gene <- safe_msigdb_collection(
    collection =
      "C2",
    subcollection =
      "CP:REACTOME",
    label =
      "Reactome"
  )

  hallmark_gsea <- run_gsea(
    gene_score,
    hallmark_term2gene,
    "Hallmark_GSEA"
  )

  gobp_gsea <- run_gsea(
    gene_score,
    gobp_term2gene,
    "GO_BP_GSEA"
  )

  reactome_gsea <- run_gsea(
    gene_score,
    reactome_term2gene,
    "Reactome_GSEA"
  )

  late_hallmark_ora <- run_ora(
    consensus_late_genes,
    common_gene_universe,
    hallmark_term2gene,
    "Persistent_late_Hallmark_ORA"
  )

  late_gobp_ora <- run_ora(
    consensus_late_genes,
    common_gene_universe,
    gobp_term2gene,
    "Persistent_late_GO_BP_ORA"
  )

  late_reactome_ora <- run_ora(
    consensus_late_genes,
    common_gene_universe,
    reactome_term2gene,
    "Persistent_late_Reactome_ORA"
  )

  skin_hallmark_ora <- run_ora(
    consensus_skin_loss_genes,
    common_gene_universe,
    hallmark_term2gene,
    "Skin_homeostasis_loss_Hallmark_ORA"
  )

  skin_gobp_ora <- run_ora(
    consensus_skin_loss_genes,
    common_gene_universe,
    gobp_term2gene,
    "Skin_homeostasis_loss_GO_BP_ORA"
  )

  skin_reactome_ora <- run_ora(
    consensus_skin_loss_genes,
    common_gene_universe,
    reactome_term2gene,
    "Skin_homeostasis_loss_Reactome_ORA"
  )
} else {
  warn_msg(
    "通路包不可用；基因一致性与状态面板分析继续完成，但GSEA/ORA为空。"
  )
}

pathway_availability$n_hallmark_gene_sets <-
  c(
    NA,
    NA,
    length(
      unique(
        hallmark_term2gene$term
      )
    )
  )

pathway_availability$n_GOBP_gene_sets <-
  c(
    NA,
    NA,
    length(
      unique(
        gobp_term2gene$term
      )
    )
  )

pathway_availability$n_Reactome_gene_sets <-
  c(
    NA,
    NA,
    length(
      unique(
        reactome_term2gene$term
      )
    )
  )

safe_write_csv(
  pathway_availability,
  file.path(
    REPORT_DIR,
    "16_pathway_resource_availability.csv"
  )
)

safe_write_csv(
  hallmark_gsea,
  file.path(
    REPORT_DIR,
    "17_Hallmark_GSEA.csv"
  )
)

safe_write_csv(
  gobp_gsea,
  file.path(
    REPORT_DIR,
    "18_GO_BP_GSEA.csv"
  )
)

safe_write_csv(
  reactome_gsea,
  file.path(
    REPORT_DIR,
    "19_Reactome_GSEA.csv"
  )
)

safe_write_csv(
  late_hallmark_ora,
  file.path(
    REPORT_DIR,
    "20_persistent_late_Hallmark_ORA.csv"
  )
)

safe_write_csv(
  late_gobp_ora,
  file.path(
    REPORT_DIR,
    "21_persistent_late_GO_BP_ORA.csv"
  )
)

safe_write_csv(
  late_reactome_ora,
  file.path(
    REPORT_DIR,
    "22_persistent_late_Reactome_ORA.csv"
  )
)

safe_write_csv(
  skin_hallmark_ora,
  file.path(
    REPORT_DIR,
    "23_skin_homeostasis_loss_Hallmark_ORA.csv"
  )
)

safe_write_csv(
  skin_gobp_ora,
  file.path(
    REPORT_DIR,
    "24_skin_homeostasis_loss_GO_BP_ORA.csv"
  )
)

safe_write_csv(
  skin_reactome_ora,
  file.path(
    REPORT_DIR,
    "25_skin_homeostasis_loss_Reactome_ORA.csv"
  )
)

# ----------------------------- 结果闸门 ----------------------------------------
get_stage_gsea_row <- function(pattern) {
  if (!nrow(
    stage_gsea
  )) {
    return(
      data.frame()
    )
  }

  stage_gsea[
    grepl(
      pattern,
      stage_gsea$ID,
      fixed = TRUE
    ),
    ,
    drop = FALSE
  ]
}

w30_gsea <- get_stage_gsea_row(
  "FROZEN_WOUND30"
)

skin_gsea <- get_stage_gsea_row(
  "FROZEN_SKIN"
)

w1_gsea <- get_stage_gsea_row(
  "FROZEN_WOUND1"
)

w7_gsea <- get_stage_gsea_row(
  "FROZEN_WOUND7"
)

w30_concordance_rate <- stage_competition_summary$proportion_expected_direction_in_both[
  stage_competition_summary$stage ==
    "Wound30"
]

skin_concordance_rate <- stage_competition_summary$proportion_expected_direction_in_both[
  stage_competition_summary$stage ==
    "Skin"
]

early_mean_concordance <- mean(
  stage_competition_summary$proportion_expected_direction_in_both[
    stage_competition_summary$stage %in%
      c(
        "Wound1",
        "Wound7"
      )
  ]
)

gene_module_available <- (
  length(
    consensus_late_genes
  ) >=
    5L &&
    length(
      consensus_skin_loss_genes
    ) >=
    5L
)

stage_direction_consistent <- (
  length(
    w30_concordance_rate
  ) ==
    1L &&
    length(
      skin_concordance_rate
    ) ==
    1L &&
    w30_concordance_rate >=
    0.40 &&
    skin_concordance_rate >=
    0.40
)

stage_gsea_consistent <- (
  nrow(
    w30_gsea
  ) ==
    1L &&
    nrow(
      skin_gsea
    ) ==
    1L &&
    w30_gsea$NES >
    0 &&
    skin_gsea$NES <
    0
)

late_more_coherent_than_early <- (
  w30_concordance_rate >=
    early_mean_concordance
)

pathway_analysis_completed <- (
  pathway_packages_available &&
    (
      nrow(
        hallmark_gsea
      ) >
        0L ||
        nrow(
          gobp_gsea
        ) >
        0L ||
        nrow(
          reactome_gsea
        ) >
        0L
    )
)

panel_analysis_completed <- (
  nrow(
    panel_score_results
  ) >=
    4L
)

if (
  gene_module_available &&
    stage_direction_consistent &&
    stage_gsea_consistent &&
    late_more_coherent_than_early &&
    pathway_analysis_completed &&
    panel_analysis_completed
) {
  project_gate <-
    "ANALYSIS_COMPLETE_GENE_PROGRAM_INTERPRETATION_STRONG"

  interpretation <- paste(
    "两个核心队列在基因层面共同支持冻结Wound30程序增强和Skin稳态程序丧失；",
    "固定阶段签名GSEA方向一致，且通路与既有成纤维细胞状态分析已完成。",
    "该结果可用于构建手稿的生物学解释层。"
  )
} else if (
  gene_module_available &&
    stage_direction_consistent &&
    stage_gsea_consistent
) {
  project_gate <-
    "ANALYSIS_COMPLETE_GENE_PROGRAM_INTERPRETATION"

  interpretation <- paste(
    "核心晚期重塑与Skin稳态丧失模块已在基因层面建立，",
    "但通路资源、早期阶段比较或状态面板中的至少一项未达到全部内部标准。"
  )
} else {
  project_gate <-
    "GENE_LEVEL_INTERPRETATION_MIXED"

  interpretation <- paste(
    "患者级状态评分仍由前面核心分析支持，",
    "但冻结阶段签名在基因层面的方向一致性不足；",
    "手稿应以状态投射为主，弱化具体共识基因和通路机制。"
  )
}

decision_criteria <- data.frame(
  criterion = c(
    "gene_module_available",
    "stage_direction_consistent",
    "stage_gsea_consistent",
    "late_more_coherent_than_early",
    "pathway_analysis_completed",
    "panel_analysis_completed"
  ),
  passed = c(
    gene_module_available,
    stage_direction_consistent,
    stage_gsea_consistent,
    late_more_coherent_than_early,
    pathway_analysis_completed,
    panel_analysis_completed
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  decision_criteria,
  file.path(
    REPORT_DIR,
    "26_STAGE12_decision_criteria.csv"
  )
)

# ----------------------------- 图形 --------------------------------------------
# 1. 两核心队列基因效应一致性散点图
pdf(
  file.path(
    FIG_DIR,
    "Figure_core_gene_effect_concordance.pdf"
  ),
  width = 8,
  height = 7
)

plot(
  all_gene_effects$mean_difference_GSE163973,
  all_gene_effects$mean_difference_GSE181316,
  pch = 1,
  cex = 0.35,
  xlab =
    "GSE163973 keloid-normal scar log2-CPM difference",
  ylab =
    "GSE181316 keloid-normal scar log2-CPM difference",
  main =
    "Cross-cohort fibroblast gene-effect concordance"
)

abline(
  h = 0,
  v = 0,
  lty = 2
)

late_idx <- match(
  consensus_late_genes,
  all_gene_effects$gene
)

late_idx <- late_idx[
  !is.na(
    late_idx
  )
]

skin_idx <- match(
  consensus_skin_loss_genes,
  all_gene_effects$gene
)

skin_idx <- skin_idx[
  !is.na(
    skin_idx
  )
]

points(
  all_gene_effects$mean_difference_GSE163973[
    late_idx
  ],
  all_gene_effects$mean_difference_GSE181316[
    late_idx
  ],
  pch = 16,
  cex = 0.7
)

points(
  all_gene_effects$mean_difference_GSE163973[
    skin_idx
  ],
  all_gene_effects$mean_difference_GSE181316[
    skin_idx
  ],
  pch = 17,
  cex = 0.7
)

legend(
  "topleft",
  legend = c(
    "All common genes",
    "Consensus Wound30-up",
    "Consensus Skin-down"
  ),
  pch = c(
    1,
    16,
    17
  ),
  bty = "n"
)

dev.off()

# 2. 固定阶段基因方向一致性
pdf(
  file.path(
    FIG_DIR,
    "Figure_frozen_stage_gene_concordance.pdf"
  ),
  width = 9,
  height = 6
)

barplot(
  stage_competition_summary$proportion_expected_direction_in_both,
  names.arg =
    stage_competition_summary$stage,
  ylim = c(
    0,
    1
  ),
  ylab =
    "Proportion concordant in expected direction",
  main =
    "Frozen wound-stage signature concordance"
)

abline(
  h = 0.5,
  lty = 2
)

dev.off()

# 3. 共识基因效应
top_late <- head(
  consensus_late_table,
  20L
)

top_skin <- head(
  consensus_skin_loss_table,
  20L
)

if (
  nrow(
    top_late
  ) >
    0L ||
    nrow(
      top_skin
    ) >
    0L
) {
  plot_genes <- rbind(
    data.frame(
      module =
        "Wound30-up",
      gene =
        top_late$gene,
      effect =
        top_late$equal_weight_mean_difference,
      stringsAsFactors = FALSE
    ),
    data.frame(
      module =
        "Skin-down",
      gene =
        top_skin$gene,
      effect =
        top_skin$equal_weight_mean_difference,
      stringsAsFactors = FALSE
    )
  )

  pdf(
    file.path(
      FIG_DIR,
      "Figure_consensus_gene_effects.pdf"
    ),
    width = 10,
    height = max(
      7,
      0.28 *
        nrow(
          plot_genes
        ) +
        2
    )
  )

  par(
    mar = c(
      5,
      9,
      3,
      2
    )
  )

  barplot(
    rev(
      plot_genes$effect
    ),
    horiz = TRUE,
    names.arg =
      rev(
        paste(
          plot_genes$gene,
          plot_genes$module,
          sep = " | "
        )
      ),
    las = 1,
    xlab =
      "Equal-weight core-cohort log2-CPM difference",
    main =
      "Cross-cohort consensus wound-state genes"
  )

  abline(
    v = 0,
    lty = 2
  )

  dev.off()
}

# 4. 固定阶段GSEA
plot_top_enrichment(
  stage_gsea,
  file.path(
    FIG_DIR,
    "Figure_frozen_stage_signature_GSEA.pdf"
  ),
  "Frozen wound-stage signature GSEA",
  top_n = 4L
)

# 5. 通路图
plot_top_enrichment(
  hallmark_gsea,
  file.path(
    FIG_DIR,
    "Figure_Hallmark_GSEA_top.pdf"
  ),
  "Hallmark GSEA",
  top_n = 15L
)

plot_top_enrichment(
  gobp_gsea,
  file.path(
    FIG_DIR,
    "Figure_GO_BP_GSEA_top.pdf"
  ),
  "GO Biological Process GSEA",
  top_n = 20L
)

plot_top_enrichment(
  reactome_gsea,
  file.path(
    FIG_DIR,
    "Figure_Reactome_GSEA_top.pdf"
  ),
  "Reactome GSEA",
  top_n = 15L
)

plot_top_enrichment(
  late_gobp_ora,
  file.path(
    FIG_DIR,
    "Figure_persistent_late_GO_BP_ORA_top.pdf"
  ),
  "Persistent late-remodeling genes: GO BP ORA",
  top_n = 15L
)

plot_top_enrichment(
  skin_gobp_ora,
  file.path(
    FIG_DIR,
    "Figure_skin_homeostasis_loss_GO_BP_ORA_top.pdf"
  ),
  "Loss-of-homeostasis genes: GO BP ORA",
  top_n = 15L
)

# 6. 状态面板患者级效应
if (nrow(
  panel_score_results
)) {
  panel_order <- order(
    panel_score_results$integrated_effect
  )

  pdf(
    file.path(
      FIG_DIR,
      "Figure_fibroblast_state_panel_effects.pdf"
    ),
    width = 10,
    height = 7
  )

  par(
    mar = c(
      5,
      12,
      3,
      2
    )
  )

  barplot(
    panel_score_results$integrated_effect[
      panel_order
    ],
    horiz = TRUE,
    names.arg =
      panel_score_results$panel[
        panel_order
      ],
    las = 1,
    xlab =
      "Integrated keloid-normal scar module-score difference",
    main =
      "Predefined fibroblast-state panel effects"
  )

  abline(
    v = 0,
    lty = 2
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
    input_registry =
      input_registry,
    signatures =
      signatures,
    core_sample_structure =
      core_sample_structure,
    all_gene_effects =
      all_gene_effects,
    signature_gene_table =
      signature_gene_table,
    consensus_late_genes =
      consensus_late_genes,
    consensus_skin_loss_genes =
      consensus_skin_loss_genes,
    stage_competition_summary =
      stage_competition_summary,
    known_state_panels =
      known_state_panels,
    panel_overlap =
      panel_overlap,
    panel_score_results =
      panel_score_results,
    panel_sample_scores =
      panel_sample_scores,
    stage_gsea =
      stage_gsea,
    hallmark_gsea =
      hallmark_gsea,
    gobp_gsea =
      gobp_gsea,
    reactome_gsea =
      reactome_gsea,
    late_ora = list(
      hallmark =
        late_hallmark_ora,
      GO_BP =
        late_gobp_ora,
      Reactome =
        late_reactome_ora
    ),
    skin_loss_ora = list(
      hallmark =
        skin_hallmark_ora,
      GO_BP =
        skin_gobp_ora,
      Reactome =
        skin_reactome_ora
    ),
    decision_criteria =
      decision_criteria
  ),
  file = file.path(
    OBJECT_DIR,
    "Keloid_WoundAge_gene_program_interpretation.rds"
  ),
  compress = "gzip"
)

# ----------------------------- 最终结论 ----------------------------------------
get_gsea_value <- function(
  x,
  field
) {
  if (
    nrow(
      x
    ) ==
      1L &&
      field %in%
      names(
        x
      )
  ) {
    x[[field]][1]
  } else {
    NA
  }
}

decision_lines <- c(
  "第十二阶段：跨核心队列一致性基因、通路与成纤维细胞状态解释（V1.2最终修正版）",
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
  "总体解释：",
  interpretation,
  "",
  "共识模块：",
  paste0(
    "- Wound30稳定签名总数：",
    length(
      signatures$Wound30
    )
  ),
  paste0(
    "- 两个核心队列均上调的Wound30基因：",
    length(
      consensus_late_genes
    )
  ),
  paste0(
    "- Skin稳定签名总数：",
    length(
      signatures$Skin
    )
  ),
  paste0(
    "- 两个核心队列均下调的Skin基因：",
    length(
      consensus_skin_loss_genes
    )
  ),
  "",
  "阶段方向一致性：",
  paste0(
    "- Wound30预期方向一致比例：",
    round(
      w30_concordance_rate,
      3
    )
  ),
  paste0(
    "- Skin预期方向一致比例：",
    round(
      skin_concordance_rate,
      3
    )
  ),
  paste0(
    "- Wound1/Wound7平均一致比例：",
    round(
      early_mean_concordance,
      3
    )
  ),
  "",
  "固定阶段签名GSEA：",
  paste0(
    "- Wound30 NES：",
    round(
      get_gsea_value(
        w30_gsea,
        "NES"
      ),
      3
    ),
    "；FDR：",
    format(
      get_gsea_value(
        w30_gsea,
        "p.adjust"
      ),
      digits = 4
    )
  ),
  paste0(
    "- Skin NES：",
    round(
      get_gsea_value(
        skin_gsea,
        "NES"
      ),
      3
    ),
    "；FDR：",
    format(
      get_gsea_value(
        skin_gsea,
        "p.adjust"
      ),
      digits = 4
    )
  ),
  paste0(
    "- Wound1 NES：",
    round(
      get_gsea_value(
        w1_gsea,
        "NES"
      ),
      3
    ),
    "；FDR：",
    format(
      get_gsea_value(
        w1_gsea,
        "p.adjust"
      ),
      digits = 4
    )
  ),
  paste0(
    "- Wound7 NES：",
    round(
      get_gsea_value(
        w7_gsea,
        "NES"
      ),
      3
    ),
    "；FDR：",
    format(
      get_gsea_value(
        w7_gsea,
        "p.adjust"
      ),
      digits = 4
    )
  ),
  "- Wound7也呈正向富集，但方向一致比例和NES均低于Wound30；",
  "  因此应解释为中晚期程序存在重叠，而不是Wound30完全排他。",
  "",
  "通路资源：",
  paste0(
    "- Pathway analysis completed：",
    pathway_analysis_completed
  ),
  paste0(
    "- Hallmark GSEA terms：",
    nrow(
      hallmark_gsea
    )
  ),
  paste0(
    "- GO BP GSEA terms：",
    nrow(
      gobp_gsea
    )
  ),
  paste0(
    "- Reactome GSEA terms：",
    nrow(
      reactome_gsea
    )
  ),
  "",
  "解释边界：",
  "- 共识基因依据冻结签名成员资格和两个核心队列方向一致定义，不依据P值挑选。",
  "- 基因层面精确P值和FDR为探索性，不替代患者级主要终点。",
  "- GSE220300只作为第三队列方向标记，不决定核心共识基因是否纳入。",
  "- 通路富集用于解释状态程序，不能证明通路因果活化。",
  "- 文献/经典成纤维细胞面板是探索性比较，不是新诊断分型。",
  "- 不根据本阶段结果修改前面已经冻结的主要分析。",
  "- 本阶段完成后停止增加分析，进入投稿级作图和手稿写作。"
)

safe_write_lines(
  decision_lines,
  file.path(
    REPORT_DIR,
    "27_STAGE12_DECISION.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    REPORT_DIR,
    "28_SESSION_INFO.txt"
  )
)

safe_write_lines(
  c(
    paste0(
      "STAGE12_V1_2_COMPLETED=",
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
      "N_CONSENSUS_LATE_GENES=",
      length(
        consensus_late_genes
      )
    ),
    paste0(
      "N_CONSENSUS_SKIN_LOSS_GENES=",
      length(
        consensus_skin_loss_genes
      )
    ),
    paste0(
      "STAGE_DIRECTION_CONSISTENT=",
      stage_direction_consistent
    ),
    paste0(
      "STAGE_GSEA_CONSISTENT=",
      stage_gsea_consistent
    ),
    paste0(
      "PATHWAY_ANALYSIS_COMPLETED=",
      pathway_analysis_completed
    ),
    paste0(
      "PANEL_ANALYSIS_COMPLETED=",
      panel_analysis_completed
    )
  ),
  file.path(
    REPORT_DIR,
    "STAGE12_COMPLETED.txt"
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
  "第十二阶段完成"
)

cat(
  "\n============================================================\n"
)

cat(
  "第十二阶段运行完成。\n"
)

cat(
  "完整分析对象保存在：",
  file.path(
    OBJECT_DIR,
    "Keloid_WoundAge_gene_program_interpretation.rds"
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
