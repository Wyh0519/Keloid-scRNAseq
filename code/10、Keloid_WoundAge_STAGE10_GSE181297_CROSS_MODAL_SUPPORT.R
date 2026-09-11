# =============================================================================
# 项目：瘢痕疙瘩“晚期重塑状态持续与伤口状态终止失败”
# 第十阶段：GSE181297单细胞—Visium跨模态固定状态支持分析
#
# 根目录：
#   C:/Users/33652/Desktop/LZZ/1、yssj
#
# 数据结构：
#   scRNA-seq:
#     Ke01  keloid
#     Ke02  keloid
#     NS02  normal scar
#
#   Visium:
#     Pt1   keloid
#     Pt2   keloid
#     NSV1  adjacent normal
#     NSV2  adjacent normal
#
# 重要限制：
#   GSE181297_RAW.tar只提供Visium表达矩阵、barcode、features和低分辨率图像，
#   未提供tissue_positions或scalefactors。因此本阶段不能：
#   - 将spot叠加到组织图像；
#   - 计算Moran's I或局部空间自相关；
#   - 命名伤口边缘、深部真皮或病理空间域；
#   - 声称精确空间定位。
#
# 本阶段可以：
#   1. 用冻结成纤维细胞规则提取scRNA fibroblast-like细胞；
#   2. 用冻结规则识别Visium fibroblast-enriched spots；
#   3. 使用固定正常伤口阶段签名进行切片/样本级pseudobulk映射；
#   4. 检验同一GEO队列内scRNA和Visium是否呈一致方向；
#   5. 以spot分布作描述性支持，但不把spots当作生物学重复。
#
# 统计边界：
#   - scRNA仅2例keloid和1例normal scar；
#   - Visium仅2例keloid和2例adjacent normal；
#   - 不进行常规显著性推断；
#   - 只报告样本级方向、效应方向和跨模态一致性；
#   - 不根据GSE181297结果修改前面冻结的主要终点。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANG = "en")

ROOT_DIR <- "C:/Users/33652/Desktop/LZZ/1、yssj"

RAW_TAR <- file.path(
  ROOT_DIR,
  "GSE181297_RAW.tar"
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

FROZEN_MARKER_RDS <- file.path(
  ROOT_DIR,
  "06A_STAGE6A_GSE181316_CELL_CALLING",
  "objects",
  "GSE181316_fibroblast_like_pseudobulk_list.rds"
)

STAGE_DIR <- file.path(
  ROOT_DIR,
  "10_STAGE10_GSE181297_CROSS_MODAL_SUPPORT"
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
  "第十阶段_GSE181297单细胞与Visium跨模态支持检查包.zip"
)

PACKAGE_TARGZ <- file.path(
  ROOT_DIR,
  "第十阶段_GSE181297单细胞与Visium跨模态支持检查包.tar.gz"
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
  "00_stage10_log.txt"
)

WARNING_FILE <- file.path(
  REPORT_DIR,
  "00_stage10_warnings.txt"
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
    "Ke01",
    "Ke02",
    "NS02",
    "Pt1",
    "Pt2",
    "NSV1",
    "NSV2"
  ),
  modality = c(
    "scRNA",
    "scRNA",
    "scRNA",
    "Visium",
    "Visium",
    "Visium",
    "Visium"
  ),
  group = c(
    "keloid",
    "keloid",
    "normal_scar",
    "keloid",
    "keloid",
    "adjacent_normal",
    "adjacent_normal"
  ),
  location = c(
    "ear",
    "back",
    "back",
    "ear",
    "abdomen",
    "leg",
    "cheek"
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

# ----------------------------- 基础读取函数 ------------------------------------
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
      ignore.case = TRUE,
      perl = TRUE
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
      ignore.case = TRUE,
      perl = TRUE
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
      ignore.case = TRUE,
      perl = TRUE
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

# ----------------------------- 基因与标志物函数 --------------------------------
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
      "没有有效基因可用于pseudobulk合并。"
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

extract_frozen_marker_sets <- function(
  marker_object
) {
  if (!"marker_definition" %in%
      names(
        marker_object
      )) {
    stop(
      "第六阶段A对象缺少marker_definition。"
    )
  }

  definition <- marker_object$marker_definition

  required_sets <- c(
    "fibro_core",
    "fibro_support",
    "immune",
    "epithelial",
    "endothelial",
    "neural"
  )

  result <- setNames(
    lapply(
      required_sets,
      function(marker_set) {
        unique(
          toupper(
            trimws(
              as.character(
                definition$gene[
                  definition$marker_set ==
                    marker_set
                ]
              )
            )
          )
        )
      }
    ),
    required_sets
  )

  if (any(
    vapply(
      result,
      length,
      integer(1)
    ) ==
      0L
  )) {
    stop(
      "至少一个冻结标志物集合为空。"
    )
  }

  result
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

  output <- matrix(
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
      output[i, ] <- Matrix::colSums(
        counts[
          idx,
          ,
          drop = FALSE
        ]
      )
    }
  }

  output
}

marker_metrics <- function(
  marker_counts,
  library_size
) {
  detected <- colSums(
    marker_counts >
      0
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
    ) *
      10000
  )

  list(
    detected = as.numeric(
      detected
    ),
    score = as.numeric(
      colMeans(
        normalized
      )
    )
  )
}

calculate_marker_results <- function(
  counts,
  gene_symbol,
  marker_sets,
  library_size
) {
  result <- list()

  for (
    marker_name in names(
      marker_sets
    )
  ) {
    marker_counts <- marker_count_matrix(
      counts,
      gene_symbol,
      marker_sets[[marker_name]]
    )

    result[[marker_name]] <- marker_metrics(
      marker_counts,
      library_size
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

# ----------------------------- 正常伤口映射函数 --------------------------------
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

  common_genes <- Reduce(
    intersect,
    lapply(
      vector_list,
      names
    )
  )

  if (
    length(
      common_genes
    ) <
      10000L
  ) {
    warn_msg(
      "样本间共同基因少于10,000：",
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

      if (
        length(
          genes
        ) <
          15L
      ) {
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

analyse_pseudobulk_set <- function(
  set_name,
  vector_list,
  required_samples,
  design,
  reference_counts,
  reference_design,
  signatures
) {
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
      "与正常参考共同背景基因少于5000：",
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

  reference_z <- standardized$reference_z
  new_z <- standardized$new_z

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
    design[
      ,
      c(
        "sample_id",
        "modality",
        "group",
        "location"
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
    coverage =
      coverage,
    calibration = list(
      common_genes =
        common_genes,
      signatures =
        signatures,
      center =
        standardized$center,
      scale =
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

cliffs_delta_descriptive <- function(
  x,
  y
) {
  pair_difference <- outer(
    x,
    y,
    "-"
  )

  mean(
    pair_difference >
      0
  ) -
    mean(
      pair_difference <
        0
    )
}

# ----------------------------- 输入检查 ----------------------------------------
log_msg(
  "第十阶段开始"
)

required_files <- c(
  RAW_TAR,
  REFERENCE_PB_RDS,
  REFERENCE_MODEL_RDS,
  FROZEN_MARKER_RDS
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

reference_pb_object <- readRDS(
  REFERENCE_PB_RDS
)

reference_model <- readRDS(
  REFERENCE_MODEL_RDS
)

marker_object <- readRDS(
  FROZEN_MARKER_RDS
)

marker_sets <- extract_frozen_marker_sets(
  marker_object
)

safe_write_csv(
  marker_object$marker_definition,
  file.path(
    REPORT_DIR,
    "02_reused_frozen_marker_definition.csv"
  )
)

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
    "03_fixed_signature_summary.csv"
  )
)

reference_design_all <- get_reference_design(
  reference_pb_object$counts,
  reference_pb_object$sample_info
)

reference_counts_all <- canonicalize_reference_counts(
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

# ----------------------------- 解包与结构审计 ----------------------------------
outer_dir <- file.path(
  WORK_DIR,
  "GSE181297_RAW"
)

dir.create(
  outer_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

log_msg(
  "解包GSE181297_RAW.tar"
)

utils::untar(
  RAW_TAR,
  exdir = outer_dir
)

all_files <- list.files(
  outer_dir,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  include.dirs = FALSE,
  no.. = TRUE
)

structure_audit <- data.frame(
  n_total_files =
    length(
      all_files
    ),
  n_matrix_files =
    sum(
      grepl(
        "_matrix\\.mtx\\.gz$",
        basename(
          all_files
        ),
        ignore.case = TRUE
      )
    ),
  n_barcode_files =
    sum(
      grepl(
        "_barcodes\\.tsv\\.gz$",
        basename(
          all_files
        ),
        ignore.case = TRUE
      )
    ),
  n_feature_files =
    sum(
      grepl(
        "_features\\.tsv\\.gz$",
        basename(
          all_files
        ),
        ignore.case = TRUE
      )
    ),
  n_lowres_images =
    sum(
      grepl(
        "tissue_lowres_image\\.png\\.gz$",
        basename(
          all_files
        ),
        ignore.case = TRUE
      )
    ),
  n_tissue_position_files =
    sum(
      grepl(
        "tissue_positions",
        basename(
          all_files
        ),
        ignore.case = TRUE
      )
    ),
  n_scalefactor_files =
    sum(
      grepl(
        "scalefactors",
        basename(
          all_files
        ),
        ignore.case = TRUE
      )
    ),
  spatial_coordinates_available =
    any(
      grepl(
        "tissue_positions",
        basename(
          all_files
        ),
        ignore.case = TRUE
      )
    ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  structure_audit,
  file.path(
    REPORT_DIR,
    "04_input_structure_and_coordinate_audit.csv"
  )
)

if (
  structure_audit$spatial_coordinates_available
) {
  warn_msg(
    "检测到tissue_positions文件；本代码仍按预设仅进行表达层面的跨模态支持分析。"
  )
} else {
  log_msg(
    "确认GSE181297无tissue_positions/scalefactors；不执行空间坐标分析。"
  )
}

# ----------------------------- 逐样本处理 --------------------------------------
scrna_high_list <- list()
scrna_broad_list <- list()
visium_high_list <- list()
visium_broad_list <- list()

scrna_qc_rows <- list()
visium_qc_rows <- list()
feature_rows <- list()
spot_distribution_rows <- list()

gene_reference_id <- NULL
gene_reference_symbol <- NULL

pdf(
  file.path(
    FIG_DIR,
    "Figure_GSE181297_QC_and_fibroblast_audit.pdf"
  ),
  width = 10,
  height = 8,
  onefile = TRUE
)

for (
  i in seq_len(
    nrow(
      sample_design
    )
  )
) {
  design_row <- sample_design[
    i,
    ,
    drop = FALSE
  ]

  sample_id <- design_row$sample_id
  modality <- design_row$modality

  log_msg(
    "处理样本：",
    sample_id,
    "|",
    modality
  )

  triplet <- locate_triplet(
    outer_dir,
    sample_id
  )

  features <- read_noheader_table(
    triplet$features,
    sep = "\t"
  )

  barcodes <- read_noheader_table(
    triplet$barcodes,
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

  barcode <- as.character(
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
        barcode
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

  colnames(
    counts
  ) <- barcode

  n_count <- Matrix::colSums(
    counts
  )

  n_feature <- Matrix::colSums(
    counts >
      0
  )

  gene_upper <- toupper(
    gene_symbol
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

  if (
    modality ==
      "scRNA"
  ) {
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
        n_count >=
        500 &
        n_feature >=
        200 &
        n_feature <=
        7500 &
        pct_mt <=
        20
    )
  } else {
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
        n_count >=
        200 &
        n_feature >=
        100 &
        pct_mt <=
        30
    )
  }

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
      "QC后细胞/spots少于100：",
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
      counts_qc >
        0
    )
  )

  pct_mt_qc <- as.numeric(
    pct_mt[
      qc_keep
    ]
  )

  marker_result <- calculate_marker_results(
    counts_qc,
    gene_reference_symbol,
    marker_sets,
    n_count_qc
  )

  max_exclusion_score <- pmax(
    marker_result$immune$score,
    marker_result$epithelial$score,
    marker_result$endothelial$score,
    marker_result$neural$score
  )

  fibro_identity_score <-
    marker_result$fibro_core$score +
    0.5 *
    marker_result$fibro_support$score -
    max_exclusion_score

  if (
    modality ==
      "scRNA"
  ) {
    broad_selected <- (
      marker_result$fibro_core$detected >=
        2 &
        marker_result$fibro_core$score >
        max_exclusion_score
    )

    high_selected <- (
      marker_result$fibro_core$detected >=
        3 &
        (
          marker_result$fibro_support$detected >=
            1 |
            marker_result$fibro_core$detected >=
            4
        ) &
        marker_result$fibro_core$score >
        max_exclusion_score +
        0.15 &
        marker_result$immune$detected <=
        1 &
        marker_result$epithelial$detected <=
        1 &
        marker_result$endothelial$detected <=
        1 &
        marker_result$neural$detected <=
        1
    )

    high_selected <- high_selected &
      broad_selected

    minimum_high <- 20L
    minimum_broad <- 50L
  } else {
    high_selected <- (
      marker_result$fibro_core$detected >=
        3 &
        marker_result$fibro_support$detected >=
        1 &
        marker_result$fibro_core$score >
        max_exclusion_score +
        0.15
    )

    broad_candidates <- (
      marker_result$fibro_core$detected >=
        2 &
        fibro_identity_score >
        0
    )

    candidate_values <-
      fibro_identity_score[
        broad_candidates &
          is.finite(
            fibro_identity_score
          )
      ]

    if (length(
      candidate_values
    )) {
      broad_cutoff <- as.numeric(
        stats::quantile(
          candidate_values,
          probs = 0.75,
          na.rm = TRUE,
          names = FALSE,
          type = 8
        )
      )
    } else {
      finite_values <-
        fibro_identity_score[
          is.finite(
            fibro_identity_score
          )
        ]

      if (!length(
        finite_values
      )) {
        stop(
          sample_id,
          "不存在有限fibro identity score。"
        )
      }

      broad_cutoff <- as.numeric(
        stats::quantile(
          finite_values,
          probs = 0.75,
          na.rm = TRUE,
          names = FALSE,
          type = 8
        )
      )

      warn_msg(
        sample_id,
        "没有spot满足broad candidate条件；使用全部有限score的75%分位数作回退。"
      )
    }

    broad_selected <- (
      broad_candidates &
        fibro_identity_score >=
        broad_cutoff
    )

    minimum_high <- 20L
    minimum_broad <- 20L
  }

  high_selected[
    is.na(
      high_selected
    )
  ] <- FALSE

  broad_selected[
    is.na(
      broad_selected
    )
  ] <- FALSE

  n_high <- sum(
    high_selected
  )

  n_broad <- sum(
    broad_selected
  )

  if (
    n_high <
      minimum_high
  ) {
    stop(
      sample_id,
      "高特异性集合少于最低分析量：",
      n_high
    )
  }

  if (
    n_broad <
      minimum_broad
  ) {
    stop(
      sample_id,
      "宽松集合少于最低分析量：",
      n_broad
    )
  }

  if (
    modality ==
      "scRNA"
  ) {
    if (
      n_high <
        50L
    ) {
      warn_msg(
        sample_id,
        "高特异性fibroblast-like细胞少于50：",
        n_high
      )
    }

    if (
      n_broad <
        100L
    ) {
      warn_msg(
        sample_id,
        "宽松fibroblast-like细胞少于100：",
        n_broad
      )
    }
  } else {
    if (
      n_high <
        50L
    ) {
      warn_msg(
        sample_id,
        "高特异性fibroblast-enriched spots少于50：",
        n_high
      )
    }

    if (
      n_broad <
        50L
    ) {
      warn_msg(
        sample_id,
        "top-quartile fibroblast-enriched spots少于50：",
        n_broad
      )
    }
  }

  high_raw <- Matrix::rowSums(
    counts_qc[
      ,
      high_selected,
      drop = FALSE
    ]
  )

  broad_raw <- Matrix::rowSums(
    counts_qc[
      ,
      broad_selected,
      drop = FALSE
    ]
  )

  high_pb <- collapse_vector_by_symbol(
    high_raw,
    gene_reference_symbol
  )

  broad_pb <- collapse_vector_by_symbol(
    broad_raw,
    gene_reference_symbol
  )

  if (
    modality ==
      "scRNA"
  ) {
    scrna_high_list[[sample_id]] <-
      high_pb

    scrna_broad_list[[sample_id]] <-
      broad_pb

    scrna_qc_rows[[length(
      scrna_qc_rows
    ) + 1L]] <- data.frame(
      sample_id =
        sample_id,
      group =
        design_row$group,
      location =
        design_row$location,
      n_input_cells =
        ncol(
          counts
        ),
      n_QC_cells =
        ncol(
          counts_qc
        ),
      QC_retention =
        ncol(
          counts_qc
        ) /
        ncol(
          counts
        ),
      n_high_fibroblast =
        n_high,
      n_broad_fibroblast =
        n_broad,
      high_fraction =
        n_high /
        ncol(
          counts_qc
        ),
      broad_fraction =
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
  } else {
    visium_high_list[[sample_id]] <-
      high_pb

    visium_broad_list[[sample_id]] <-
      broad_pb

    stage_module <- matrix(
      NA_real_,
      nrow = ncol(
        counts_qc
      ),
      ncol = 4,
      dimnames = list(
        colnames(
          counts_qc
        ),
        names(
          signatures
        )
      )
    )

    for (
      stage in names(
        signatures
      )
    ) {
      signature_counts <- marker_count_matrix(
        counts_qc,
        gene_reference_symbol,
        signatures[[stage]]
      )

      normalized <- log1p(
        sweep(
          signature_counts,
          2,
          pmax(
            n_count_qc,
            1
          ),
          "/"
        ) *
          10000
      )

      stage_module[, stage] <- colMeans(
        normalized
      )
    }

    spot_late_remodeling <-
      stage_module[, "Wound30"] -
      (
        stage_module[, "Wound1"] +
          stage_module[, "Wound7"]
      ) /
      2

    spot_distribution_rows[[length(
      spot_distribution_rows
    ) + 1L]] <- data.frame(
      sample_id =
        sample_id,
      group =
        design_row$group,
      n_QC_spots =
        ncol(
          counts_qc
        ),
      high_n =
        n_high,
      broad_n =
        n_broad,
      high_median_fibro_identity =
        stats::median(
          fibro_identity_score[
            high_selected
          ]
        ),
      broad_median_fibro_identity =
        stats::median(
          fibro_identity_score[
            broad_selected
          ]
        ),
      high_median_spot_late_remodeling =
        stats::median(
          spot_late_remodeling[
            high_selected
          ]
        ),
      broad_median_spot_late_remodeling =
        stats::median(
          spot_late_remodeling[
            broad_selected
          ]
        ),
      high_fraction_spot_late_positive =
        mean(
          spot_late_remodeling[
            high_selected
          ] >
            0
        ),
      broad_fraction_spot_late_positive =
        mean(
          spot_late_remodeling[
            broad_selected
          ] >
            0
        ),
      stringsAsFactors = FALSE
    )

    visium_qc_rows[[length(
      visium_qc_rows
    ) + 1L]] <- data.frame(
      sample_id =
        sample_id,
      group =
        design_row$group,
      location =
        design_row$location,
      n_input_spots =
        ncol(
          counts
        ),
      n_QC_spots =
        ncol(
          counts_qc
        ),
      QC_retention =
        ncol(
          counts_qc
        ) /
        ncol(
          counts
        ),
      n_high_fibroblast_enriched =
        n_high,
      n_broad_fibroblast_enriched =
        n_broad,
      high_fraction =
        n_high /
        ncol(
          counts_qc
        ),
      broad_fraction =
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
  }

  feature_rows[[length(
    feature_rows
  ) + 1L]] <- data.frame(
    sample_id =
      sample_id,
    modality =
      modality,
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
      n_count_qc +
        1
    ),
    breaks = 50,
    xlab =
      "log10 nCount + 1",
    main = paste0(
      sample_id,
      " QC counts"
    )
  )

  hist(
    n_feature_qc,
    breaks = 50,
    xlab =
      "Detected genes",
    main = paste0(
      sample_id,
      " QC features"
    )
  )

  hist(
    fibro_identity_score,
    breaks = 50,
    xlab =
      "Fibroblast identity score",
    main = paste0(
      sample_id,
      " fibroblast identity"
    )
  )

  plot(
    marker_result$fibro_core$score,
    max_exclusion_score,
    pch = ifelse(
      high_selected,
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
      " frozen classification"
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
    marker_result
  )

  gc(
    verbose = FALSE
  )
}

dev.off()

scrna_qc_summary <- do.call(
  rbind,
  scrna_qc_rows
)

visium_qc_summary <- do.call(
  rbind,
  visium_qc_rows
)

feature_summary <- do.call(
  rbind,
  feature_rows
)

spot_distribution_summary <- do.call(
  rbind,
  spot_distribution_rows
)

rownames(
  scrna_qc_summary
) <- NULL

rownames(
  visium_qc_summary
) <- NULL

rownames(
  feature_summary
) <- NULL

rownames(
  spot_distribution_summary
) <- NULL

safe_write_csv(
  scrna_qc_summary,
  file.path(
    REPORT_DIR,
    "05_scRNA_QC_fibroblast_summary.csv"
  )
)

safe_write_csv(
  visium_qc_summary,
  file.path(
    REPORT_DIR,
    "06_Visium_QC_fibroblast_spot_summary.csv"
  )
)

safe_write_csv(
  feature_summary,
  file.path(
    REPORT_DIR,
    "07_feature_consistency.csv"
  )
)

safe_write_csv(
  spot_distribution_summary,
  file.path(
    REPORT_DIR,
    "08_Visium_spot_distribution_summary.csv"
  )
)

# ----------------------------- 固定状态映射 ------------------------------------
scrna_design <- sample_design[
  sample_design$modality ==
    "scRNA",
  ,
  drop = FALSE
]

visium_design <- sample_design[
  sample_design$modality ==
    "Visium",
  ,
  drop = FALSE
]

scrna_high_result <- analyse_pseudobulk_set(
  set_name =
    "scRNA_high_specificity",
  vector_list =
    scrna_high_list,
  required_samples =
    scrna_design$sample_id,
  design =
    scrna_design,
  reference_counts =
    reference_counts,
  reference_design =
    reference_design,
  signatures =
    signatures
)

scrna_broad_result <- analyse_pseudobulk_set(
  set_name =
    "scRNA_broad_sensitivity",
  vector_list =
    scrna_broad_list,
  required_samples =
    scrna_design$sample_id,
  design =
    scrna_design,
  reference_counts =
    reference_counts,
  reference_design =
    reference_design,
  signatures =
    signatures
)

visium_high_result <- analyse_pseudobulk_set(
  set_name =
    "Visium_high_specificity_spots",
  vector_list =
    visium_high_list,
  required_samples =
    visium_design$sample_id,
  design =
    visium_design,
  reference_counts =
    reference_counts,
  reference_design =
    reference_design,
  signatures =
    signatures
)

visium_broad_result <- analyse_pseudobulk_set(
  set_name =
    "Visium_broad_top_quartile_spots",
  vector_list =
    visium_broad_list,
  required_samples =
    visium_design$sample_id,
  design =
    visium_design,
  reference_counts =
    reference_counts,
  reference_design =
    reference_design,
  signatures =
    signatures
)

all_scores <- rbind(
  scrna_high_result$mapping,
  scrna_broad_result$mapping,
  visium_high_result$mapping,
  visium_broad_result$mapping
)

coverage <- rbind(
  scrna_high_result$coverage,
  scrna_broad_result$coverage,
  visium_high_result$coverage,
  visium_broad_result$coverage
)

rownames(
  all_scores
) <- NULL

rownames(
  coverage
) <- NULL

safe_write_csv(
  coverage,
  file.path(
    REPORT_DIR,
    "09_signature_gene_coverage.csv"
  )
)

safe_write_csv(
  all_scores,
  file.path(
    REPORT_DIR,
    "10_GSE181297_cross_modal_wound_state_scores.csv"
  )
)

# ----------------------------- 跨模态方向判定 ----------------------------------
get_score_set <- function(
  analysis_set
) {
  x <- all_scores[
    all_scores$analysis_set ==
      analysis_set,
    ,
    drop = FALSE
  ]

  if (!nrow(
    x
  )) {
    stop(
      "找不到分析集合：",
      analysis_set
    )
  }

  x
}

direction_summary_rows <- list()

for (
  analysis_set in unique(
    all_scores$analysis_set
  )
) {
  x <- get_score_set(
    analysis_set
  )

  if (
    unique(
      x$modality
    ) ==
      "scRNA"
  ) {
    disease_value <- x$late_remodeling_state[
      x$group ==
        "keloid"
    ]

    control_value <- x$late_remodeling_state[
      x$group ==
        "normal_scar"
    ]

    disease_ordinal <- x$ordinal_wound_state_position[
      x$group ==
        "keloid"
    ]

    control_ordinal <- x$ordinal_wound_state_position[
      x$group ==
        "normal_scar"
    ]

    late_mean_difference <-
      mean(
        disease_value
      ) -
      mean(
        control_value
      )

    ordinal_mean_difference <-
      mean(
        disease_ordinal
      ) -
      mean(
        control_ordinal
      )

    complete_late_separation <-
      min(
        disease_value
      ) >
      max(
        control_value
      )

    complete_ordinal_separation <-
      min(
        disease_ordinal
      ) >
      max(
        control_ordinal
      )

    late_cliffs_delta <-
      cliffs_delta_descriptive(
        disease_value,
        control_value
      )

    ordinal_cliffs_delta <-
      cliffs_delta_descriptive(
        disease_ordinal,
        control_ordinal
      )

    n_disease <-
      length(
        disease_value
      )

    n_control <-
      length(
        control_value
      )
  } else {
    disease_value <- x$late_remodeling_state[
      x$group ==
        "keloid"
    ]

    control_value <- x$late_remodeling_state[
      x$group ==
        "adjacent_normal"
    ]

    disease_ordinal <- x$ordinal_wound_state_position[
      x$group ==
        "keloid"
    ]

    control_ordinal <- x$ordinal_wound_state_position[
      x$group ==
        "adjacent_normal"
    ]

    late_mean_difference <-
      mean(
        disease_value
      ) -
      mean(
        control_value
      )

    ordinal_mean_difference <-
      mean(
        disease_ordinal
      ) -
      mean(
        control_ordinal
      )

    complete_late_separation <-
      min(
        disease_value
      ) >
      max(
        control_value
      )

    complete_ordinal_separation <-
      min(
        disease_ordinal
      ) >
      max(
        control_ordinal
      )

    late_cliffs_delta <-
      cliffs_delta_descriptive(
        disease_value,
        control_value
      )

    ordinal_cliffs_delta <-
      cliffs_delta_descriptive(
        disease_ordinal,
        control_ordinal
      )

    n_disease <-
      length(
        disease_value
      )

    n_control <-
      length(
        control_value
      )
  }

  direction_summary_rows[[length(
    direction_summary_rows
  ) + 1L]] <- data.frame(
    analysis_set =
      analysis_set,
    modality =
      unique(
        x$modality
      ),
    n_keloid =
      n_disease,
    n_control =
      n_control,
    late_remodeling_mean_difference =
      late_mean_difference,
    late_remodeling_cliffs_delta =
      late_cliffs_delta,
    late_complete_separation =
      complete_late_separation,
    ordinal_mean_difference =
      ordinal_mean_difference,
    ordinal_cliffs_delta =
      ordinal_cliffs_delta,
    ordinal_complete_separation =
      complete_ordinal_separation,
    stringsAsFactors = FALSE
  )
}

direction_summary <- do.call(
  rbind,
  direction_summary_rows
)

rownames(
  direction_summary
) <- NULL

safe_write_csv(
  direction_summary,
  file.path(
    REPORT_DIR,
    "11_cross_modal_direction_summary.csv"
  )
)

get_direction_row <- function(
  analysis_set
) {
  x <- direction_summary[
    direction_summary$analysis_set ==
      analysis_set,
    ,
    drop = FALSE
  ]

  if (nrow(
    x
  ) !=
      1L) {
    stop(
      "无法唯一获得方向结果：",
      analysis_set
    )
  }

  x
}

scrna_high_direction <- get_direction_row(
  "scRNA_high_specificity"
)

scrna_broad_direction <- get_direction_row(
  "scRNA_broad_sensitivity"
)

visium_high_direction <- get_direction_row(
  "Visium_high_specificity_spots"
)

visium_broad_direction <- get_direction_row(
  "Visium_broad_top_quartile_spots"
)

scrna_high_support <- (
  scrna_high_direction$late_remodeling_mean_difference >
    0 &&
    scrna_high_direction$ordinal_mean_difference >
    0 &&
    scrna_high_direction$late_remodeling_cliffs_delta >
    0 &&
    scrna_high_direction$ordinal_cliffs_delta >
    0
)

scrna_broad_support <- (
  scrna_broad_direction$late_remodeling_mean_difference >
    0 &&
    scrna_broad_direction$ordinal_mean_difference >
    0
)

visium_high_support <- (
  visium_high_direction$late_remodeling_mean_difference >
    0 &&
    visium_high_direction$ordinal_mean_difference >
    0 &&
    visium_high_direction$late_remodeling_cliffs_delta >=
    0.5 &&
    visium_high_direction$ordinal_cliffs_delta >=
    0.5
)

visium_broad_support <- (
  visium_broad_direction$late_remodeling_mean_difference >
    0 &&
    visium_broad_direction$ordinal_mean_difference >
    0
)

cross_modal_consistency <- (
  scrna_high_support &&
    scrna_broad_support &&
    visium_high_support &&
    visium_broad_support
)

strong_complete_separation <- (
  scrna_high_direction$late_complete_separation &&
    visium_high_direction$late_complete_separation
)

if (
  cross_modal_consistency &&
    strong_complete_separation
) {
  project_gate <-
    "PASS_STRONG_CROSS_MODAL_SUPPORT"

  interpretation <- paste(
    "GSE181297的scRNA和Visium切片均显示瘢痕疙瘩",
    "late-remodeling state及ordinal wound-state position高于各自对照，",
    "且主要集合形成完整样本级分离。"
  )
} else if (
  cross_modal_consistency
) {
  project_gate <-
    "PASS_CROSS_MODAL_DIRECTION_SUPPORT"

  interpretation <- paste(
    "GSE181297的scRNA和Visium两种模态均呈一致的晚期重塑方向，",
    "但至少一种模态未形成完整样本级分离。"
  )
} else if (
  scrna_high_support ||
    visium_high_support
) {
  project_gate <-
    "PARTIAL_GSE181297_SUPPORT"

  interpretation <- paste(
    "GSE181297仅有一种模态或一种成纤维细胞定义支持晚期重塑方向，",
    "只能作为部分描述性证据。"
  )
} else {
  project_gate <-
    "NO_GSE181297_SUPPORT"

  interpretation <- paste(
    "GSE181297未提供一致的scRNA—Visium晚期重塑方向支持。"
  )
}

criteria <- data.frame(
  criterion = c(
    "scRNA_high_support",
    "scRNA_broad_support",
    "Visium_high_support",
    "Visium_broad_support",
    "cross_modal_consistency",
    "strong_complete_separation",
    "spatial_coordinates_available"
  ),
  passed = c(
    scrna_high_support,
    scrna_broad_support,
    visium_high_support,
    visium_broad_support,
    cross_modal_consistency,
    strong_complete_separation,
    structure_audit$spatial_coordinates_available
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  criteria,
  file.path(
    REPORT_DIR,
    "12_cross_modal_support_criteria.csv"
  )
)

# ----------------------------- 图形 --------------------------------------------
plot_heatmap <- function(
  data,
  file_path,
  title
) {
  heat_matrix <- as.matrix(
    data[
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
    data$sample_key,
    data$group,
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
    file_path,
    width = 7,
    height = max(
      5,
      0.7 *
        nrow(
          heat_matrix
        ) +
        2
    )
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
      title
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
}

plot_heatmap(
  get_score_set(
    "scRNA_high_specificity"
  ),
  file.path(
    FIG_DIR,
    "Figure_scRNA_high_specificity_stage_heatmap.pdf"
  ),
  "GSE181297 scRNA fixed wound-state scores"
)

plot_heatmap(
  get_score_set(
    "Visium_high_specificity_spots"
  ),
  file.path(
    FIG_DIR,
    "Figure_Visium_high_specificity_stage_heatmap.pdf"
  ),
  "GSE181297 Visium section-level wound-state scores"
)

pdf(
  file.path(
    FIG_DIR,
    "Figure_cross_modal_late_remodeling_support.pdf"
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
  analysis_set in c(
    "scRNA_high_specificity",
    "scRNA_broad_sensitivity",
    "Visium_high_specificity_spots",
    "Visium_broad_top_quartile_spots"
  )
) {
  x <- get_score_set(
    analysis_set
  )

  control_group <- if (
    unique(
      x$modality
    ) ==
      "scRNA"
  ) {
    "normal_scar"
  } else {
    "adjacent_normal"
  }

  position <- ifelse(
    x$group ==
      control_group,
    1,
    2
  )

  y_range <- range(
    x$late_remodeling_state,
    finite = TRUE
  )

  padding <- max(
    0.1,
    diff(
      y_range
    ) *
      0.2
  )

  plot(
    jitter(
      position,
      amount = 0.04
    ),
    x$late_remodeling_state,
    pch = ifelse(
      x$group ==
        "keloid",
      16,
      1
    ),
    xlim = c(
      0.5,
      2.5
    ),
    ylim = c(
      y_range[1] -
        padding,
      y_range[2] +
        padding
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
      control_group,
      "keloid"
    ),
    las = 2
  )

  text(
    jitter(
      position,
      amount = 0.04
    ),
    x$late_remodeling_state,
    labels =
      x$sample_key,
    pos = 3,
    cex = 0.75
  )
}

par(
  op
)

dev.off()

pdf(
  file.path(
    FIG_DIR,
    "Figure_Visium_spot_distribution_descriptive.pdf"
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
    7,
    4,
    3,
    1
  )
)

barplot(
  spot_distribution_summary$high_median_spot_late_remodeling,
  names.arg =
    spot_distribution_summary$sample_id,
  las = 2,
  ylab =
    "Median spot late-remodeling module",
  main =
    "High-specificity fibroblast-enriched spots"
)

barplot(
  spot_distribution_summary$high_fraction_spot_late_positive,
  names.arg =
    spot_distribution_summary$sample_id,
  las = 2,
  ylab =
    "Fraction of selected spots > 0",
  main =
    "Positive late-remodeling spot fraction"
)

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
    structure_audit =
      structure_audit,
    scrna_qc_summary =
      scrna_qc_summary,
    visium_qc_summary =
      visium_qc_summary,
    spot_distribution_summary =
      spot_distribution_summary,
    scRNA_high_specificity =
      scrna_high_result,
    scRNA_broad_sensitivity =
      scrna_broad_result,
    Visium_high_specificity =
      visium_high_result,
    Visium_broad_top_quartile =
      visium_broad_result,
    all_scores =
      all_scores,
    direction_summary =
      direction_summary,
    criteria =
      criteria,
    fixed_signatures =
      signatures,
    frozen_marker_sets =
      marker_sets
  ),
  file = file.path(
    OBJECT_DIR,
    "GSE181297_cross_modal_fixed_state_support.rds"
  ),
  compress = "gzip"
)

# ----------------------------- 最终结论 ----------------------------------------
decision_lines <- c(
  "第十阶段：GSE181297单细胞—Visium跨模态支持结论",
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
  "scRNA高特异性集合：",
  paste0(
    "- late-remodeling mean difference = ",
    round(
      scrna_high_direction$late_remodeling_mean_difference,
      4
    )
  ),
  paste0(
    "- late-remodeling Cliff delta = ",
    round(
      scrna_high_direction$late_remodeling_cliffs_delta,
      3
    )
  ),
  paste0(
    "- ordinal mean difference = ",
    round(
      scrna_high_direction$ordinal_mean_difference,
      4
    )
  ),
  paste0(
    "- complete late separation = ",
    scrna_high_direction$late_complete_separation
  ),
  "",
  "Visium高特异性集合：",
  paste0(
    "- late-remodeling mean difference = ",
    round(
      visium_high_direction$late_remodeling_mean_difference,
      4
    )
  ),
  paste0(
    "- late-remodeling Cliff delta = ",
    round(
      visium_high_direction$late_remodeling_cliffs_delta,
      3
    )
  ),
  paste0(
    "- ordinal mean difference = ",
    round(
      visium_high_direction$ordinal_mean_difference,
      4
    )
  ),
  paste0(
    "- complete late separation = ",
    visium_high_direction$late_complete_separation
  ),
  "",
  "结构限制：",
  paste0(
    "- tissue_positions/scalefactors available = ",
    structure_audit$spatial_coordinates_available
  ),
  "- Visium结果是切片级和spot分布级表达支持，不是空间定位证据。",
  "- scRNA仅2例keloid和1例normal scar，不进行显著性推断。",
  "- Visium仅2例keloid和2例adjacent normal，不把spot当独立重复。",
  "- Pt1/Pt2与NSV1/NSV2不是患者匹配样本。",
  "- 低分辨率图像没有barcode坐标，不能进行组织图像叠加。",
  "- 所有标志物和伤口阶段签名均在本队列分析前冻结。",
  "- ordinal wound-state position不是精确术后日龄。"
)

safe_write_lines(
  decision_lines,
  file.path(
    REPORT_DIR,
    "13_STAGE10_DECISION.txt"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    REPORT_DIR,
    "14_SESSION_INFO.txt"
  )
)

safe_write_lines(
  c(
    paste0(
      "STAGE10_COMPLETED=",
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
      "SCRNA_HIGH_SUPPORT=",
      scrna_high_support
    ),
    paste0(
      "SCRNA_BROAD_SUPPORT=",
      scrna_broad_support
    ),
    paste0(
      "VISIUM_HIGH_SUPPORT=",
      visium_high_support
    ),
    paste0(
      "VISIUM_BROAD_SUPPORT=",
      visium_broad_support
    ),
    paste0(
      "CROSS_MODAL_CONSISTENCY=",
      cross_modal_consistency
    ),
    paste0(
      "SPATIAL_COORDINATES_AVAILABLE=",
      structure_audit$spatial_coordinates_available
    )
  ),
  file.path(
    REPORT_DIR,
    "STAGE10_COMPLETED.txt"
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
  "第十阶段完成"
)

cat(
  "\n============================================================\n"
)

cat(
  "第十阶段运行完成。\n"
)

cat(
  "本地分析对象保存在：",
  file.path(
    OBJECT_DIR,
    "GSE181297_cross_modal_fixed_state_support.rds"
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
