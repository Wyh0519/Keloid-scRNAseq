# =============================================================================
# 项目：瘢痕疙瘩“晚期重塑状态持续与伤口状态终止失败”
# 第八阶段：GSE241124正常伤口空间转录组的固定签名外部验证（V1.1修正版）
#
# 根目录：
#   C:/Users/33652/Desktop/LZZ/1、yssj
#
# 核心目的：
#   1. 读取4名独立供者、Skin/Day1/Day7/Day30共16张Visium切片；
#   2. 使用第六阶段A已经冻结的成纤维细胞标志物定义空间中的
#      fibroblast-enriched spots；
#   3. 使用第四阶段固定、且3次LOODO均稳定的Fibroblast阶段签名；
#   4. 将每张空间切片的fibroblast-enriched pseudobulk投射至
#      正常单细胞伤口状态参照；
#   5. 评价独立空间模态能否恢复正常伤口的序数过程；
#   6. 生成每张切片的fibroblast enrichment与late-remodeling空间图。
#
# 统计边界：
#   - 本阶段是独立供者和独立模态的外部验证；
#   - 不在空间数据中重新选择伤口阶段签名；
#   - 不用空间结果调整前面已经冻结的疾病终点；
#   - Visium spot为混合细胞，结果解释为fibroblast-enriched spatial domain，
#     不能解释为纯成纤维细胞；
#   - 空间中的ordinal position不是精确术后日龄。
#   - V1.1兼容同一空间ZIP内同时存在tissue_positions.csv与
#     tissue_positions_list.csv的标准Space Ranger目录结构。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANG = "en")

ROOT_DIR <- "C:/Users/33652/Desktop/LZZ/1、yssj"

RAW_TAR <- file.path(
  ROOT_DIR,
  "GSE241124_RAW.tar"
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
  "08_STAGE8_GSE241124_SPATIAL_REFERENCE"
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
  "第八阶段_GSE241124正常伤口空间固定签名验证检查包.zip"
)

PACKAGE_TARGZ <- file.path(
  ROOT_DIR,
  "第八阶段_GSE241124正常伤口空间固定签名验证检查包.tar.gz"
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
  "00_stage8_log.txt"
)

WARNING_FILE <- file.path(
  REPORT_DIR,
  "00_stage8_warnings.txt"
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
  "rhdf5",
  quietly = TRUE
)) {
  BiocManager::install(
    "rhdf5",
    ask = FALSE,
    update = FALSE
  )
}

required_packages <- c(
  "Matrix",
  "rhdf5"
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
    )
  )
}

# ----------------------------- 固定样本设计 ------------------------------------
sample_design <- data.frame(
  sample_id = c(
    "Donor4_Skin",
    "Donor4_Wound1",
    "Donor4_Wound7",
    "Donor4_Wound30",
    "Donor1_Skin",
    "Donor1_Wound1",
    "Donor1_Wound7",
    "Donor1_Wound30",
    "Donor2_Skin",
    "Donor2_Wound1",
    "Donor2_Wound7",
    "Donor2_Wound30",
    "Donor3_Skin",
    "Donor3_Wound1",
    "Donor3_Wound7",
    "Donor3_Wound30"
  ),
  file_prefix = c(
    "GSM7717016_Skin",
    "GSM7717017_Wound1",
    "GSM7717018_Wound7",
    "GSM7717019_Wound30",
    "GSM8238462_P17401_1001",
    "GSM8238463_P17401_1002",
    "GSM8238464_P17401_1003",
    "GSM8238465_P17401_1004",
    "GSM8238466_P20063_101",
    "GSM8238467_P20063_102",
    "GSM8238468_P20063_103",
    "GSM8238469_P20063_104",
    "GSM8238470_P20063_105",
    "GSM8238471_P20063_106",
    "GSM8238472_P20063_107",
    "GSM8238473_P20063_108"
  ),
  donor = rep(
    c(
      "Donor4",
      "Donor1",
      "Donor2",
      "Donor3"
    ),
    each = 4
  ),
  condition = rep(
    c(
      "Skin",
      "Wound1",
      "Wound7",
      "Wound30"
    ),
    times = 4
  ),
  stage_code = rep(
    0:3,
    times = 4
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  sample_design,
  file.path(
    REPORT_DIR,
    "01_fixed_spatial_sample_design.csv"
  )
)

# ----------------------------- 基础函数 ----------------------------------------
read_10x_h5 <- function(path) {
  required_nodes <- c(
    "matrix/data",
    "matrix/indices",
    "matrix/indptr",
    "matrix/shape",
    "matrix/barcodes",
    "matrix/features/name",
    "matrix/features/id"
  )

  h5_listing <- rhdf5::h5ls(
    path,
    recursive = TRUE
  )

  available_nodes <- paste0(
    h5_listing$group,
    "/",
    h5_listing$name
  )

  available_nodes <- sub(
    "^//",
    "/",
    available_nodes
  )

  available_nodes <- sub(
    "^/",
    "",
    available_nodes
  )

  missing_nodes <- setdiff(
    required_nodes,
    available_nodes
  )

  if (length(missing_nodes)) {
    stop(
      "H5缺少10x节点：",
      paste(
        missing_nodes,
        collapse = ", "
      ),
      " | 文件：",
      basename(path)
    )
  }

  data_values <- as.numeric(
    rhdf5::h5read(
      path,
      "matrix/data"
    )
  )

  row_indices <- as.integer(
    rhdf5::h5read(
      path,
      "matrix/indices"
    )
  )

  column_pointer <- as.integer(
    rhdf5::h5read(
      path,
      "matrix/indptr"
    )
  )

  matrix_shape <- as.integer(
    rhdf5::h5read(
      path,
      "matrix/shape"
    )
  )

  barcodes <- as.character(
    rhdf5::h5read(
      path,
      "matrix/barcodes"
    )
  )

  gene_symbol <- as.character(
    rhdf5::h5read(
      path,
      "matrix/features/name"
    )
  )

  gene_id <- as.character(
    rhdf5::h5read(
      path,
      "matrix/features/id"
    )
  )

  if (
    length(matrix_shape) != 2L ||
      length(column_pointer) != matrix_shape[2] + 1L ||
      length(barcodes) != matrix_shape[2] ||
      length(gene_symbol) != matrix_shape[1] ||
      length(gene_id) != matrix_shape[1] ||
      length(data_values) != length(row_indices)
  ) {
    stop(
      "H5内部维度不一致：",
      basename(path)
    )
  }

  matrix_result <- methods::new(
    "dgCMatrix",
    i = row_indices,
    p = column_pointer,
    x = data_values,
    Dim = matrix_shape,
    Dimnames = list(
      make.unique(
        gene_symbol
      ),
      barcodes
    )
  )

  list(
    counts = matrix_result,
    gene_id = gene_id,
    gene_symbol = gene_symbol,
    barcodes = barcodes
  )
}

read_tissue_positions <- function(path) {
  first_line <- readLines(
    path,
    n = 1L,
    warn = FALSE
  )

  has_header <- grepl(
    "barcode",
    first_line,
    ignore.case = TRUE
  )

  if (has_header) {
    positions <- utils::read.csv(
      path,
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    normalized_names <- tolower(
      gsub(
        "[^a-z0-9]",
        "",
        names(
          positions
        )
      )
    )

    find_column <- function(candidates) {
      candidate_normalized <- tolower(
        gsub(
          "[^a-z0-9]",
          "",
          candidates
        )
      )

      idx <- match(
        candidate_normalized,
        normalized_names
      )

      idx <- idx[
        !is.na(
          idx
        )
      ]

      if (!length(idx)) {
        return(
          NA_integer_
        )
      }

      idx[1]
    }

    barcode_col <- find_column(
      c(
        "barcode"
      )
    )

    tissue_col <- find_column(
      c(
        "in_tissue",
        "intissue"
      )
    )

    array_row_col <- find_column(
      c(
        "array_row",
        "arrayrow"
      )
    )

    array_col_col <- find_column(
      c(
        "array_col",
        "arraycol"
      )
    )

    pixel_row_col <- find_column(
      c(
        "pxl_row_in_fullres",
        "pxlrowinfullres",
        "imagerow"
      )
    )

    pixel_col_col <- find_column(
      c(
        "pxl_col_in_fullres",
        "pxlcolinfullres",
        "imagecol"
      )
    )

    required_idx <- c(
      barcode_col,
      tissue_col,
      array_row_col,
      array_col_col,
      pixel_row_col,
      pixel_col_col
    )

    if (anyNA(required_idx)) {
      stop(
        "无法解析带表头的tissue_positions文件：",
        basename(path)
      )
    }

    result <- data.frame(
      barcode = as.character(
        positions[[barcode_col]]
      ),
      in_tissue = as.integer(
        positions[[tissue_col]]
      ),
      array_row = as.integer(
        positions[[array_row_col]]
      ),
      array_col = as.integer(
        positions[[array_col_col]]
      ),
      pxl_row_in_fullres = as.numeric(
        positions[[pixel_row_col]]
      ),
      pxl_col_in_fullres = as.numeric(
        positions[[pixel_col_col]]
      ),
      stringsAsFactors = FALSE
    )
  } else {
    positions <- utils::read.csv(
      path,
      header = FALSE,
      stringsAsFactors = FALSE
    )

    if (ncol(positions) < 6L) {
      stop(
        "无表头tissue_positions文件列数少于6：",
        basename(path)
      )
    }

    result <- data.frame(
      barcode = as.character(
        positions[[1]]
      ),
      in_tissue = as.integer(
        positions[[2]]
      ),
      array_row = as.integer(
        positions[[3]]
      ),
      array_col = as.integer(
        positions[[4]]
      ),
      pxl_row_in_fullres = as.numeric(
        positions[[5]]
      ),
      pxl_col_in_fullres = as.numeric(
        positions[[6]]
      ),
      stringsAsFactors = FALSE
    )
  }

  result
}

locate_sample_files <- function(
  outer_dir,
  file_prefix
) {
  files <- list.files(
    outer_dir,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  )

  h5_candidates <- files[
    grepl(
      file_prefix,
      basename(
        files
      ),
      fixed = TRUE
    ) &
      grepl(
        "\\.h5$",
        files,
        ignore.case = TRUE
      )
  ]

  zip_candidates <- files[
    grepl(
      file_prefix,
      basename(
        files
      ),
      fixed = TRUE
    ) &
      grepl(
        "\\.zip$",
        files,
        ignore.case = TRUE
      )
  ]

  if (
    length(h5_candidates) != 1L ||
      length(zip_candidates) != 1L
  ) {
    stop(
      "无法唯一定位空间文件：",
      file_prefix,
      " | H5=",
      length(h5_candidates),
      " | ZIP=",
      length(zip_candidates)
    )
  }

  list(
    h5 = h5_candidates,
    spatial_zip = zip_candidates
  )
}

extract_spatial_zip <- function(
  zip_path,
  output_dir
) {
  if (dir.exists(output_dir)) {
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

  # 排除macOS压缩包可能附带的隐藏副本。
  valid_files <- files[
    !grepl(
      "(^|[/\\\\])__MACOSX([/\\\\]|$)",
      files,
      ignore.case = TRUE,
      perl = TRUE
    ) &
      !startsWith(
        basename(files),
        "._"
      )
  ]

  position_candidates <- valid_files[
    grepl(
      "^tissue_positions(_list)?\\.csv$",
      basename(
        valid_files
      ),
      ignore.case = TRUE,
      perl = TRUE
    )
  ]

  if (!length(position_candidates)) {
    # 对少数带前缀的文件名保留一次宽松回退。
    position_candidates <- valid_files[
      grepl(
        "tissue_positions.*\\.csv$",
        basename(
          valid_files
        ),
        ignore.case = TRUE,
        perl = TRUE
      )
    ]
  }

  if (!length(position_candidates)) {
    stop(
      "空间ZIP中未找到tissue_positions坐标文件：",
      basename(zip_path)
    )
  }

  candidate_basenames <- tolower(
    basename(
      position_candidates
    )
  )

  # Space Ranger新版本通常同时包含：
  #   tissue_positions.csv        （带表头，优先）
  #   tissue_positions_list.csv   （旧格式兼容文件）
  # 两者并存是合法结构，不应报错。
  exact_header <- which(
    candidate_basenames ==
      "tissue_positions.csv"
  )

  exact_legacy <- which(
    candidate_basenames ==
      "tissue_positions_list.csv"
  )

  if (length(exact_header)) {
    selected_candidates <-
      position_candidates[
        exact_header
      ]

    selection_strategy <-
      "preferred_tissue_positions.csv"
  } else if (length(exact_legacy)) {
    selected_candidates <-
      position_candidates[
        exact_legacy
      ]

    selection_strategy <-
      "fallback_tissue_positions_list.csv"
  } else {
    selected_candidates <-
      position_candidates

    selection_strategy <-
      "fallback_pattern_match"
  }

  # 若ZIP内存在同名嵌套副本，优先路径层级最浅、路径最短者。
  path_depth <- lengths(
    strsplit(
      normalizePath(
        selected_candidates,
        winslash = "/",
        mustWork = FALSE
      ),
      "/",
      fixed = TRUE
    )
  )

  path_length <- nchar(
    normalizePath(
      selected_candidates,
      winslash = "/",
      mustWork = FALSE
    )
  )

  selected_order <- order(
    path_depth,
    path_length,
    normalizePath(
      selected_candidates,
      winslash = "/",
      mustWork = FALSE
    )
  )

  selected_position <-
    selected_candidates[
      selected_order[1]
    ]

  if (length(selected_candidates) > 1L) {
    warning(
      "空间ZIP含多个同优先级坐标文件；已按最浅路径确定性选择：",
      basename(selected_position),
      " | ZIP=",
      basename(zip_path),
      call. = FALSE
    )
  }

  list(
    positions =
      selected_position,
    all_files =
      valid_files,
    position_candidates =
      position_candidates,
    selected_position =
      selected_position,
    selection_strategy =
      selection_strategy
  )
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

  if (length(missing_samples)) {
    stop(
      "空间pseudobulk缺少样本：",
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

  if (length(common_genes) < 10000L) {
    warn_msg(
      "空间切片间共同基因少于10,000：",
      length(common_genes)
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

  if (anyNA(idx)) {
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

  if (any(n_genes < 20L)) {
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

sum_gene_counts <- function(
  counts,
  gene_symbol,
  genes
) {
  gene_upper <- toupper(
    trimws(
      as.character(
        gene_symbol
      )
    )
  )

  genes <- toupper(
    genes
  )

  idx <- which(
    gene_upper %in%
      genes
  )

  if (!length(idx)) {
    return(
      matrix(
        0,
        nrow = length(
          genes
        ),
        ncol = ncol(
          counts
        ),
        dimnames = list(
          genes,
          colnames(
            counts
          )
        )
      )
    )
  }

  output <- matrix(
    0,
    nrow = length(
      genes
    ),
    ncol = ncol(
      counts
    ),
    dimnames = list(
      genes,
      colnames(
        counts
      )
    )
  )

  for (
    i in seq_along(
      genes
    )
  ) {
    gene_idx <- which(
      gene_upper ==
        genes[i]
    )

    if (length(gene_idx)) {
      output[i, ] <- Matrix::colSums(
        counts[
          gene_idx,
          ,
          drop = FALSE
        ]
      )
    }
  }

  output
}

calculate_marker_metrics <- function(
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
    marker_counts <- sum_gene_counts(
      counts,
      gene_symbol,
      marker_sets[[marker_name]]
    )

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

    result[[marker_name]] <- list(
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

  result
}

calculate_spot_stage_modules <- function(
  counts,
  gene_symbol,
  signatures,
  library_size
) {
  stage_scores <- matrix(
    NA_real_,
    nrow = ncol(
      counts
    ),
    ncol = length(
      signatures
    ),
    dimnames = list(
      colnames(
        counts
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
    signature_counts <- sum_gene_counts(
      counts,
      gene_symbol,
      signatures[[stage]]
    )

    normalized <- log1p(
      sweep(
        signature_counts,
        2,
        pmax(
          library_size,
          1
        ),
        "/"
      ) *
        10000
    )

    stage_scores[, stage] <- colMeans(
      normalized
    )
  }

  stage_scores
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
        ) ^ 2
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

  if (anyNA(idx)) {
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

calculate_validation_metrics <- function(
  mapping,
  sample_design,
  analysis_set
) {
  x <- merge(
    mapping,
    sample_design[
      ,
      c(
        "sample_id",
        "donor",
        "condition",
        "stage_code"
      )
    ],
    by.x = "sample_key",
    by.y = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )

  x <- x[
    match(
      mapping$sample_key,
      x$sample_key
    ),
    ,
    drop = FALSE
  ]

  predicted_code <- c(
    Skin = 0,
    Wound1 = 1,
    Wound7 = 2,
    Wound30 = 3
  )[
    x$nearest_stage
  ]

  exact_correct <-
    predicted_code ==
    x$stage_code

  adjacent_correct <-
    abs(
      predicted_code -
        x$stage_code
    ) <=
    1

  overall <- data.frame(
    analysis_set =
      analysis_set,
    n_sections =
      nrow(
        x
      ),
    exact_accuracy =
      mean(
        exact_correct
      ),
    adjacent_accuracy =
      mean(
        adjacent_correct
      ),
    mean_ordinal_error =
      mean(
        abs(
          predicted_code -
            x$stage_code
        )
      ),
    pooled_spearman =
      suppressWarnings(
        stats::cor(
          x$stage_code,
          x$ordinal_wound_state_position,
          method = "spearman"
        )
      ),
    stringsAsFactors = FALSE
  )

  donor_rows <- list()

  for (
    donor in unique(
      x$donor
    )
  ) {
    xd <- x[
      x$donor ==
        donor,
      ,
      drop = FALSE
    ]

    predicted_donor_code <- c(
      Skin = 0,
      Wound1 = 1,
      Wound7 = 2,
      Wound30 = 3
    )[
      xd$nearest_stage
    ]

    donor_rows[[length(
      donor_rows
    ) + 1L]] <- data.frame(
      analysis_set =
        analysis_set,
      donor =
        donor,
      exact_accuracy =
        mean(
          predicted_donor_code ==
            xd$stage_code
        ),
      adjacent_accuracy =
        mean(
          abs(
            predicted_donor_code -
              xd$stage_code
          ) <=
            1
        ),
      mean_ordinal_error =
        mean(
          abs(
            predicted_donor_code -
              xd$stage_code
          )
        ),
      spearman =
        suppressWarnings(
          stats::cor(
            xd$stage_code,
            xd$ordinal_wound_state_position,
            method =
              "spearman"
          )
        ),
      wound30_has_highest_late_remodeling =
        xd$late_remodeling_state[
          xd$condition ==
            "Wound30"
        ] ==
        max(
          xd$late_remodeling_state,
          na.rm = TRUE
        ),
      skin_has_highest_skin_return =
        xd$skin_return_state[
          xd$condition ==
            "Skin"
        ] ==
        max(
          xd$skin_return_state,
          na.rm = TRUE
        ),
      stringsAsFactors = FALSE
    )
  }

  donor_metrics <- do.call(
    rbind,
    donor_rows
  )

  list(
    mapping_with_truth =
      x,
    overall =
      overall,
    donor =
      donor_metrics
  )
}

# ----------------------------- 输入检查 ----------------------------------------
log_msg(
  "第八阶段开始"
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

if (length(missing_files)) {
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

marker_definition <- marker_object$marker_definition

safe_write_csv(
  marker_definition,
  file.path(
    REPORT_DIR,
    "02_reused_frozen_marker_definition.csv"
  )
)

signatures <- get_stable_fibroblast_signatures(
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
    "正常单细胞参考的Fibroblast pseudobulk不是12个。"
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

# ----------------------------- 解包RAW -----------------------------------------
outer_dir <- file.path(
  WORK_DIR,
  "GSE241124_RAW"
)

dir.create(
  outer_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

log_msg(
  "解包GSE241124_RAW.tar"
)

utils::untar(
  RAW_TAR,
  exdir = outer_dir
)

# ----------------------------- 逐切片处理 --------------------------------------
qc_rows <- list()
spot_selection_rows <- list()
feature_rows <- list()
position_file_rows <- list()
high_pb_list <- list()
broad_pb_list <- list()
spot_score_objects <- list()

pdf(
  file.path(
    FIG_DIR,
    "Figure_spatial_QC_and_fibroblast_enrichment.pdf"
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
  sample_row <- sample_design[
    i,
    ,
    drop = FALSE
  ]

  sample_id <- sample_row$sample_id
  file_prefix <- sample_row$file_prefix

  log_msg(
    "处理空间切片：",
    sample_id
  )

  sample_files <- locate_sample_files(
    outer_dir,
    file_prefix
  )

  sample_spatial_dir <- file.path(
    WORK_DIR,
    paste0(
      "spatial_",
      sample_id
    )
  )

  spatial_files <- extract_spatial_zip(
    sample_files$spatial_zip,
    sample_spatial_dir
  )

  position_file_rows[[length(
    position_file_rows
  ) + 1L]] <- data.frame(
    sample_id =
      sample_id,
    spatial_zip =
      basename(
        sample_files$spatial_zip
      ),
    n_position_candidates =
      length(
        spatial_files$position_candidates
      ),
    position_candidates =
      paste(
        basename(
          spatial_files$position_candidates
        ),
        collapse = " | "
      ),
    selected_position_file =
      basename(
        spatial_files$selected_position
      ),
    selection_strategy =
      spatial_files$selection_strategy,
    stringsAsFactors = FALSE
  )

  h5_object <- read_10x_h5(
    sample_files$h5
  )

  counts <- h5_object$counts
  gene_id <- h5_object$gene_id
  gene_symbol <- h5_object$gene_symbol

  positions <- read_tissue_positions(
    spatial_files$positions
  )

  position_idx <- match(
    colnames(
      counts
    ),
    positions$barcode
  )

  if (anyNA(position_idx)) {
    warn_msg(
      sample_id,
      "有",
      sum(
        is.na(
          position_idx
        )
      ),
      "个H5 barcode未匹配tissue_positions。"
    )
  }

  position_aligned <- positions[
    position_idx,
    ,
    drop = FALSE
  ]

  in_tissue <- (
    !is.na(
      position_aligned$in_tissue
    ) &
      position_aligned$in_tissue ==
      1L
  )

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

  qc_keep <- (
    in_tissue &
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
      "QC后in-tissue spots少于100：",
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

  positions_qc <- position_aligned[
    qc_keep,
    ,
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

  marker_metrics <- calculate_marker_metrics(
    counts_qc,
    gene_symbol,
    marker_sets,
    n_count_qc
  )

  max_exclusion_score <- pmax(
    marker_metrics$immune$score,
    marker_metrics$epithelial$score,
    marker_metrics$endothelial$score,
    marker_metrics$neural$score
  )

  fibro_identity_score <-
    marker_metrics$fibro_core$score +
    0.5 *
    marker_metrics$fibro_support$score -
    max_exclusion_score

  high_spot <- (
    marker_metrics$fibro_core$detected >=
      3 &
      marker_metrics$fibro_support$detected >=
      1 &
      marker_metrics$fibro_core$score >
      max_exclusion_score +
      0.15
  )

  high_spot[
    is.na(
      high_spot
    )
  ] <- FALSE

  broad_candidates <- (
    marker_metrics$fibro_core$detected >=
      2 &
      fibro_identity_score >
      0
  )

  broad_candidates[
    is.na(
      broad_candidates
    )
  ] <- FALSE

  candidate_identity_values <-
    fibro_identity_score[
      broad_candidates &
        is.finite(
          fibro_identity_score
        )
    ]

  if (length(
    candidate_identity_values
  )) {
    broad_cutoff <- as.numeric(
      stats::quantile(
        candidate_identity_values,
        probs = 0.75,
        na.rm = TRUE,
        names = FALSE,
        type = 8
      )
    )
  } else {
    finite_identity_values <-
      fibro_identity_score[
        is.finite(
          fibro_identity_score
        )
      ]

    if (!length(
      finite_identity_values
    )) {
      stop(
        sample_id,
        "没有可用于空间成纤维细胞富集的有限marker score。"
      )
    }

    broad_cutoff <- as.numeric(
      stats::quantile(
        finite_identity_values,
        probs = 0.75,
        na.rm = TRUE,
        names = FALSE,
        type = 8
      )
    )

    warn_msg(
      sample_id,
      "没有spot满足预设broad candidate条件；",
      "broad cutoff改用全部有限fibro identity score的75%分位数，",
      "但仍保留fibro_core_detected>=2和identity>0条件。"
    )
  }

  broad_spot <- (
    broad_candidates &
      fibro_identity_score >=
      broad_cutoff
  )

  broad_spot[
    is.na(
      broad_spot
    )
  ] <- FALSE

  if (
    sum(
      high_spot
    ) <
      50L
  ) {
    warn_msg(
      sample_id,
      "高特异性fibroblast-enriched spots少于50：",
      sum(
        high_spot
      )
    )
  }

  if (
    sum(
      broad_spot
    ) <
      50L
  ) {
    warn_msg(
      sample_id,
      "broad fibroblast-enriched spots少于50：",
      sum(
        broad_spot
      )
    )
  }

  if (
    sum(
      high_spot
    ) <
      20L
  ) {
    stop(
      sample_id,
      "高特异性fibroblast-enriched spots少于20，无法建立稳定pseudobulk。"
    )
  }

  if (
    sum(
      broad_spot
    ) <
      20L
  ) {
    stop(
      sample_id,
      "broad fibroblast-enriched spots少于20，无法建立稳定pseudobulk。"
    )
  }

  high_pb_raw <- Matrix::rowSums(
    counts_qc[
      ,
      high_spot,
      drop = FALSE
    ]
  )

  broad_pb_raw <- Matrix::rowSums(
    counts_qc[
      ,
      broad_spot,
      drop = FALSE
    ]
  )

  high_pb_list[[sample_id]] <-
    collapse_vector_by_symbol(
      high_pb_raw,
      gene_symbol
    )

  broad_pb_list[[sample_id]] <-
    collapse_vector_by_symbol(
      broad_pb_raw,
      gene_symbol
    )

  stage_modules <- calculate_spot_stage_modules(
    counts_qc,
    gene_symbol,
    signatures,
    n_count_qc
  )

  spot_late_remodeling <-
    stage_modules[, "Wound30"] -
    (
      stage_modules[, "Wound1"] +
        stage_modules[, "Wound7"]
    ) /
    2

  spot_score_objects[[sample_id]] <- data.frame(
    sample_id = sample_id,
    donor = sample_row$donor,
    condition = sample_row$condition,
    barcode = colnames(
      counts_qc
    ),
    pxl_row_in_fullres =
      positions_qc$pxl_row_in_fullres,
    pxl_col_in_fullres =
      positions_qc$pxl_col_in_fullres,
    nCount = n_count_qc,
    nFeature = n_feature_qc,
    percent_mt = pct_mt_qc,
    fibro_identity_score =
      fibro_identity_score,
    high_fibroblast_enriched =
      high_spot,
    broad_fibroblast_enriched =
      broad_spot,
    module_Skin =
      stage_modules[, "Skin"],
    module_Wound1 =
      stage_modules[, "Wound1"],
    module_Wound7 =
      stage_modules[, "Wound7"],
    module_Wound30 =
      stage_modules[, "Wound30"],
    spot_late_remodeling =
      spot_late_remodeling,
    stringsAsFactors = FALSE
  )

  qc_rows[[length(
    qc_rows
  ) + 1L]] <- data.frame(
    sample_id = sample_id,
    donor = sample_row$donor,
    condition = sample_row$condition,
    n_h5_spots = ncol(
      counts
    ),
    n_in_tissue_spots = sum(
      in_tissue,
      na.rm = TRUE
    ),
    n_QC_spots = ncol(
      counts_qc
    ),
    QC_retention_from_in_tissue =
      ncol(
        counts_qc
      ) /
      sum(
        in_tissue,
        na.rm = TRUE
      ),
    n_high_fibroblast_enriched =
      sum(
        high_spot
      ),
    n_broad_fibroblast_enriched =
      sum(
        broad_spot
      ),
    high_fraction =
      mean(
        high_spot
      ),
    broad_fraction =
      mean(
        broad_spot
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
    median_fibro_identity =
      stats::median(
        fibro_identity_score
      ),
    broad_cutoff =
      broad_cutoff,
    stringsAsFactors = FALSE
  )

  spot_selection_rows[[length(
    spot_selection_rows
  ) + 1L]] <- data.frame(
    sample_id =
      sample_id,
    donor =
      sample_row$donor,
    condition =
      sample_row$condition,
    selection_set = c(
      "high_specificity",
      "broad_top_quartile"
    ),
    n_selected = c(
      sum(
        high_spot
      ),
      sum(
        broad_spot
      )
    ),
    median_fibro_identity = c(
      stats::median(
        fibro_identity_score[
          high_spot
        ]
      ),
      stats::median(
        fibro_identity_score[
          broad_spot
        ]
      )
    ),
    median_spot_late_remodeling = c(
      stats::median(
        spot_late_remodeling[
          high_spot
        ]
      ),
      stats::median(
        spot_late_remodeling[
          broad_spot
        ]
      )
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
        gene_symbol
      ),
    n_unique_gene_symbols =
      length(
        unique(
          gene_symbol
        )
      ),
    stringsAsFactors = FALSE
  )

  # 审计图：QC、fibroblast identity和late-remodeling空间分布。
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
    fibro_identity_score,
    breaks = 50,
    xlab =
      "Fibroblast identity score",
    main = paste0(
      sample_id,
      " fibroblast enrichment"
    )
  )

  plot(
    positions_qc$pxl_col_in_fullres,
    -positions_qc$pxl_row_in_fullres,
    pch = 16,
    cex = 0.35,
    xlab = "Spatial x",
    ylab = "Spatial y",
    main = paste0(
      sample_id,
      " high fibroblast spots"
    )
  )

  points(
    positions_qc$pxl_col_in_fullres[
      high_spot
    ],
    -positions_qc$pxl_row_in_fullres[
      high_spot
    ],
    pch = 16,
    cex = 0.45
  )

  plot(
    positions_qc$pxl_col_in_fullres,
    -positions_qc$pxl_row_in_fullres,
    pch = 1,
    cex = 0.25,
    xlab = "Spatial x",
    ylab = "Spatial y",
    main = paste0(
      sample_id,
      " late-remodeling module"
    )
  )

  if (sum(
    broad_spot
  ) >
      0L) {
    late_values <-
      spot_late_remodeling[
        broad_spot
      ]

    late_rank <- rank(
      late_values,
      ties.method = "average"
    ) /
      length(
        late_values
      )

    point_cex <- 0.35 +
      0.65 *
      late_rank

    points(
      positions_qc$pxl_col_in_fullres[
        broad_spot
      ],
      -positions_qc$pxl_row_in_fullres[
        broad_spot
      ],
      pch = 16,
      cex = point_cex
    )
  }

  par(
    op
  )

  rm(
    counts,
    counts_qc,
    stage_modules
  )

  gc(
    verbose = FALSE
  )

  unlink(
    sample_spatial_dir,
    recursive = TRUE,
    force = TRUE
  )
}

dev.off()

if (!length(
  position_file_rows
)) {
  stop(
    "未生成任何空间坐标文件选择记录。"
  )
}

position_file_summary <- do.call(
  rbind,
  position_file_rows
)

rownames(
  position_file_summary
) <- NULL

safe_write_csv(
  position_file_summary,
  file.path(
    REPORT_DIR,
    "04A_spatial_position_file_selection.csv"
  )
)

qc_summary <- do.call(
  rbind,
  qc_rows
)

spot_selection_summary <- do.call(
  rbind,
  spot_selection_rows
)

feature_summary <- do.call(
  rbind,
  feature_rows
)

spot_scores <- do.call(
  rbind,
  spot_score_objects
)

rownames(
  qc_summary
) <- NULL

rownames(
  spot_selection_summary
) <- NULL

rownames(
  feature_summary
) <- NULL

rownames(
  spot_scores
) <- NULL

safe_write_csv(
  qc_summary,
  file.path(
    REPORT_DIR,
    "04_spatial_QC_summary.csv"
  )
)

safe_write_csv(
  spot_selection_summary,
  file.path(
    REPORT_DIR,
    "05_fibroblast_enriched_spot_summary.csv"
  )
)

safe_write_csv(
  feature_summary,
  file.path(
    REPORT_DIR,
    "06_spatial_feature_summary.csv"
  )
)

# spot-level文件较大，但检查包需要审计；
# 仅保存每个样本前2000个spots，完整对象保存在RDS中。
spot_scores_preview <- do.call(
  rbind,
  lapply(
    split(
      spot_scores,
      spot_scores$sample_id
    ),
    function(x) {
      head(
        x,
        2000L
      )
    }
  )
)

rownames(
  spot_scores_preview
) <- NULL

safe_write_csv(
  spot_scores_preview,
  file.path(
    REPORT_DIR,
    "07_spot_score_preview_max2000_per_section.csv"
  )
)

# ----------------------------- 空间pseudobulk映射 ------------------------------
analyse_spatial_set <- function(
  set_name,
  pb_list,
  reference_counts,
  reference_design,
  signatures,
  sample_design
) {
  spatial_counts <- merge_named_vector_list(
    pb_list,
    required_names =
      sample_design$sample_id
  )

  common_genes <- intersect(
    rownames(
      reference_counts
    ),
    rownames(
      spatial_counts
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
      "与单细胞参考的共同背景基因少于5000：",
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

  spatial_expression <- log_cpm(
    spatial_counts[
      common_genes,
      ,
      drop = FALSE
    ]
  )

  reference_scores <- score_samples_by_rank(
    reference_expression,
    signatures
  )

  spatial_scores <- score_samples_by_rank(
    spatial_expression,
    signatures
  )

  standardized <- standardize_by_reference(
    reference_scores,
    spatial_scores
  )

  reference_z <- standardized$reference_z
  spatial_z <- standardized$new_z

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
    spatial_z,
    centroids,
    reference_distance_limit
  )

  mapping <- add_fixed_contrasts(
    mapping,
    spatial_z
  )

  mapping$analysis_set <-
    set_name

  validation <- calculate_validation_metrics(
    mapping,
    sample_design,
    set_name
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
      spatial_counts,
    mapping =
      validation$mapping_with_truth,
    overall =
      validation$overall,
    donor =
      validation$donor,
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

high_result <- analyse_spatial_set(
  set_name =
    "high_specificity_spots",
  pb_list =
    high_pb_list,
  reference_counts =
    reference_counts,
  reference_design =
    reference_design,
  signatures =
    signatures,
  sample_design =
    sample_design
)

broad_result <- analyse_spatial_set(
  set_name =
    "broad_top_quartile_spots",
  pb_list =
    broad_pb_list,
  reference_counts =
    reference_counts,
  reference_design =
    reference_design,
  signatures =
    signatures,
  sample_design =
    sample_design
)

section_scores <- rbind(
  high_result$mapping,
  broad_result$mapping
)

overall_metrics <- rbind(
  high_result$overall,
  broad_result$overall
)

donor_metrics <- rbind(
  high_result$donor,
  broad_result$donor
)

coverage <- rbind(
  high_result$coverage,
  broad_result$coverage
)

rownames(
  section_scores
) <- NULL

rownames(
  overall_metrics
) <- NULL

rownames(
  donor_metrics
) <- NULL

rownames(
  coverage
) <- NULL

safe_write_csv(
  coverage,
  file.path(
    REPORT_DIR,
    "08_signature_gene_coverage.csv"
  )
)

safe_write_csv(
  section_scores,
  file.path(
    REPORT_DIR,
    "09_spatial_section_wound_state_scores.csv"
  )
)

safe_write_csv(
  overall_metrics,
  file.path(
    REPORT_DIR,
    "10_spatial_validation_overall_metrics.csv"
  )
)

safe_write_csv(
  donor_metrics,
  file.path(
    REPORT_DIR,
    "11_spatial_validation_donor_metrics.csv"
  )
)

# ----------------------------- 阶段描述统计 ------------------------------------
stage_summary <- do.call(
  rbind,
  lapply(
    split(
      section_scores,
      section_scores$analysis_set
    ),
    function(x) {
      do.call(
        rbind,
        lapply(
          c(
            "Skin",
            "Wound1",
            "Wound7",
            "Wound30"
          ),
          function(stage) {
            xs <- x[
              x$condition ==
                stage,
              ,
              drop = FALSE
            ]

            data.frame(
              analysis_set =
                unique(
                  x$analysis_set
                ),
              condition =
                stage,
              n_sections =
                nrow(
                  xs
                ),
              mean_ordinal_position =
                mean(
                  xs$ordinal_wound_state_position
                ),
              median_ordinal_position =
                stats::median(
                  xs$ordinal_wound_state_position
                ),
              mean_late_remodeling =
                mean(
                  xs$late_remodeling_state
                ),
              median_late_remodeling =
                stats::median(
                  xs$late_remodeling_state
                ),
              mean_z_Wound30 =
                mean(
                  xs$z_Wound30
                ),
              mean_skin_return =
                mean(
                  xs$skin_return_state
                ),
              stringsAsFactors = FALSE
            )
          }
        )
      )
    }
  )
)

rownames(
  stage_summary
) <- NULL

safe_write_csv(
  stage_summary,
  file.path(
    REPORT_DIR,
    "12_stage_level_spatial_summary.csv"
  )
)

# ----------------------------- 决策闸门 ----------------------------------------
get_overall_row <- function(
  analysis_set
) {
  x <- overall_metrics[
    overall_metrics$analysis_set ==
      analysis_set,
    ,
    drop = FALSE
  ]

  if (nrow(
    x
  ) !=
      1L) {
    stop(
      "无法唯一获得空间验证指标：",
      analysis_set
    )
  }

  x
}

high_overall <- get_overall_row(
  "high_specificity_spots"
)

broad_overall <- get_overall_row(
  "broad_top_quartile_spots"
)

high_donor <- donor_metrics[
  donor_metrics$analysis_set ==
    "high_specificity_spots",
  ,
  drop = FALSE
]

broad_donor <- donor_metrics[
  donor_metrics$analysis_set ==
    "broad_top_quartile_spots",
  ,
  drop = FALSE
]

high_spot_count_pass <- all(
  qc_summary$n_high_fibroblast_enriched >=
    50
)

broad_spot_count_pass <- all(
  qc_summary$n_broad_fibroblast_enriched >=
    50
)

high_strong <- (
  high_overall$exact_accuracy >=
    0.625 &&
    high_overall$adjacent_accuracy >=
    0.875 &&
    high_overall$pooled_spearman >=
    0.70 &&
    sum(
      high_donor$wound30_has_highest_late_remodeling
    ) >=
    3L &&
    sum(
      high_donor$spearman >
        0
    ) >=
    3L
)

broad_concordant <- (
  broad_overall$exact_accuracy >=
    0.50 &&
    broad_overall$adjacent_accuracy >=
    0.75 &&
    broad_overall$pooled_spearman >=
    0.50 &&
    sum(
      broad_donor$wound30_has_highest_late_remodeling
    ) >=
    2L
)

high_moderate <- (
  high_overall$adjacent_accuracy >=
    0.75 &&
    high_overall$pooled_spearman >=
    0.50 &&
    sum(
      high_donor$wound30_has_highest_late_remodeling
    ) >=
    2L
)

if (
  high_spot_count_pass &&
    broad_spot_count_pass &&
    high_strong &&
    broad_concordant
) {
  project_gate <-
    "PASS_STRONG_SPATIAL_REFERENCE_VALIDATION"

  interpretation <- paste(
    "固定Fibroblast伤口状态签名在4名独立供者和空间转录组模态中",
    "恢复了正常伤口Skin至Day30的序数过程，",
    "且Day30切片在多数供者中具有最高late-remodeling state。"
  )
} else if (
  high_spot_count_pass &&
    high_moderate &&
    broad_concordant
) {
  project_gate <-
    "PASS_ORDINAL_SPATIAL_VALIDATION"

  interpretation <- paste(
    "空间数据支持正常伤口状态的总体序数变化，",
    "但精确四阶段分类未达到强验证标准。"
  )
} else if (
  high_spot_count_pass &&
    high_moderate
) {
  project_gate <-
    "PARTIAL_SPATIAL_SUPPORT"

  interpretation <- paste(
    "高特异性fibroblast-enriched spots提供部分序数支持，",
    "但宽松空间定义或Day30峰值一致性不足。"
  )
} else {
  project_gate <-
    "FAIL_SPATIAL_REFERENCE_VALIDATION"

  interpretation <- paste(
    "固定单细胞Fibroblast阶段签名未在空间模态中稳定恢复正常伤口序数过程。"
  )
}

decision_criteria <- data.frame(
  criterion = c(
    "high_spot_count_pass",
    "broad_spot_count_pass",
    "high_strong_validation",
    "broad_concordant_validation",
    "high_moderate_validation"
  ),
  passed = c(
    high_spot_count_pass,
    broad_spot_count_pass,
    high_strong,
    broad_concordant,
    high_moderate
  ),
  stringsAsFactors = FALSE
)

safe_write_csv(
  decision_criteria,
  file.path(
    REPORT_DIR,
    "13_spatial_validation_criteria.csv"
  )
)

# ----------------------------- 汇总图 ------------------------------------------
pdf(
  file.path(
    FIG_DIR,
    "Figure_spatial_section_ordinal_validation.pdf"
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
    5,
    4,
    3,
    1
  )
)

for (
  analysis_set in c(
    "high_specificity_spots",
    "broad_top_quartile_spots"
  )
) {
  x <- section_scores[
    section_scores$analysis_set ==
      analysis_set,
    ,
    drop = FALSE
  ]

  plot(
    x$stage_code,
    x$ordinal_wound_state_position,
    pch = 16,
    xlab =
      "True wound stage code",
    ylab =
      "Predicted ordinal position",
    main =
      analysis_set,
    xaxt = "n"
  )

  axis(
    1,
    at = 0:3,
    labels = c(
      "Skin",
      "D1",
      "D7",
      "D30"
    )
  )

  text(
    x$stage_code,
    x$ordinal_wound_state_position,
    labels =
      x$donor,
    pos = 3,
    cex = 0.65
  )

  abline(
    a = 0,
    b = 1,
    lty = 2
  )
}

for (
  analysis_set in c(
    "high_specificity_spots",
    "broad_top_quartile_spots"
  )
) {
  x <- section_scores[
    section_scores$analysis_set ==
      analysis_set,
    ,
    drop = FALSE
  ]

  boxplot(
    late_remodeling_state ~ condition,
    data = x,
    order = c(
      "Skin",
      "Wound1",
      "Wound7",
      "Wound30"
    ),
    xlab = "",
    ylab =
      "Late-remodeling state",
    main =
      analysis_set,
    las = 2
  )

  stripchart(
    late_remodeling_state ~ condition,
    data = x,
    vertical = TRUE,
    method = "jitter",
    add = TRUE,
    pch = 16
  )
}

par(
  op
)

dev.off()

# 供者内连线
pdf(
  file.path(
    FIG_DIR,
    "Figure_spatial_donor_trajectories.pdf"
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
    5,
    4,
    3,
    1
  )
)

for (
  analysis_set in c(
    "high_specificity_spots",
    "broad_top_quartile_spots"
  )
) {
  x <- section_scores[
    section_scores$analysis_set ==
      analysis_set,
    ,
    drop = FALSE
  ]

  plot(
    NA,
    xlim = c(
      0,
      3
    ),
    ylim = range(
      x$ordinal_wound_state_position,
      finite = TRUE
    ),
    xaxt = "n",
    xlab = "",
    ylab =
      "Ordinal wound-state position",
    main = paste0(
      analysis_set,
      "\nordinal trajectory"
    )
  )

  axis(
    1,
    at = 0:3,
    labels = c(
      "Skin",
      "D1",
      "D7",
      "D30"
    )
  )

  for (
    donor in unique(
      x$donor
    )
  ) {
    xd <- x[
      x$donor ==
        donor,
      ,
      drop = FALSE
    ]

    xd <- xd[
      order(
        xd$stage_code
      ),
      ,
      drop = FALSE
    ]

    lines(
      xd$stage_code,
      xd$ordinal_wound_state_position,
      type = "b",
      pch = 16
    )

    text(
      xd$stage_code[4],
      xd$ordinal_wound_state_position[4],
      labels = donor,
      pos = 4,
      cex = 0.7
    )
  }
}

for (
  analysis_set in c(
    "high_specificity_spots",
    "broad_top_quartile_spots"
  )
) {
  x <- section_scores[
    section_scores$analysis_set ==
      analysis_set,
    ,
    drop = FALSE
  ]

  plot(
    NA,
    xlim = c(
      0,
      3
    ),
    ylim = range(
      x$late_remodeling_state,
      finite = TRUE
    ),
    xaxt = "n",
    xlab = "",
    ylab =
      "Late-remodeling state",
    main = paste0(
      analysis_set,
      "\nlate-remodeling trajectory"
    )
  )

  axis(
    1,
    at = 0:3,
    labels = c(
      "Skin",
      "D1",
      "D7",
      "D30"
    )
  )

  for (
    donor in unique(
      x$donor
    )
  ) {
    xd <- x[
      x$donor ==
        donor,
      ,
      drop = FALSE
    ]

    xd <- xd[
      order(
        xd$stage_code
      ),
      ,
      drop = FALSE
    ]

    lines(
      xd$stage_code,
      xd$late_remodeling_state,
      type = "b",
      pch = 16
    )

    text(
      xd$stage_code[4],
      xd$late_remodeling_state[4],
      labels = donor,
      pos = 4,
      cex = 0.7
    )
  }
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
    position_file_summary =
      position_file_summary,
    spot_selection_summary =
      spot_selection_summary,
    spot_scores =
      spot_scores,
    high_specificity =
      high_result,
    broad_top_quartile =
      broad_result,
    section_scores =
      section_scores,
    overall_metrics =
      overall_metrics,
    donor_metrics =
      donor_metrics,
    stage_summary =
      stage_summary,
    decision_criteria =
      decision_criteria,
    fixed_signatures =
      signatures,
    frozen_marker_sets =
      marker_sets
  ),
  file = file.path(
    OBJECT_DIR,
    "GSE241124_spatial_fixed_signature_validation.rds"
  ),
  compress = "gzip"
)

# ----------------------------- 最终结论 ----------------------------------------
decision_lines <- c(
  "第八阶段：GSE241124正常伤口空间固定签名验证结论（V1.1修正版）",
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
  "高特异性空间集合：",
  paste0(
    "- exact accuracy = ",
    round(
      high_overall$exact_accuracy,
      3
    )
  ),
  paste0(
    "- adjacent accuracy = ",
    round(
      high_overall$adjacent_accuracy,
      3
    )
  ),
  paste0(
    "- mean ordinal error = ",
    round(
      high_overall$mean_ordinal_error,
      3
    )
  ),
  paste0(
    "- pooled Spearman = ",
    round(
      high_overall$pooled_spearman,
      3
    )
  ),
  paste0(
    "- donors with D30 highest late-remodeling = ",
    sum(
      high_donor$wound30_has_highest_late_remodeling
    ),
    "/4"
  ),
  paste0(
    "- donors with positive ordinal Spearman = ",
    sum(
      high_donor$spearman >
        0
    ),
    "/4"
  ),
  "",
  "宽松top-quartile空间集合：",
  paste0(
    "- exact accuracy = ",
    round(
      broad_overall$exact_accuracy,
      3
    )
  ),
  paste0(
    "- adjacent accuracy = ",
    round(
      broad_overall$adjacent_accuracy,
      3
    )
  ),
  paste0(
    "- pooled Spearman = ",
    round(
      broad_overall$pooled_spearman,
      3
    )
  ),
  "",
  "解释边界：",
  "- GSE241124与GSE241132来自不同供者，但属于同一正常伤口研究框架。",
  "- 空间切片为独立模态验证，不是新的疾病队列。",
  "- Visium spots为混合细胞，本文只能称fibroblast-enriched spatial domains。",
  "- high-specificity spots为主要空间集合；broad top-quartile为敏感性集合。",
  "- spot级模块用于空间分布展示，正式阶段映射基于切片级pseudobulk。",
  "- 未在空间数据中重新选择伤口阶段基因。",
  "- ordinal position不是精确术后日龄。"
)

safe_write_lines(
  decision_lines,
  file.path(
    REPORT_DIR,
    "14_STAGE8_DECISION.txt"
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
      "STAGE8_COMPLETED=",
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
      "HIGH_SPOT_COUNT_PASS=",
      high_spot_count_pass
    ),
    paste0(
      "BROAD_SPOT_COUNT_PASS=",
      broad_spot_count_pass
    ),
    paste0(
      "HIGH_STRONG_VALIDATION=",
      high_strong
    ),
    paste0(
      "BROAD_CONCORDANT_VALIDATION=",
      broad_concordant
    ),
    paste0(
      "HIGH_MODERATE_VALIDATION=",
      high_moderate
    )
  ),
  file.path(
    REPORT_DIR,
    "STAGE8_COMPLETED.txt"
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
  "第八阶段完成"
)

cat(
  "\n============================================================\n"
)

cat(
  "第八阶段运行完成。\n"
)

cat(
  "本地分析对象保存在：",
  file.path(
    OBJECT_DIR,
    "GSE241124_spatial_fixed_signature_validation.rds"
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
