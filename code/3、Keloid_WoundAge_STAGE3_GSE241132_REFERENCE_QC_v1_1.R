# =============================================================================
# 项目：瘢痕疙瘩“伤口状态持续与修复终止失败”
# 第三阶段：正常人体伤口单细胞参考队列（GSE241132）正式QC与可建模性审计（V1.1修正版）
#
# 根目录：
#   C:/Users/33652/Desktop/LZZ/1、yssj
#
# 本阶段只处理GSE241132，原因：
#   正常伤口时间参照能否可靠建立，是整个课题最关键的停止/继续闸门。
#
# 本阶段完成：
#   1. 读取12个单细胞矩阵；
#   2. 将矩阵barcode与作者细胞元数据逐样本匹配；
#   3. 计算nCount、nFeature、线粒体/核糖体/血红蛋白比例；
#   4. 采用保守QC，不重新把所有细胞混合聚类；
#   5. 按“供者×阶段×主细胞类型”生成pseudobulk计数；
#   6. 判断成纤维细胞及其他主要细胞是否具备时间建模条件；
#   7. 自动生成第三阶段检查包供上传。
#
# 注意：
#   - 不进行正式差异表达或时间模型；
#   - 不随机拆分细胞；
#   - 不删除或覆盖原始数据；
#   - 本阶段生成的RDS对象保留在本机，不放入上传ZIP。
#   - V1.1修复table()长度不一致及Matrix 1.7+稀疏矩阵转换警告，并增加空集保护。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANG = "en")

ROOT_DIR <- "C:/Users/33652/Desktop/LZZ/1、yssj"
STAGE_DIR <- file.path(ROOT_DIR, "03_STAGE3_REFERENCE_QC")
REPORT_DIR <- file.path(STAGE_DIR, "report")
FIG_DIR <- file.path(STAGE_DIR, "figures")
OBJECT_DIR <- file.path(STAGE_DIR, "objects")
WORK_DIR <- file.path(STAGE_DIR, "_temporary_work")

PACKAGE_ZIP <- file.path(ROOT_DIR, "第三阶段_GSE241132正式QC与可建模性检查包.zip")
PACKAGE_TARGZ <- file.path(ROOT_DIR, "第三阶段_GSE241132正式QC与可建模性检查包.tar.gz")

if (dir.exists(STAGE_DIR)) unlink(STAGE_DIR, recursive = TRUE, force = TRUE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(REPORT_DIR, "00_stage3_log.txt")
WARNING_FILE <- file.path(REPORT_DIR, "00_stage3_warnings.txt")

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                " | ", paste(..., collapse = " "))
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
    utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8"),
    error = function(e) warn_msg("CSV写入失败:", basename(path), "|", conditionMessage(e))
  )
}

safe_write_lines <- function(x, path) {
  tryCatch(
    writeLines(x, path, useBytes = TRUE),
    error = function(e) warn_msg("文本写入失败:", basename(path), "|", conditionMessage(e))
  )
}

normalize_slash <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

find_file_exact <- function(fname) {
  x <- list.files(ROOT_DIR, recursive = TRUE, full.names = TRUE,
                  all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  x <- x[
    !grepl("/03_STAGE3_REFERENCE_QC/", normalize_slash(x), fixed = TRUE) &
    basename(x) != basename(PACKAGE_ZIP) &
    basename(x) != basename(PACKAGE_TARGZ)
  ]
  hit <- which(tolower(basename(x)) == tolower(fname))
  if (!length(hit)) return(NA_character_)
  x[hit[1]]
}

# ----------------------------- 依赖包 ------------------------------------------
if (!requireNamespace("Matrix", quietly = TRUE)) {
  install.packages("Matrix", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("无法安装或加载Matrix包。")
}

# ----------------------------- 基础函数 ----------------------------------------
read_text_table <- function(path, sep = "\t") {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  utils::read.table(
    con,
    header = TRUE,
    sep = sep,
    quote = "\"",
    comment.char = "",
    fill = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

read_noheader_table <- function(path, sep = "\t") {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
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

read_mtx_gz <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  x <- Matrix::readMM(con)

  # Matrix 1.7+不再推荐直接由dgTMatrix转换为dgCMatrix。
  # 按官方提示转换为CsparseMatrix，实际返回可用于后续计算的列压缩稀疏矩阵。
  if (!inherits(x, "CsparseMatrix")) {
    x <- methods::as(x, "CsparseMatrix")
  }
  x
}

extract_10x_barcode <- function(x) {
  x <- as.character(x)
  out <- rep(NA_character_, length(x))
  hit <- grepl("([ACGTN]{12,}-[0-9]+)$", x, ignore.case = TRUE, perl = TRUE)
  out[hit] <- toupper(sub("^.*?([ACGTN]{12,}-[0-9]+)$", "\\1",
                          x[hit], ignore.case = TRUE, perl = TRUE))
  out
}

select_metadata_barcode <- function(md) {
  candidates <- list()

  rn <- rownames(md)
  if (!is.null(rn) && length(rn) == nrow(md) &&
      !identical(rn, as.character(seq_len(nrow(md))))) {
    candidates[[".rownames"]] <- rn
  }

  for (nm in names(md)) {
    v <- md[[nm]]
    if (is.character(v) || is.factor(v)) {
      candidates[[nm]] <- as.character(v)
    }
  }

  if (!length(candidates)) {
    return(list(name = NA_character_, values = rep(NA_character_, nrow(md)), rate = 0))
  }

  rates <- vapply(candidates, function(v) {
    mean(!is.na(extract_10x_barcode(v)))
  }, numeric(1))

  best <- which.max(rates)
  list(
    name = names(candidates)[best],
    values = candidates[[best]],
    rate = rates[best]
  )
}

extract_nested_zip_samples <- function(outer_tar, out_dir) {
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE, force = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  utils::untar(outer_tar, exdir = out_dir)

  zips <- list.files(
    out_dir,
    recursive = TRUE,
    full.names = TRUE,
    pattern = "\\.zip$",
    ignore.case = TRUE
  )

  if (length(zips) != 12L) {
    warn_msg("GSE241132嵌套ZIP数量不是12，实际为:", length(zips))
  }

  zips
}

sample_id_from_zip <- function(path) {
  x <- basename(path)
  x <- sub("\\.zip$", "", x, ignore.case = TRUE)
  x <- sub("^GSM[0-9]+_", "", x)
  x
}

extract_triplet_from_zip <- function(zip_path, sample_tmp) {
  if (dir.exists(sample_tmp)) unlink(sample_tmp, recursive = TRUE, force = TRUE)
  dir.create(sample_tmp, recursive = TRUE, showWarnings = FALSE)

  listing <- utils::unzip(zip_path, list = TRUE)
  members <- as.character(listing$Name)

  m <- members[grepl("matrix\\.mtx(\\.gz)?$", members, ignore.case = TRUE)]
  b <- members[grepl("barcodes?\\.tsv(\\.gz)?$", members, ignore.case = TRUE)]
  f <- members[grepl("(features|genes)\\.tsv(\\.gz)?$", members, ignore.case = TRUE)]

  if (!length(m) || !length(b) || !length(f)) {
    stop("嵌套ZIP缺少matrix/barcodes/features：", basename(zip_path))
  }

  utils::unzip(
    zip_path,
    files = c(m[1], b[1], f[1]),
    exdir = sample_tmp,
    junkpaths = FALSE
  )

  files <- list.files(sample_tmp, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, include.dirs = FALSE, no.. = TRUE)

  locate <- function(regex) {
    hit <- files[grepl(regex, basename(files), ignore.case = TRUE, perl = TRUE)]
    if (!length(hit)) return(NA_character_)
    hit[1]
  }

  list(
    matrix = locate("matrix\\.mtx(\\.gz)?$"),
    barcodes = locate("barcodes?\\.tsv(\\.gz)?$"),
    features = locate("(features|genes)\\.tsv(\\.gz)?$")
  )
}

safe_quantile <- function(x, probs) {
  if (!length(x) || all(is.na(x))) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE))
}

first_unique_nonempty <- function(x, label, sample_id) {
  x <- trimws(as.character(x))
  x <- unique(x[!is.na(x) & nzchar(x)])
  if (!length(x)) {
    warn_msg(sample_id, label, "为空，写入NA")
    return(NA_character_)
  }
  if (length(x) > 1L) {
    warn_msg(sample_id, label, "存在多个取值：", paste(x, collapse = " | "),
             "；使用第一个值")
  }
  x[1]
}

safe_rbind <- function(x, template = NULL) {
  if (!length(x)) {
    if (is.null(template)) return(data.frame())
    return(template[0, , drop = FALSE])
  }
  do.call(rbind, x)
}

# ----------------------------- 输入检查 ----------------------------------------
log_msg("第三阶段开始")

raw_tar <- find_file_exact("GSE241132_RAW.tar")
meta_file <- find_file_exact("GSE241132_cell_metadata.txt.gz")

if (is.na(raw_tar)) stop("未找到GSE241132_RAW.tar")
if (is.na(meta_file)) stop("未找到GSE241132_cell_metadata.txt.gz")

log_msg("读取作者细胞元数据:", meta_file)
metadata <- read_text_table(meta_file, sep = "\t")

required_cols <- c(
  "orig.ident", "Patient", "Condition", "Doublet",
  "newMainCellTypes", "newCellTypes"
)
missing_cols <- setdiff(required_cols, names(metadata))
if (length(missing_cols)) {
  stop("GSE241132元数据缺少关键字段：", paste(missing_cols, collapse = ", "))
}

bc_info <- select_metadata_barcode(metadata)
metadata$barcode_source <- bc_info$name
metadata$barcode_norm <- extract_10x_barcode(bc_info$values)

if (bc_info$rate < 0.90) {
  warn_msg("元数据barcode自动识别率低于90%，来源字段:",
           bc_info$name, "，识别率:", round(bc_info$rate, 4))
}

metadata$orig.ident <- trimws(as.character(metadata$orig.ident))
metadata$Patient <- trimws(as.character(metadata$Patient))
metadata$Condition <- trimws(as.character(metadata$Condition))
metadata$Doublet <- trimws(as.character(metadata$Doublet))
metadata$newMainCellTypes <- trimws(as.character(metadata$newMainCellTypes))
metadata$newCellTypes <- trimws(as.character(metadata$newCellTypes))

# ----------------------------- 逐样本处理 --------------------------------------
outer_dir <- file.path(WORK_DIR, "outer")
sample_zips <- extract_nested_zip_samples(raw_tar, outer_dir)

sample_order <- c(
  "PWH26D0", "PWH26D1", "PWH26D7", "PWH26D30",
  "PWH27D0", "PWH27D1", "PWH27D7", "PWH27D30",
  "PWH28D0", "PWH28D1", "PWH28D7", "PWH28D30"
)

zip_ids <- vapply(sample_zips, sample_id_from_zip, character(1))
sample_zips <- sample_zips[match(sample_order, zip_ids)]
if (any(is.na(sample_zips))) {
  stop("无法按预期找到全部12个样本ZIP。")
}

qc_sample_rows <- list()
barcode_match_rows <- list()
celltype_rows <- list()
feasibility_rows <- list()
pseudobulk_list <- list()
gene_reference <- NULL
gene_reference_id <- NULL
feature_consistency <- list()

pdf(file.path(FIG_DIR, "Figure_QC_distributions_by_sample.pdf"),
    width = 10, height = 8, onefile = TRUE)

for (zip_path in sample_zips) {
  sid <- sample_id_from_zip(zip_path)
  log_msg("处理样本:", sid)

  sample_tmp <- file.path(WORK_DIR, paste0("sample_", sid))
  triplet <- extract_triplet_from_zip(zip_path, sample_tmp)

  barcodes_df <- read_noheader_table(triplet$barcodes, sep = "\t")
  features_df <- read_noheader_table(triplet$features, sep = "\t")
  counts <- read_mtx_gz(triplet$matrix)

  raw_barcodes <- toupper(as.character(barcodes_df[[1]]))
  gene_id <- as.character(features_df[[1]])
  gene_symbol <- if (ncol(features_df) >= 2L) as.character(features_df[[2]]) else gene_id

  if (nrow(counts) != length(gene_symbol) || ncol(counts) != length(raw_barcodes)) {
    stop("矩阵维度与features/barcodes不一致：", sid)
  }

  # 先以稳定的gene_id核对并对齐各样本，避免按行号错误合并pseudobulk。
  if (is.null(gene_reference)) {
    gene_reference <- gene_symbol
    gene_reference_id <- gene_id
    feature_id_identical <- TRUE
    feature_symbol_identical <- TRUE
  } else {
    feature_id_identical <- identical(gene_id, gene_reference_id)
    feature_symbol_identical <- identical(gene_symbol, gene_reference)

    if (!feature_id_identical) {
      if (anyDuplicated(gene_id) || anyDuplicated(gene_reference_id)) {
        stop("样本", sid, "的gene_id顺序不同且存在重复gene_id，无法安全对齐。")
      }
      reorder_idx <- match(gene_reference_id, gene_id)
      if (anyNA(reorder_idx)) {
        stop("样本", sid, "缺少参考样本中的部分gene_id，无法安全合并pseudobulk。")
      }
      counts <- counts[reorder_idx, , drop = FALSE]
      gene_id <- gene_id[reorder_idx]
      gene_symbol <- gene_symbol[reorder_idx]
      feature_id_identical <- identical(gene_id, gene_reference_id)
      feature_symbol_identical <- identical(gene_symbol, gene_reference)
    }
  }

  rownames(counts) <- make.unique(gene_reference)
  colnames(counts) <- raw_barcodes

  feature_consistency[[length(feature_consistency) + 1L]] <- data.frame(
    sample_id = sid,
    n_features = length(gene_symbol),
    identical_gene_ids_to_first = feature_id_identical,
    identical_gene_symbols_to_first = feature_symbol_identical,
    stringsAsFactors = FALSE
  )

  md_s <- metadata[metadata$orig.ident == sid, , drop = FALSE]
  if (!nrow(md_s)) stop("元数据中没有样本：", sid)

  md_s <- md_s[!is.na(md_s$barcode_norm), , drop = FALSE]
  md_s <- md_s[!duplicated(md_s$barcode_norm), , drop = FALSE]

  idx <- match(md_s$barcode_norm, raw_barcodes)
  matched <- !is.na(idx)
  match_rate <- if (length(matched)) mean(matched) else 0

  barcode_match_rows[[length(barcode_match_rows) + 1L]] <- data.frame(
    sample_id = sid,
    n_matrix_barcodes = length(raw_barcodes),
    n_metadata_rows = nrow(md_s),
    n_metadata_matched = sum(matched),
    metadata_match_rate = match_rate,
    n_matrix_not_in_metadata = length(raw_barcodes) - length(unique(idx[matched])),
    barcode_source = bc_info$name,
    stringsAsFactors = FALSE
  )

  if (match_rate < 0.90) {
    warn_msg(sid, "元数据与矩阵barcode匹配率低于90%:",
             round(match_rate, 4))
  }
  if (!any(matched)) {
    stop("样本", sid, "没有任何metadata barcode与表达矩阵匹配。")
  }

  md_s <- md_s[matched, , drop = FALSE]
  counts <- counts[, idx[matched], drop = FALSE]
  colnames(counts) <- md_s$barcode_norm

  n_count <- Matrix::colSums(counts)
  n_feature <- Matrix::colSums(counts > 0)

  gene_upper <- toupper(rownames(counts))
  mt_idx <- grepl("^MT-", gene_upper)
  ribo_idx <- grepl("^RP[SL][0-9]", gene_upper)
  hb_idx <- grepl("^HB[ABDEGMQZ][0-9A-Z]*$", gene_upper)

  pct_mt <- if (any(mt_idx)) Matrix::colSums(counts[mt_idx, , drop = FALSE]) /
    pmax(n_count, 1) * 100 else rep(0, ncol(counts))

  pct_ribo <- if (any(ribo_idx)) Matrix::colSums(counts[ribo_idx, , drop = FALSE]) /
    pmax(n_count, 1) * 100 else rep(0, ncol(counts))

  pct_hb <- if (any(hb_idx)) Matrix::colSums(counts[hb_idx, , drop = FALSE]) /
    pmax(n_count, 1) * 100 else rep(0, ncol(counts))

  md_s$nCount_calc <- as.numeric(n_count)
  md_s$nFeature_calc <- as.numeric(n_feature)
  md_s$percent_mt_calc <- as.numeric(pct_mt)
  md_s$percent_ribo_calc <- as.numeric(pct_ribo)
  md_s$percent_hb_calc <- as.numeric(pct_hb)

  # 作者已经提供Doublet字段，本阶段采用保守QC，不再用细胞级随机算法制造额外删失。
  md_s$qc_keep <- (
    tolower(trimws(md_s$Doublet)) == "singlet" &
      is.finite(md_s$nFeature_calc) &
      is.finite(md_s$nCount_calc) &
      is.finite(md_s$percent_mt_calc) &
      md_s$nFeature_calc >= 200 &
      md_s$nFeature_calc <= 7500 &
      md_s$nCount_calc >= 500 &
      md_s$percent_mt_calc <= 20
  )
  md_s$qc_keep[is.na(md_s$qc_keep)] <- FALSE

  donor_s <- first_unique_nonempty(md_s$Patient, "Patient", sid)
  condition_s <- first_unique_nonempty(md_s$Condition, "Condition", sid)

  q_nf <- safe_quantile(md_s$nFeature_calc, c(0.01, 0.25, 0.50, 0.75, 0.99))
  q_nc <- safe_quantile(md_s$nCount_calc, c(0.01, 0.25, 0.50, 0.75, 0.99))
  q_mt <- safe_quantile(md_s$percent_mt_calc, c(0.01, 0.25, 0.50, 0.75, 0.99))

  qc_sample_rows[[length(qc_sample_rows) + 1L]] <- data.frame(
    sample_id = sid,
    donor = donor_s,
    condition = condition_s,
    n_matrix_barcodes = length(raw_barcodes),
    n_author_metadata_cells = nrow(md_s),
    n_qc_keep = sum(md_s$qc_keep),
    qc_retention_rate = mean(md_s$qc_keep),
    nFeature_p01 = q_nf[1],
    nFeature_p25 = q_nf[2],
    nFeature_median = q_nf[3],
    nFeature_p75 = q_nf[4],
    nFeature_p99 = q_nf[5],
    nCount_p01 = q_nc[1],
    nCount_p25 = q_nc[2],
    nCount_median = q_nc[3],
    nCount_p75 = q_nc[4],
    nCount_p99 = q_nc[5],
    pctMT_p01 = q_mt[1],
    pctMT_p25 = q_mt[2],
    pctMT_median = q_mt[3],
    pctMT_p75 = q_mt[4],
    pctMT_p99 = q_mt[5],
    stringsAsFactors = FALSE
  )

  # 逐样本QC图
  op <- par(no.readonly = TRUE)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  hist(md_s$nFeature_calc, breaks = 50,
       main = paste0(sid, " nFeature"), xlab = "Detected genes")
  abline(v = c(200, 7500), lty = 2)
  hist(log10(md_s$nCount_calc + 1), breaks = 50,
       main = paste0(sid, " log10 nCount"), xlab = "log10 counts + 1")
  abline(v = log10(501), lty = 2)
  hist(md_s$percent_mt_calc, breaks = 50,
       main = paste0(sid, " mitochondrial %"), xlab = "percent.mt")
  abline(v = 20, lty = 2)
  plot(log10(md_s$nCount_calc + 1), md_s$nFeature_calc,
       pch = 16, cex = 0.35,
       xlab = "log10 counts + 1", ylab = "Detected genes",
       main = paste0(sid, " QC scatter"))
  par(op)

  keep_idx <- which(md_s$qc_keep)
  counts_qc <- counts[, keep_idx, drop = FALSE]
  md_qc <- md_s[keep_idx, , drop = FALSE]

  # 每个样本的本地对象，不进入上传检查包
  saveRDS(
    list(counts = counts_qc, metadata = md_qc),
    file = file.path(OBJECT_DIR, paste0("GSE241132_", sid, "_QC.rds")),
    compress = "gzip"
  )

  if (!nrow(md_qc)) {
    warn_msg(sid, "QC后没有保留细胞；该样本不会生成细胞类型计数或pseudobulk。")
  } else {
    # 修复V1.0错误：table()的每个参数必须等长，不能把标量sid与细胞向量直接混用。
    ct_input <- data.frame(
      sample_id = rep(sid, nrow(md_qc)),
      donor = rep(donor_s, nrow(md_qc)),
      condition = rep(condition_s, nrow(md_qc)),
      main_cell_type = trimws(as.character(md_qc$newMainCellTypes)),
      stringsAsFactors = FALSE
    )
    ct_input <- ct_input[
      !is.na(ct_input$main_cell_type) & nzchar(ct_input$main_cell_type),
      ,
      drop = FALSE
    ]

    if (nrow(ct_input)) {
      ct_df <- stats::aggregate(
        x = list(n_cells = rep(1L, nrow(ct_input))),
        by = ct_input,
        FUN = sum
      )
      celltype_rows[[length(celltype_rows) + 1L]] <- ct_df
    } else {
      warn_msg(sid, "QC后细胞缺少有效newMainCellTypes标签。")
    }

    # pseudobulk：每个样本×主细胞类型
    cell_types <- sort(unique(trimws(as.character(md_qc$newMainCellTypes))))
    cell_types <- cell_types[!is.na(cell_types) & nzchar(cell_types)]

    for (cell_type in cell_types) {
      j <- which(trimws(as.character(md_qc$newMainCellTypes)) == cell_type)
      if (length(j) < 20L) next

      pb <- Matrix::rowSums(counts_qc[, j, drop = FALSE])
      key <- paste(sid, cell_type, sep = "||")
      names(pb) <- make.unique(gene_reference)
      pseudobulk_list[[key]] <- pb

      feasibility_rows[[length(feasibility_rows) + 1L]] <- data.frame(
        sample_id = sid,
        donor = donor_s,
        condition = condition_s,
        main_cell_type = cell_type,
        n_cells = length(j),
        pseudobulk_library_size = sum(pb),
        pseudobulk_detected_genes = sum(pb > 0),
        stringsAsFactors = FALSE
      )
    }
  }

  rm(counts, counts_qc, md_s, md_qc)
  gc(verbose = FALSE)
  unlink(sample_tmp, recursive = TRUE, force = TRUE)
}

dev.off()

# ----------------------------- 汇总与对象输出 ----------------------------------
qc_sample_summary <- safe_rbind(qc_sample_rows)
barcode_match_summary <- safe_rbind(barcode_match_rows)
celltype_counts <- safe_rbind(celltype_rows)
pb_sample_info <- safe_rbind(feasibility_rows)
feature_consistency_df <- safe_rbind(feature_consistency)

if (!nrow(qc_sample_summary)) stop("未生成样本QC汇总。")
if (!nrow(barcode_match_summary)) stop("未生成barcode匹配汇总。")
if (!nrow(celltype_counts)) stop("未生成任何主细胞类型计数。")
if (!nrow(pb_sample_info)) stop("未生成任何pseudobulk样本信息。")

safe_write_csv(qc_sample_summary,
               file.path(REPORT_DIR, "01_sample_QC_summary.csv"))
safe_write_csv(barcode_match_summary,
               file.path(REPORT_DIR, "02_barcode_matching_summary.csv"))
safe_write_csv(celltype_counts,
               file.path(REPORT_DIR, "03_main_celltype_counts_after_QC.csv"))
safe_write_csv(pb_sample_info,
               file.path(REPORT_DIR, "04_pseudobulk_sample_information.csv"))
safe_write_csv(feature_consistency_df,
               file.path(REPORT_DIR, "05_feature_consistency.csv"))

# 合并pseudobulk矩阵并保存在本机objects目录
if (!length(pseudobulk_list)) {
  stop("未生成任何pseudobulk。")
}

pb_lengths <- vapply(pseudobulk_list, length, integer(1))
if (length(unique(pb_lengths)) != 1L || unique(pb_lengths) != length(gene_reference)) {
  stop("pseudobulk向量长度不一致，停止合并。")
}
pb_matrix <- do.call(cbind, pseudobulk_list)
rownames(pb_matrix) <- make.unique(gene_reference)
saveRDS(
  list(
    counts = pb_matrix,
    sample_info = pb_sample_info,
    gene_symbol = gene_reference,
    gene_id = gene_reference_id
  ),
  file = file.path(OBJECT_DIR, "GSE241132_mainCellType_pseudobulk_counts.rds"),
  compress = "gzip"
)

# 细胞类型时间建模可行性
all_design <- unique(qc_sample_summary[, c("sample_id", "donor", "condition")])

ct_levels <- sort(unique(celltype_counts$main_cell_type))
model_feasibility <- do.call(rbind, lapply(ct_levels, function(ct) {
  sub <- celltype_counts[celltype_counts$main_cell_type == ct, , drop = FALSE]
  m <- merge(all_design, sub, by = c("sample_id", "donor", "condition"), all.x = TRUE)
  m$n_cells[is.na(m$n_cells)] <- 0L

  data.frame(
    main_cell_type = ct,
    n_donor_stage_groups = nrow(m),
    min_cells = min(m$n_cells),
    median_cells = stats::median(m$n_cells),
    groups_ge_30 = sum(m$n_cells >= 30),
    groups_ge_50 = sum(m$n_cells >= 50),
    groups_ge_100 = sum(m$n_cells >= 100),
    all_12_ge_30 = all(m$n_cells >= 30),
    all_12_ge_50 = all(m$n_cells >= 50),
    all_12_ge_100 = all(m$n_cells >= 100),
    stringsAsFactors = FALSE
  )
}))

safe_write_csv(model_feasibility,
               file.path(REPORT_DIR, "06_celltype_time_model_feasibility.csv"))

# 简单热图：主细胞类型在12个样本中的细胞数
heat_tab <- xtabs(n_cells ~ main_cell_type + sample_id, data = celltype_counts)
pdf(file.path(FIG_DIR, "Figure_celltype_counts_heatmap.pdf"),
    width = 12, height = max(6, nrow(heat_tab) * 0.42))
par(mar = c(9, 12, 3, 2))
image(
  x = seq_len(ncol(heat_tab)),
  y = seq_len(nrow(heat_tab)),
  z = t(log10(heat_tab + 1)),
  axes = FALSE,
  xlab = "", ylab = "",
  main = "GSE241132 main cell-type counts after QC (log10[n+1])"
)
axis(1, at = seq_len(ncol(heat_tab)), labels = colnames(heat_tab),
     las = 2, cex.axis = 0.8)
axis(2, at = seq_len(nrow(heat_tab)), labels = rownames(heat_tab),
     las = 2, cex.axis = 0.8)
box()
dev.off()

# ----------------------------- 闸门判断 ----------------------------------------
barcode_pass <- all(barcode_match_summary$metadata_match_rate >= 0.90)
qc_retention_pass <- all(qc_sample_summary$qc_retention_rate >= 0.70)

fib_row <- model_feasibility[
  grepl("^fibro", tolower(model_feasibility$main_cell_type)),
  ,
  drop = FALSE
]

fibroblast_pass <- nrow(fib_row) >= 1L &&
  any(fib_row$all_12_ge_100, na.rm = TRUE)

n_celltypes_all12_ge50 <- sum(model_feasibility$all_12_ge_50, na.rm = TRUE)
reference_ready <- barcode_pass &&
  qc_retention_pass &&
  fibroblast_pass &&
  n_celltypes_all12_ge50 >= 2L

decision_lines <- c(
  "第三阶段：GSE241132正常伤口参考队列正式QC结论（V1.1修正版）",
  "============================================================",
  paste0("运行时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("12个样本是否全部完成：", nrow(qc_sample_summary) == 12L),
  paste0("所有样本barcode匹配率≥90%：", barcode_pass),
  paste0("所有样本QC保留率≥70%：", qc_retention_pass),
  paste0("成纤维细胞是否在全部12个供者-阶段组≥100个细胞：", fibroblast_pass),
  paste0("全部12组均≥50个细胞的主细胞类型数量：", n_celltypes_all12_ge50),
  paste0("是否允许进入第四阶段时间模型：", reference_ready),
  "",
  "解释：",
  "- reference_ready=TRUE：可进入供者留一的伤口阶段模型。",
  "- reference_ready=FALSE但成纤维细胞通过：可缩小为成纤维细胞主导模型。",
  "- 成纤维细胞不通过：停止“分子伤口年龄”主线，改为更保守的阶段签名研究。",
  "",
  "本阶段未进行：",
  "- 未训练时间分类器；",
  "- 未进行差异表达；",
  "- 未把细胞随机分割为训练/测试集；",
  "- 未分析瘢痕疙瘩队列。"
)

safe_write_lines(decision_lines,
                 file.path(REPORT_DIR, "07_STAGE3_DECISION.txt"))

capture.output(sessionInfo(),
               file = file.path(REPORT_DIR, "08_SESSION_INFO.txt"))

safe_write_lines(
  c(
    paste0("STAGE3_COMPLETED=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    paste0("BARCODE_PASS=", barcode_pass),
    paste0("QC_RETENTION_PASS=", qc_retention_pass),
    paste0("FIBROBLAST_PASS=", fibroblast_pass),
    paste0("REFERENCE_READY=", reference_ready)
  ),
  file.path(REPORT_DIR, "STAGE3_COMPLETED.txt")
)

# 临时文件删除；本地objects保留
if (dir.exists(WORK_DIR)) unlink(WORK_DIR, recursive = TRUE, force = TRUE)

# ----------------------------- 自动打包 ----------------------------------------
package_dir <- file.path(STAGE_DIR, "_package_for_review")
if (dir.exists(package_dir)) unlink(package_dir, recursive = TRUE, force = TRUE)
dir.create(package_dir, recursive = TRUE, showWarnings = FALSE)

file.copy(REPORT_DIR, package_dir, recursive = TRUE)
file.copy(FIG_DIR, package_dir, recursive = TRUE)

create_zip_windows <- function(source_dir, zip_path) {
  if (file.exists(zip_path)) file.remove(zip_path)

  esc <- function(x) gsub("'", "''", normalize_slash(x), fixed = TRUE)
  cmd <- paste0(
    "$ErrorActionPreference='Stop'; ",
    "$items=Get-ChildItem -LiteralPath '", esc(source_dir), "'; ",
    "Compress-Archive -Path $items.FullName -DestinationPath '",
    esc(zip_path), "' -Force"
  )

  try(
    system2(
      "powershell.exe",
      c("-NoProfile", "-ExecutionPolicy", "Bypass", "-Command", cmd),
      stdout = TRUE, stderr = TRUE
    ),
    silent = TRUE
  )

  file.exists(zip_path) && file.info(zip_path)$size > 0
}

zip_ok <- create_zip_windows(package_dir, PACKAGE_ZIP)

if (!zip_ok) {
  if (file.exists(PACKAGE_TARGZ)) file.remove(PACKAGE_TARGZ)
  old <- getwd()
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
  setwd(old)
}

unlink(package_dir, recursive = TRUE, force = TRUE)

log_msg("第三阶段完成")

cat("\n============================================================\n")
cat("第三阶段运行完成。\n")
cat("本地分析对象保存在：", OBJECT_DIR, "\n", sep = "")
if (file.exists(PACKAGE_ZIP)) {
  cat("请上传：", PACKAGE_ZIP, "\n", sep = "")
} else if (file.exists(PACKAGE_TARGZ)) {
  cat("请上传备用包：", PACKAGE_TARGZ, "\n", sep = "")
} else {
  cat("自动压缩失败，请手动压缩report和figures文件夹。\n")
}
cat("============================================================\n")
