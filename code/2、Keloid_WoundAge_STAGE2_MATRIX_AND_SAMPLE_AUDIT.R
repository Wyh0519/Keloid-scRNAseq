# =============================================================================
# 项目：瘢痕疙瘩“伤口状态持续与修复终止失败”
# 第二阶段：矩阵结构、样本对应、细胞/spot数量与可分析性深度审计
#
# 根目录：
#   C:/Users/33652/Desktop/LZZ/1、yssj
#
# 运行前建议补充：
#   GSE181316_RAW.tar
#   GSE181316_family.soft.gz
#   （GSE181316_series_matrix.txt.gz 可选）
#
# 目的：
#   1. 顺序解包各数据集，不把全部大型数据同时展开；
#   2. 检查嵌套ZIP/TAR.GZ内的10x矩阵；
#   3. 读取Matrix Market头部，核对基因数、细胞数、非零元素；
#   4. 核对barcodes/features行数是否与矩阵维度一致；
#   5. 汇总现有细胞元数据和空间spot元数据；
#   6. 确认GSE220300真实样本设计；
#   7. 确认GSE181297空间文件是否具有坐标/比例因子；
#   8. 最终确认GSE274709能否用于表达与复发分析；
#   9. 自动生成第二阶段检查包供上传。
#
# 特点：
#   - 只使用base R，不安装第三方包；
#   - 不把完整表达矩阵载入内存；
#   - 每个数据集审计完成后删除临时解包文件；
#   - 不删除、移动或覆盖原始下载数据。
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANG = "en")

ROOT_DIR <- "C:/Users/33652/Desktop/LZZ/1、yssj"
STAGE_DIR <- file.path(ROOT_DIR, "02_STAGE2_MATRIX_AUDIT")
REPORT_DIR <- file.path(STAGE_DIR, "report")
WORK_DIR <- file.path(STAGE_DIR, "_temporary_work")
PACKAGE_ZIP <- file.path(ROOT_DIR, "第二阶段_矩阵结构与样本映射检查包.zip")
PACKAGE_TARGZ <- file.path(ROOT_DIR, "第二阶段_矩阵结构与样本映射检查包.tar.gz")

# 仅清理本阶段旧输出，不接触原始数据
if (dir.exists(STAGE_DIR)) unlink(STAGE_DIR, recursive = TRUE, force = TRUE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(WORK_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(REPORT_DIR, "00_stage2_log.txt")
WARNING_FILE <- file.path(REPORT_DIR, "00_stage2_warnings.txt")

log_msg <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
                paste(..., collapse = " "))
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
    error = function(e) warn_msg("写CSV失败:", basename(path), "|", conditionMessage(e))
  )
}

safe_write_lines <- function(x, path) {
  tryCatch(
    writeLines(x, path, useBytes = TRUE),
    error = function(e) warn_msg("写文本失败:", basename(path), "|", conditionMessage(e))
  )
}

normalize_slash <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

extract_gse <- function(x) {
  b <- basename(x)
  out <- rep(NA_character_, length(b))
  hit <- grepl("GSE[0-9]+", b, ignore.case = TRUE, perl = TRUE)
  if (any(hit)) {
    out[hit] <- toupper(sub("^.*?(GSE[0-9]+).*$", "\\1", b[hit],
                           ignore.case = TRUE, perl = TRUE))
  }
  out
}

find_file_exact <- function(fname) {
  x <- list.files(ROOT_DIR, recursive = TRUE, full.names = TRUE,
                  all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  # 排除本阶段输出
  x <- x[!grepl("/02_STAGE2_MATRIX_AUDIT/", normalize_slash(x), fixed = TRUE)]
  hit <- which(tolower(basename(x)) == tolower(fname))
  if (!length(hit)) return(NA_character_)
  x[hit[1]]
}

find_file_pattern <- function(pattern) {
  x <- list.files(ROOT_DIR, recursive = TRUE, full.names = TRUE,
                  all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  x <- x[!grepl("/02_STAGE2_MATRIX_AUDIT/", normalize_slash(x), fixed = TRUE)]
  x[grepl(pattern, basename(x), ignore.case = TRUE, perl = TRUE)]
}

count_lines_stream <- function(path, chunk = 100000L) {
  con <- NULL
  total <- 0L
  tryCatch({
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
    repeat {
      x <- readLines(con, n = chunk, warn = FALSE)
      if (!length(x)) break
      total <- total + length(x)
    }
    total
  }, error = function(e) {
    warn_msg("行数统计失败:", basename(path), "|", conditionMessage(e))
    NA_integer_
  }, finally = {
    if (!is.null(con)) try(close(con), silent = TRUE)
  })
}

read_mtx_header <- function(path) {
  con <- NULL
  ans <- data.frame(
    n_features_matrix = NA_real_,
    n_barcodes_matrix = NA_real_,
    n_nonzero = NA_real_,
    mtx_header_ok = FALSE,
    stringsAsFactors = FALSE
  )
  tryCatch({
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
    first <- readLines(con, n = 1L, warn = FALSE)
    if (!length(first) || !grepl("^%%MatrixMarket", first)) {
      warn_msg("不是标准Matrix Market文件:", basename(path))
      return(ans)
    }
    repeat {
      line <- readLines(con, n = 1L, warn = FALSE)
      if (!length(line)) break
      if (!startsWith(line, "%")) {
        vals <- suppressWarnings(as.numeric(strsplit(trimws(line), "\\s+")[[1]]))
        if (length(vals) >= 3L && all(is.finite(vals[1:3]))) {
          ans$n_features_matrix <- vals[1]
          ans$n_barcodes_matrix <- vals[2]
          ans$n_nonzero <- vals[3]
          ans$mtx_header_ok <- TRUE
        }
        break
      }
    }
    ans
  }, error = function(e) {
    warn_msg("MTX头读取失败:", basename(path), "|", conditionMessage(e))
    ans
  }, finally = {
    if (!is.null(con)) try(close(con), silent = TRUE)
  })
}

list_archive <- function(path) {
  tryCatch({
    if (grepl("\\.zip$", path, ignore.case = TRUE)) {
      z <- utils::unzip(path, list = TRUE)
      as.character(z$Name)
    } else {
      as.character(utils::untar(path, list = TRUE))
    }
  }, error = function(e) {
    warn_msg("无法列出嵌套压缩包:", basename(path), "|", conditionMessage(e))
    character(0)
  })
}

extract_members <- function(archive, members, exdir) {
  dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    if (grepl("\\.zip$", archive, ignore.case = TRUE)) {
      utils::unzip(archive, files = members, exdir = exdir, junkpaths = FALSE)
    } else {
      utils::untar(archive, files = members, exdir = exdir)
    }
    TRUE
  }, error = function(e) {
    warn_msg("提取成员失败:", basename(archive), "|", conditionMessage(e))
    FALSE
  })
}

locate_one <- function(root, regex) {
  x <- list.files(root, recursive = TRUE, full.names = TRUE,
                  all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  hit <- x[grepl(regex, basename(x), ignore.case = TRUE, perl = TRUE)]
  if (!length(hit)) return(NA_character_)
  hit[1]
}

audit_triplet_paths <- function(gse, sample_id, matrix_path, barcode_path, feature_path, source_archive) {
  m <- if (!is.na(matrix_path) && file.exists(matrix_path)) read_mtx_header(matrix_path) else
    data.frame(n_features_matrix=NA_real_, n_barcodes_matrix=NA_real_,
               n_nonzero=NA_real_, mtx_header_ok=FALSE)

  nbc <- if (!is.na(barcode_path) && file.exists(barcode_path)) count_lines_stream(barcode_path) else NA_integer_
  nfeat <- if (!is.na(feature_path) && file.exists(feature_path)) count_lines_stream(feature_path) else NA_integer_

  valid <- isTRUE(m$mtx_header_ok[1]) &&
    !is.na(nbc) && !is.na(nfeat) &&
    m$n_barcodes_matrix[1] == nbc &&
    m$n_features_matrix[1] == nfeat

  data.frame(
    gse = gse,
    sample_id = sample_id,
    source_archive = basename(source_archive),
    matrix_file = ifelse(is.na(matrix_path), NA_character_, basename(matrix_path)),
    barcode_file = ifelse(is.na(barcode_path), NA_character_, basename(barcode_path)),
    feature_file = ifelse(is.na(feature_path), NA_character_, basename(feature_path)),
    n_features_matrix = m$n_features_matrix[1],
    n_barcodes_matrix = m$n_barcodes_matrix[1],
    n_nonzero = m$n_nonzero[1],
    n_barcode_lines = nbc,
    n_feature_lines = nfeat,
    mtx_header_ok = m$mtx_header_ok[1],
    dimensions_match = valid,
    stringsAsFactors = FALSE
  )
}

audit_nested_10x_archive <- function(gse, archive, sample_id, work_parent) {
  members <- list_archive(archive)
  if (!length(members)) {
    return(data.frame(
      gse=gse, sample_id=sample_id, source_archive=basename(archive),
      matrix_file=NA, barcode_file=NA, feature_file=NA,
      n_features_matrix=NA, n_barcodes_matrix=NA, n_nonzero=NA,
      n_barcode_lines=NA, n_feature_lines=NA,
      mtx_header_ok=FALSE, dimensions_match=FALSE
    ))
  }

  matrix_member <- members[grepl("matrix\\.mtx(\\.gz)?$", members, ignore.case = TRUE)]
  barcode_member <- members[grepl("barcodes?\\.tsv(\\.gz)?$", members, ignore.case = TRUE)]
  feature_member <- members[grepl("(features|genes)\\.tsv(\\.gz)?$", members, ignore.case = TRUE)]

  chosen <- unique(c(matrix_member[1], barcode_member[1], feature_member[1]))
  chosen <- chosen[!is.na(chosen) & nzchar(chosen)]

  sample_tmp <- file.path(work_parent, paste0("inner_", gsub("[^A-Za-z0-9_-]", "_", sample_id)))
  if (dir.exists(sample_tmp)) unlink(sample_tmp, recursive = TRUE, force = TRUE)
  dir.create(sample_tmp, recursive = TRUE, showWarnings = FALSE)

  ok <- length(chosen) == 3L && extract_members(archive, chosen, sample_tmp)
  if (!ok) {
    unlink(sample_tmp, recursive = TRUE, force = TRUE)
    return(data.frame(
      gse=gse, sample_id=sample_id, source_archive=basename(archive),
      matrix_file=ifelse(length(matrix_member), matrix_member[1], NA),
      barcode_file=ifelse(length(barcode_member), barcode_member[1], NA),
      feature_file=ifelse(length(feature_member), feature_member[1], NA),
      n_features_matrix=NA, n_barcodes_matrix=NA, n_nonzero=NA,
      n_barcode_lines=NA, n_feature_lines=NA,
      mtx_header_ok=FALSE, dimensions_match=FALSE
    ))
  }

  matrix_path <- locate_one(sample_tmp, "matrix\\.mtx(\\.gz)?$")
  barcode_path <- locate_one(sample_tmp, "barcodes?\\.tsv(\\.gz)?$")
  feature_path <- locate_one(sample_tmp, "(features|genes)\\.tsv(\\.gz)?$")

  res <- audit_triplet_paths(gse, sample_id, matrix_path, barcode_path,
                             feature_path, archive)
  unlink(sample_tmp, recursive = TRUE, force = TRUE)
  res
}

audit_outer_with_nested <- function(gse, outer_file, nested_regex) {
  ds_work <- file.path(WORK_DIR, gse)
  if (dir.exists(ds_work)) unlink(ds_work, recursive = TRUE, force = TRUE)
  dir.create(ds_work, recursive = TRUE, showWarnings = FALSE)

  log_msg("解包外层文件:", basename(outer_file))
  ok <- tryCatch({
    utils::untar(outer_file, exdir = ds_work)
    TRUE
  }, error = function(e) {
    warn_msg("外层解包失败:", basename(outer_file), "|", conditionMessage(e))
    FALSE
  })
  if (!ok) return(data.frame())

  nested <- list.files(ds_work, recursive = TRUE, full.names = TRUE,
                       all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  nested <- nested[grepl(nested_regex, basename(nested), ignore.case = TRUE, perl = TRUE)]

  out <- list()
  for (a in nested) {
    sid <- basename(a)
    sid <- sub("\\.(zip|tar\\.gz|tgz)$", "", sid, ignore.case = TRUE)
    sid <- sub("^GSM[0-9]+_", "", sid)
    sid <- sub("_matrix$", "", sid, ignore.case = TRUE)
    log_msg(gse, "检查嵌套样本:", sid)
    out[[length(out)+1L]] <- audit_nested_10x_archive(gse, a, sid, ds_work)
  }

  unlink(ds_work, recursive = TRUE, force = TRUE)
  if (length(out)) do.call(rbind, out) else data.frame()
}

audit_outer_direct_triplets <- function(gse, outer_file) {
  ds_work <- file.path(WORK_DIR, gse)
  if (dir.exists(ds_work)) unlink(ds_work, recursive = TRUE, force = TRUE)
  dir.create(ds_work, recursive = TRUE, showWarnings = FALSE)

  log_msg("解包直接MTX数据:", basename(outer_file))
  ok <- tryCatch({
    utils::untar(outer_file, exdir = ds_work)
    TRUE
  }, error = function(e) {
    warn_msg("外层解包失败:", basename(outer_file), "|", conditionMessage(e))
    FALSE
  })
  if (!ok) return(data.frame())

  files <- list.files(ds_work, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  matrix_files <- files[grepl("matrix\\.mtx(\\.gz)?$", basename(files), ignore.case = TRUE)]

  out <- list()
  for (mp in matrix_files) {
    stem <- sub("_matrix\\.mtx(\\.gz)?$", "", basename(mp), ignore.case = TRUE)
    bp <- files[grepl(paste0("^", gsub("([.()+*?^$|{}\\[\\]\\\\])", "\\\\\\1", stem),
                              "_barcodes?\\.tsv(\\.gz)?$"),
                      basename(files), ignore.case = TRUE, perl = TRUE)]
    fp <- files[grepl(paste0("^", gsub("([.()+*?^$|{}\\[\\]\\\\])", "\\\\\\1", stem),
                              "_(features|genes)\\.tsv(\\.gz)?$"),
                      basename(files), ignore.case = TRUE, perl = TRUE)]
    sid <- sub("^GSM[0-9]+_", "", stem)
    log_msg(gse, "检查直接样本:", sid)
    out[[length(out)+1L]] <- audit_triplet_paths(
      gse, sid, mp,
      if (length(bp)) bp[1] else NA_character_,
      if (length(fp)) fp[1] else NA_character_,
      outer_file
    )
  }

  unlink(ds_work, recursive = TRUE, force = TRUE)
  if (length(out)) do.call(rbind, out) else data.frame()
}

# ---------------------------- 元数据汇总函数 ----------------------------------
safe_read_table <- function(path, sep = "\t") {
  tryCatch({
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    utils::read.table(con, header = TRUE, sep = sep, quote = "\"",
                      comment.char = "", fill = TRUE, check.names = FALSE)
  }, error = function(e) {
    warn_msg("元数据读取失败:", basename(path), "|", conditionMessage(e))
    NULL
  })
}

table_df <- function(...) {
  x <- as.data.frame(table(..., useNA = "ifany"), stringsAsFactors = FALSE)
  names(x)[ncol(x)] <- "n"
  x
}

parse_soft_samples <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE, encoding = "UTF-8")
  starts <- grep("^\\^SAMPLE\\s*=", lines)
  if (!length(starts)) return(data.frame())
  ends <- c(starts[-1]-1L, length(lines))

  get1 <- function(block, key) {
    h <- grep(key, block, ignore.case = TRUE)
    if (!length(h)) return(NA_character_)
    sub("^[^=]*=\\s*", "", block[h[1]])
  }
  getall <- function(block, key) {
    h <- grep(key, block, ignore.case = TRUE)
    if (!length(h)) return(NA_character_)
    paste(sub("^[^=]*=\\s*", "", block[h]), collapse = " | ")
  }

  rows <- vector("list", length(starts))
  for (i in seq_along(starts)) {
    b <- lines[starts[i]:ends[i]]
    rows[[i]] <- data.frame(
      gsm = sub("^\\^SAMPLE\\s*=\\s*", "", b[1]),
      title = get1(b, "^!Sample_title\\s*="),
      characteristics = getall(b, "^!Sample_characteristics_ch1"),
      description = getall(b, "^!Sample_description"),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# ------------------------------- 开始审计 --------------------------------------
log_msg("第二阶段开始")
if (!dir.exists(ROOT_DIR)) stop("根目录不存在：", ROOT_DIR)

matrix_results <- list()

# GSE241132：12个嵌套ZIP
f <- find_file_exact("GSE241132_RAW.tar")
if (!is.na(f)) matrix_results[["GSE241132"]] <- audit_outer_with_nested(
  "GSE241132", f, "\\.zip$"
) else warn_msg("缺少GSE241132_RAW.tar")

# GSE163973：6个嵌套TAR.GZ
f <- find_file_exact("GSE163973_RAW.tar")
if (!is.na(f)) matrix_results[["GSE163973"]] <- audit_outer_with_nested(
  "GSE163973", f, "\\.tar\\.gz$"
) else warn_msg("缺少GSE163973_RAW.tar")

# GSE265972：9个嵌套ZIP
f <- find_file_exact("GSE265972_processed.tar.gz")
if (!is.na(f)) matrix_results[["GSE265972"]] <- audit_outer_with_nested(
  "GSE265972", f, "\\.zip$"
) else warn_msg("缺少GSE265972_processed.tar.gz")

# GSE220300：11套直接MTX
f <- find_file_exact("GSE220300_RAW.tar")
if (!is.na(f)) matrix_results[["GSE220300"]] <- audit_outer_direct_triplets(
  "GSE220300", f
) else warn_msg("缺少GSE220300_RAW.tar")

# GSE181297：直接MTX（scRNA+空间）
f <- find_file_exact("GSE181297_RAW.tar")
if (!is.na(f)) matrix_results[["GSE181297"]] <- audit_outer_direct_triplets(
  "GSE181297", f
) else warn_msg("缺少GSE181297_RAW.tar")

# 新增推荐队列GSE181316：可选但强烈建议
f <- find_file_exact("GSE181316_RAW.tar")
if (!is.na(f)) matrix_results[["GSE181316"]] <- audit_outer_direct_triplets(
  "GSE181316", f
) else warn_msg("尚未发现GSE181316_RAW.tar；第二阶段仍可完成，但建议补充后重跑。")

matrix_audit <- if (length(matrix_results)) do.call(rbind, matrix_results) else data.frame()
safe_write_csv(matrix_audit, file.path(REPORT_DIR, "01_matrix_dimension_audit.csv"))

# 按数据集汇总矩阵通过率
if (nrow(matrix_audit)) {
  matrix_summary <- do.call(rbind, lapply(split(matrix_audit, matrix_audit$gse), function(x) {
    data.frame(
      gse = x$gse[1],
      n_samples_checked = nrow(x),
      n_valid_triplets = sum(x$dimensions_match, na.rm = TRUE),
      all_samples_valid = all(x$dimensions_match),
      min_cells_or_spots = suppressWarnings(min(x$n_barcodes_matrix, na.rm = TRUE)),
      max_cells_or_spots = suppressWarnings(max(x$n_barcodes_matrix, na.rm = TRUE)),
      min_features = suppressWarnings(min(x$n_features_matrix, na.rm = TRUE)),
      max_features = suppressWarnings(max(x$n_features_matrix, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(matrix_summary) <- NULL
} else {
  matrix_summary <- data.frame()
}
safe_write_csv(matrix_summary, file.path(REPORT_DIR, "02_matrix_dataset_summary.csv"))

# ------------------------ GSE241124空间文件审计 -------------------------------
spatial_outer <- find_file_exact("GSE241124_RAW.tar")
spatial_files <- data.frame()
if (!is.na(spatial_outer)) {
  ds_work <- file.path(WORK_DIR, "GSE241124")
  dir.create(ds_work, recursive = TRUE, showWarnings = FALSE)
  log_msg("解包GSE241124空间文件")
  tryCatch(utils::untar(spatial_outer, exdir = ds_work),
           error = function(e) warn_msg("GSE241124解包失败:", conditionMessage(e)))

  files <- list.files(ds_work, recursive = TRUE, full.names = TRUE,
                      all.files = TRUE, include.dirs = FALSE, no.. = TRUE)
  h5 <- files[grepl("\\.h5$", files, ignore.case = TRUE)]
  zips <- files[grepl("spatial.*\\.zip$", basename(files), ignore.case = TRUE)]

  rows <- list()
  for (zp in zips) {
    mem <- list_archive(zp)
    b <- sub("_spatial_images\\.zip$", "", basename(zp), ignore.case = TRUE)
    h5_match <- h5[startsWith(tolower(basename(h5)), tolower(b))]
    rows[[length(rows)+1L]] <- data.frame(
      sample_prefix = b,
      h5_present = length(h5_match) > 0L,
      h5_size_mb = if (length(h5_match)) round(file.info(h5_match[1])$size/1024^2, 3) else NA_real_,
      spatial_zip_present = TRUE,
      n_zip_members = length(mem),
      has_tissue_positions = any(grepl("tissue_positions", mem, ignore.case = TRUE)),
      has_scalefactors = any(grepl("scalefactors.*json", mem, ignore.case = TRUE)),
      has_hires_image = any(grepl("hires.*\\.(png|jpg|jpeg)$", mem, ignore.case = TRUE)),
      has_lowres_image = any(grepl("lowres.*\\.(png|jpg|jpeg)$", mem, ignore.case = TRUE)),
      stringsAsFactors = FALSE
    )
  }
  spatial_files <- if (length(rows)) do.call(rbind, rows) else data.frame()
  unlink(ds_work, recursive = TRUE, force = TRUE)
} else {
  warn_msg("缺少GSE241124_RAW.tar")
}
safe_write_csv(spatial_files, file.path(REPORT_DIR, "03_GSE241124_spatial_structure.csv"))

# ------------------------ GSE181297空间完整性专项 ------------------------------
g181_members <- data.frame()
f181 <- find_file_exact("GSE181297_RAW.tar")
if (!is.na(f181)) {
  mem <- tryCatch(utils::untar(f181, list = TRUE), error = function(e) character(0))
  g181_members <- data.frame(
    member = mem,
    is_h5 = grepl("\\.h5$", mem, ignore.case = TRUE),
    is_image = grepl("\\.(png|jpg|jpeg)(\\.gz)?$", mem, ignore.case = TRUE),
    is_tissue_positions = grepl("tissue_positions", mem, ignore.case = TRUE),
    is_scalefactors = grepl("scalefactors", mem, ignore.case = TRUE),
    stringsAsFactors = FALSE
  )
}
safe_write_csv(g181_members, file.path(REPORT_DIR, "04_GSE181297_spatial_member_check.csv"))

# ---------------------------- 元数据完整汇总 ----------------------------------
# GSE241132
f <- find_file_exact("GSE241132_cell_metadata.txt.gz")
if (!is.na(f)) {
  d <- safe_read_table(f, "\t")
  if (!is.null(d)) {
    req <- intersect(c("orig.ident","Patient","Gender","Age","Condition",
                       "Doublet","newMainCellTypes","newCellTypes"), names(d))
    safe_write_csv(d[seq_len(min(1000L, nrow(d))), req, drop=FALSE],
                   file.path(REPORT_DIR, "05_GSE241132_metadata_first1000.csv"))
    if (all(c("Patient","Condition","orig.ident") %in% names(d))) {
      safe_write_csv(table_df(d$Patient, d$Condition, d$orig.ident),
                     file.path(REPORT_DIR, "05_GSE241132_sample_cell_counts.csv"))
    }
    if (all(c("Patient","Condition","newMainCellTypes") %in% names(d))) {
      safe_write_csv(table_df(d$Patient, d$Condition, d$newMainCellTypes),
                     file.path(REPORT_DIR, "05_GSE241132_celltype_counts.csv"))
    }
    if ("Doublet" %in% names(d)) {
      safe_write_csv(table_df(d$orig.ident, d$Doublet),
                     file.path(REPORT_DIR, "05_GSE241132_doublet_counts.csv"))
    }
  }
}

# GSE241124
f <- find_file_exact("GSE241124_spatialseq_metadata_acutewound.txt.gz")
if (!is.na(f)) {
  d <- safe_read_table(f, "\t")
  if (!is.null(d)) {
    req <- intersect(c("barcode","orig.ident","Patient","Donor","Gender","Age",
                       "Seq_Batch","Condition","Sample_name","AnnoType",
                       "nCount_Spatial","nFeature_Spatial"), names(d))
    safe_write_csv(d[seq_len(min(1000L,nrow(d))), req, drop=FALSE],
                   file.path(REPORT_DIR, "06_GSE241124_metadata_first1000.csv"))
    grp <- intersect(c("Donor","Condition","Sample_name"), names(d))
    if (length(grp) == 3L) {
      safe_write_csv(table_df(d$Donor, d$Condition, d$Sample_name),
                     file.path(REPORT_DIR, "06_GSE241124_spot_counts.csv"))
    }
  }
}

# GSE163973
f <- find_file_exact("GSE163973_integrate.all.NS.all.KL_cell.meta.data.csv.gz")
if (!is.na(f)) {
  d <- safe_read_table(f, ",")
  if (!is.null(d)) {
    req <- intersect(c("orig.ident","dataset","condition","cellType",
                       "nCount_RNA","nFeature_RNA","percent.mt"), names(d))
    safe_write_csv(d[seq_len(min(1000L,nrow(d))), req, drop=FALSE],
                   file.path(REPORT_DIR, "07_GSE163973_metadata_first1000.csv"))
    if (all(c("orig.ident","condition") %in% names(d))) {
      safe_write_csv(table_df(d$orig.ident, d$condition),
                     file.path(REPORT_DIR, "07_GSE163973_sample_counts.csv"))
    }
    if (all(c("orig.ident","condition","cellType") %in% names(d))) {
      safe_write_csv(table_df(d$orig.ident, d$condition, d$cellType),
                     file.path(REPORT_DIR, "07_GSE163973_celltype_counts.csv"))
    }
  }
}

# GSE265972
f <- find_file_exact("GSE265972_VU_anno_metadata.txt.gz")
if (!is.na(f)) {
  d <- safe_read_table(f, "\t")
  if (!is.null(d)) {
    req <- intersect(c("orig.ident","Condition","samples","barcode","ct","keep",
                       "CellType","nCount_RNA","nFeature_RNA","percent.mt"), names(d))
    safe_write_csv(d[seq_len(min(1000L,nrow(d))), req, drop=FALSE],
                   file.path(REPORT_DIR, "08_GSE265972_metadata_first1000.csv"))
    if (all(c("samples","Condition") %in% names(d))) {
      safe_write_csv(table_df(d$samples, d$Condition),
                     file.path(REPORT_DIR, "08_GSE265972_sample_counts.csv"))
    }
    if (all(c("samples","Condition","CellType") %in% names(d))) {
      safe_write_csv(table_df(d$samples, d$Condition, d$CellType),
                     file.path(REPORT_DIR, "08_GSE265972_celltype_counts.csv"))
    }
  }
}

# -------------------------- SOFT与样本字典 ------------------------------------
soft_paths <- find_file_pattern("family\\.soft\\.gz$")
sample_dict_list <- list()
for (sp in soft_paths) {
  g <- extract_gse(sp)
  x <- tryCatch(parse_soft_samples(sp), error=function(e) {
    warn_msg("SOFT解析失败:", basename(sp), "|", conditionMessage(e))
    data.frame()
  })
  if (!nrow(x)) next
  x$gse <- g

  # 预设字段
  x$group <- NA_character_
  x$patient_id <- NA_character_
  x$region <- NA_character_
  x$activity <- NA_character_

  if (g == "GSE220300") {
    x$group <- ifelse(grepl("^AC|^AP", x$title), "active_keloid",
               ifelse(grepl("^IC|^IP", x$title), "inactive_keloid",
               ifelse(grepl("^MS", x$title), "mature_scar",
               ifelse(grepl("^NS", x$title), "normal_scar", NA))))
    x$patient_id <- ifelse(grepl("^[AI][CP][0-9]+$", x$title),
                           paste0(substr(x$title,1,1), sub("^[AI][CP]", "", x$title)),
                           x$title)
    x$region <- ifelse(grepl("^[AI]C", x$title), "center",
                ifelse(grepl("^[AI]P", x$title), "periphery", NA))
    x$activity <- ifelse(grepl("^A", x$title), "active",
                  ifelse(grepl("^I", x$title), "inactive", NA))
  } else if (g == "GSE163973") {
    x$group <- ifelse(grepl("keloid", x$title, ignore.case=TRUE), "keloid", "normal_scar")
    x$patient_id <- x$title
  } else if (g == "GSE181297") {
    x$group <- ifelse(x$title %in% c("Ke01","Ke02","Pt1","Pt2"), "keloid",
               ifelse(x$title=="NS02", "normal_scar", "adjacent_normal"))
    x$patient_id <- x$title
  } else if (g == "GSE265972") {
    x$group <- ifelse(grepl("^VU", x$title), "venous_ulcer", "healthy_skin")
    x$patient_id <- x$title
  } else if (g == "GSE241132" || g == "GSE241124") {
    x$group <- ifelse(grepl("Wound1", x$title), "day1",
               ifelse(grepl("Wound7", x$title), "day7",
               ifelse(grepl("Wound30", x$title), "day30", "skin")))
    x$patient_id <- sub("[, -].*$", "", x$title)
  } else if (g == "GSE181316") {
    x$group <- ifelse(grepl("^keloid", x$title, ignore.case=TRUE), "keloid",
               ifelse(grepl("^scar", x$title, ignore.case=TRUE), "normal_scar", "healthy_skin"))
    x$patient_id <- ifelse(grepl("^keloid_3", x$title, ignore.case=TRUE),
                           "keloid_patient3", x$title)
  } else if (g == "GSE274709") {
    x$group <- ifelse(grepl("Keloid", x$title, ignore.case=TRUE), "keloid",
               ifelse(grepl("_ml", x$title, ignore.case=TRUE), "control_group1", "control_group2"))
    x$patient_id <- x$title
  }

  sample_dict_list[[length(sample_dict_list)+1L]] <- x
}
sample_dictionary <- if (length(sample_dict_list)) do.call(rbind, sample_dict_list) else data.frame()
safe_write_csv(sample_dictionary, file.path(REPORT_DIR, "09_master_sample_dictionary.csv"))

# -------------------------- GSE274709最终判断 ---------------------------------
g274_outer <- find_file_exact("GSE274709_RAW.tar")
g274_series <- find_file_pattern("^GSE274709_series_matrix.*\\.gz$")
g274_result <- data.frame(
  outer_present = !is.na(g274_outer),
  n_outer_members = NA_integer_,
  n_idxstats = NA_integer_,
  has_non_idxstats_expression = FALSE,
  series_matrix_present = length(g274_series) > 0L,
  series_matrix_data_rows = NA_integer_,
  usable_for_expression = FALSE,
  stringsAsFactors = FALSE
)
if (!is.na(g274_outer)) {
  mem <- tryCatch(utils::untar(g274_outer, list=TRUE), error=function(e) character(0))
  g274_result$n_outer_members <- length(mem)
  g274_result$n_idxstats <- sum(grepl("idxstats", mem, ignore.case=TRUE))
  g274_result$has_non_idxstats_expression <- any(
    grepl("count|expression|tpm|fpkm|matrix", mem, ignore.case=TRUE) &
      !grepl("idxstats", mem, ignore.case=TRUE)
  )
}
if (length(g274_series)) {
  con <- gzfile(g274_series[1], "rt")
  ln <- readLines(con, warn=FALSE)
  close(con)
  b <- grep("^!series_matrix_table_begin", ln)
  e <- grep("^!series_matrix_table_end", ln)
  if (length(b) && length(e)) {
    g274_result$series_matrix_data_rows <- max(0L, e[1] - b[1] - 2L)
  }
}
g274_result$usable_for_expression <- isTRUE(g274_result$has_non_idxstats_expression) ||
  (!is.na(g274_result$series_matrix_data_rows) && g274_result$series_matrix_data_rows > 0L)
safe_write_csv(g274_result, file.path(REPORT_DIR, "10_GSE274709_usability.csv"))

# ---------------------------- 第二阶段决策 ------------------------------------
get_summary_row <- function(g) {
  if (!nrow(matrix_summary) || !g %in% matrix_summary$gse) return(NULL)
  matrix_summary[matrix_summary$gse==g,,drop=FALSE]
}
s241132 <- get_summary_row("GSE241132")
s220300 <- get_summary_row("GSE220300")
s163973 <- get_summary_row("GSE163973")

core_matrix_pass <- !is.null(s241132) && isTRUE(s241132$all_samples_valid) &&
                    !is.null(s220300) && isTRUE(s220300$all_samples_valid) &&
                    !is.null(s163973) && isTRUE(s163973$all_samples_valid)

spatial_pass <- nrow(spatial_files)==16L &&
  all(spatial_files$h5_present) &&
  all(spatial_files$has_tissue_positions) &&
  all(spatial_files$has_scalefactors) &&
  all(spatial_files$has_lowres_image | spatial_files$has_hires_image)

g181_spatial_complete <- nrow(g181_members)>0 &&
  any(g181_members$is_tissue_positions) &&
  any(g181_members$is_scalefactors)

has_g181316 <- "GSE181316" %in% matrix_summary$gse &&
  isTRUE(matrix_summary$all_samples_valid[matrix_summary$gse=="GSE181316"][1])

decision_lines <- c(
  "第二阶段：矩阵结构与样本映射审计结论",
  "============================================================",
  paste0("运行时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("核心单细胞矩阵是否全部通过：", core_matrix_pass),
  paste0("正常伤口空间文件是否完整：", spatial_pass),
  paste0("GSE181297是否具备标准空间坐标与比例因子：", g181_spatial_complete),
  paste0("是否已经补充并通过GSE181316：", has_g181316),
  paste0("GSE274709能否直接做表达/复发分析：", g274_result$usable_for_expression),
  "",
  "设计层面的固定结论：",
  "- GSE220300只有2名活动性患者和2名非活动性患者的中心/周边样本；活动度分析只能作为探索性证据。",
  "- GSE163973为3例keloid与3例normal scar，可作为独立验证之一。",
  "- GSE274709若usable_for_expression=FALSE，应从计算性复发验证中移除。",
  "- 论文主要确认性问题应改为：正常伤口状态在多个独立瘢痕疙瘩队列中的持续与偏轨。",
  "- 活动度、周边/中央和复发只作为临床支持，不作为唯一主终点。",
  "",
  "进入第三阶段（正式QC）的条件：",
  "1. GSE241132、GSE220300、GSE163973矩阵维度全部匹配；",
  "2. GSE241124的16个H5与空间坐标/图像结构完整；",
  "3. 元数据可恢复供者、阶段和分组；",
  "4. 强烈建议补充GSE181316，以增加独立瘢痕疙瘩患者级证据。"
)
safe_write_lines(decision_lines, file.path(REPORT_DIR, "11_STAGE2_DECISION.txt"))

capture.output(sessionInfo(), file=file.path(REPORT_DIR, "12_SESSION_INFO.txt"))
safe_write_lines(
  c(paste0("STAGE2_COMPLETED=",format(Sys.time(),"%Y-%m-%d %H:%M:%S")),
    paste0("CORE_MATRIX_PASS=",core_matrix_pass),
    paste0("NORMAL_SPATIAL_PASS=",spatial_pass),
    paste0("GSE181316_PRESENT_AND_VALID=",has_g181316),
    paste0("GSE274709_USABLE=",g274_result$usable_for_expression)),
  file.path(REPORT_DIR, "STAGE2_COMPLETED.txt")
)

# 删除临时工作目录，检查包仅保留报告
if (dir.exists(WORK_DIR)) unlink(WORK_DIR, recursive=TRUE, force=TRUE)

# ------------------------------- 自动压缩 --------------------------------------
create_zip_windows <- function(source_dir, zip_path) {
  if (file.exists(zip_path)) file.remove(zip_path)
  esc <- function(x) gsub("'", "''", normalize_slash(x), fixed=TRUE)
  cmd <- paste0(
    "$ErrorActionPreference='Stop'; ",
    "$items=Get-ChildItem -LiteralPath '",esc(source_dir),"'; ",
    "Compress-Archive -Path $items.FullName -DestinationPath '",esc(zip_path),"' -Force"
  )
  try(system2("powershell.exe",
              c("-NoProfile","-ExecutionPolicy","Bypass","-Command",cmd),
              stdout=TRUE, stderr=TRUE), silent=TRUE)
  file.exists(zip_path) && file.info(zip_path)$size>0
}

zip_ok <- create_zip_windows(STAGE_DIR, PACKAGE_ZIP)
if (!zip_ok) {
  if (file.exists(PACKAGE_TARGZ)) file.remove(PACKAGE_TARGZ)
  old <- getwd()
  setwd(ROOT_DIR)
  try(utils::tar(basename(PACKAGE_TARGZ), basename(STAGE_DIR),
                 compression="gzip", tar="internal"), silent=TRUE)
  setwd(old)
}

log_msg("第二阶段完成")
cat("\n============================================================\n")
cat("第二阶段运行完成。\n")
if (file.exists(PACKAGE_ZIP)) {
  cat("请上传：", PACKAGE_ZIP, "\n", sep="")
} else if (file.exists(PACKAGE_TARGZ)) {
  cat("请上传备用包：", PACKAGE_TARGZ, "\n", sep="")
} else {
  cat("请手动压缩：", STAGE_DIR, "\n", sep="")
}
cat("============================================================\n")
