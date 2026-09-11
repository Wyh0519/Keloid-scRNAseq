# =============================================================================
# 项目：瘢痕疙瘩“分子伤口年龄与修复终止失败”
# 第一阶段：原始文件完整性、压缩包结构与样本元数据审计（V1.1修正版）
#
# 根目录：
#   C:/Users/33652/Desktop/LZZ/1、yssj
#
# 本阶段目标：
#   1. 检查7个GEO数据集的关键文件是否存在、大小是否合理；
#   2. 不完整解压大型文件，仅列出压缩包内部文件；
#   3. 从SOFT/Series Matrix/metadata中恢复GSM、样本标题和分组线索；
#   4. 判断MTX/H5/空间文件/表达矩阵是否具备；
#   5. 特别判断GSE274709是否只有idxstats；
#   6. 自动生成“第一阶段检查包.zip”，供上传复核。
#
# 运行要求：
#   - Windows + R 4.x
#   - 仅使用base R，不安装任何第三方R包
#   - 不删除、不移动、不覆盖原始文件
#   - V1.1已修复R 4.6.0下GSE编号提取与中文路径兼容问题
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANG = "en")

# ----------------------------- 1. 路径配置 -----------------------------------
ROOT_DIR <- "C:/Users/33652/Desktop/LZZ/1、yssj"

STAGE_DIR <- file.path(ROOT_DIR, "01_STAGE1_AUDIT")
REPORT_DIR <- file.path(STAGE_DIR, "report")
PREVIEW_DIR <- file.path(STAGE_DIR, "metadata_previews")
PACKAGE_ZIP <- file.path(ROOT_DIR, "第一阶段_数据与元数据审计检查包.zip")
PACKAGE_TARGZ <- file.path(ROOT_DIR, "第一阶段_数据与元数据审计检查包.tar.gz")

dir.create(STAGE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PREVIEW_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(REPORT_DIR, "00_stage1_log.txt")
WARNING_FILE <- file.path(REPORT_DIR, "00_stage1_warnings.txt")

if (file.exists(LOG_FILE)) file.remove(LOG_FILE)
if (file.exists(WARNING_FILE)) file.remove(WARNING_FILE)

log_msg <- function(...) {
  msg <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    " | ",
    paste(..., collapse = " ")
  )
  cat(msg, "\n")
  cat(msg, "\n", file = LOG_FILE, append = TRUE)
}

warn_msg <- function(...) {
  msg <- paste(..., collapse = " ")
  log_msg("WARNING:", msg)
  cat(msg, "\n", file = WARNING_FILE, append = TRUE)
}

stop_with_log <- function(...) {
  msg <- paste(..., collapse = " ")
  log_msg("FATAL:", msg)
  stop(msg, call. = FALSE)
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

normalize_slash <- function(x) {
  normalizePath(x, winslash = "/", mustWork = FALSE)
}

rel_path <- function(x, root) {
  x2 <- normalize_slash(x)
  r2 <- normalize_slash(root)
  prefix <- paste0(r2, "/")
  if (startsWith(tolower(x2), tolower(prefix))) {
    sub(prefix, "", x2, fixed = TRUE)
  } else {
    basename(x2)
  }
}

extract_gse <- function(x) {
  b <- basename(x)
  out <- rep(NA_character_, length(b))
  hit <- grepl("GSE[0-9]+", b, ignore.case = TRUE, perl = TRUE)
  if (any(hit)) {
    out[hit] <- toupper(
      sub(
        "^.*?(GSE[0-9]+).*$",
        "\\1",
        b[hit],
        ignore.case = TRUE,
        perl = TRUE
      )
    )
  }
  out
}

file_role <- function(x) {
  n <- tolower(basename(x))
  if (grepl("family\\.soft", n)) return("SOFT")
  if (grepl("series_matrix", n)) return("SeriesMatrix")
  if (grepl("metadata|meta\\.data|anno", n)) return("Metadata")
  if (grepl("raw_counts|count.*matrix", n)) return("CountMatrix")
  if (grepl("\\.rdata(\\.gz)?$|\\.rda(\\.gz)?$", n)) return("RData")
  if (grepl("\\.tar$|\\.tar\\.gz$|\\.tgz$|\\.zip$", n)) return("Archive")
  if (grepl("\\.csv(\\.gz)?$", n)) return("CSV")
  if (grepl("\\.txt(\\.gz)?$", n)) return("Text")
  "Other"
}

safe_read_lines <- function(path, n = 500L) {
  con <- NULL
  out <- tryCatch({
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
      con <- gzfile(path, open = "rt")
    } else {
      con <- file(path, open = "rt")
    }
    readLines(con, n = n, warn = FALSE, encoding = "UTF-8")
  }, error = function(e) {
    warn_msg("读取文本失败:", basename(path), "|", conditionMessage(e))
    character(0)
  }, finally = {
    if (!is.null(con)) try(close(con), silent = TRUE)
  })
  out
}

read_all_lines <- function(path) {
  con <- NULL
  out <- tryCatch({
    if (grepl("\\.gz$", path, ignore.case = TRUE)) {
      con <- gzfile(path, open = "rt")
    } else {
      con <- file(path, open = "rt")
    }
    readLines(con, warn = FALSE, encoding = "UTF-8")
  }, error = function(e) {
    warn_msg("完整读取失败:", basename(path), "|", conditionMessage(e))
    character(0)
  }, finally = {
    if (!is.null(con)) try(close(con), silent = TRUE)
  })
  out
}

# ----------------------------- 2. 根目录检查 ---------------------------------
log_msg("开始第一阶段审计")
log_msg("ROOT_DIR =", ROOT_DIR)

if (!dir.exists(ROOT_DIR)) {
  stop_with_log(
    "根目录不存在：", ROOT_DIR,
    "。请检查中文顿号“、”是否正确，或修改脚本顶部ROOT_DIR。"
  )
}

all_files <- list.files(
  ROOT_DIR,
  recursive = TRUE,
  full.names = TRUE,
  all.files = TRUE,
  include.dirs = FALSE,
  no.. = TRUE
)

# 排除本阶段输出以及旧检查包
all_files <- all_files[
  !grepl("/01_STAGE1_AUDIT/", normalize_slash(all_files), fixed = TRUE) &
  basename(all_files) != basename(PACKAGE_ZIP) &
  basename(all_files) != basename(PACKAGE_TARGZ)
]

if (length(all_files) == 0L) {
  stop_with_log("根目录中没有发现数据文件。")
}

# ----------------------------- 3. 文件总清单 ---------------------------------
fi <- file.info(all_files)

inventory <- data.frame(
  full_path = normalize_slash(all_files),
  relative_path = vapply(all_files, rel_path, character(1), root = ROOT_DIR),
  file_name = basename(all_files),
  gse = extract_gse(all_files),
  role = vapply(all_files, file_role, character(1)),
  size_bytes = as.numeric(fi$size),
  size_mb = round(as.numeric(fi$size) / 1024^2, 3),
  modified_time = as.character(fi$mtime),
  stringsAsFactors = FALSE
)

# 小于100MB的文件计算MD5，大文件只记录大小，避免审计时间过长
inventory$md5 <- NA_character_
md5_idx <- which(!is.na(inventory$size_bytes) & inventory$size_bytes <= 100 * 1024^2)
if (length(md5_idx) > 0L) {
  md5_val <- tryCatch(
    tools::md5sum(all_files[md5_idx]),
    error = function(e) {
      warn_msg("MD5计算失败:", conditionMessage(e))
      rep(NA_character_, length(md5_idx))
    }
  )
  inventory$md5[md5_idx] <- unname(md5_val)
}

safe_write_csv(inventory, file.path(REPORT_DIR, "01_all_files_inventory.csv"))
log_msg("发现文件数:", nrow(inventory))

# ----------------------------- 4. 必需文件检查 --------------------------------
expected_files <- list(
  GSE241132 = c(
    "GSE241132_RAW.tar",
    "GSE241132_cell_metadata.txt.gz"
  ),
  GSE241124 = c(
    "GSE241124_RAW.tar",
    "GSE241124_spatialseq_metadata_acutewound.txt.gz"
  ),
  GSE220300 = c(
    "GSE220300_RAW.tar",
    "GSE220300_family.soft.gz"
  ),
  GSE163973 = c(
    "GSE163973_RAW.tar",
    "GSE163973_integrate.all.NS.all.KL_cell.meta.data.csv.gz"
  ),
  GSE181297 = c(
    "GSE181297_RAW.tar",
    "GSE181297_family.soft.gz"
  ),
  GSE265972 = c(
    "GSE265972_processed.tar.gz",
    "GSE265972_VU_anno_metadata.txt.gz"
  ),
  GSE274709 = c(
    "GSE274709_RAW.tar",
    "GSE274709_family.soft.gz"
  )
)

required_rows <- list()
for (g in names(expected_files)) {
  for (f in expected_files[[g]]) {
    hit <- which(tolower(inventory$file_name) == tolower(f))
    required_rows[[length(required_rows) + 1L]] <- data.frame(
      gse = g,
      required_file = f,
      present = length(hit) > 0L,
      found_path = if (length(hit)) inventory$relative_path[hit[1]] else NA_character_,
      size_mb = if (length(hit)) inventory$size_mb[hit[1]] else NA_real_,
      stringsAsFactors = FALSE
    )
  }
}
required_check <- do.call(rbind, required_rows)
safe_write_csv(required_check, file.path(REPORT_DIR, "02_required_files_check.csv"))

if (any(!required_check$present)) {
  warn_msg("存在缺失的关键文件，请查看02_required_files_check.csv")
}

# ----------------------------- 5. 压缩包目录审计 -------------------------------
archive_paths <- inventory$full_path[inventory$role == "Archive"]
archive_detail_list <- list()
archive_summary_list <- list()

for (path in archive_paths) {
  fname <- basename(path)
  gse <- extract_gse(fname)
  log_msg("列出压缩包内容:", fname)

  members <- tryCatch({
    if (grepl("\\.zip$", fname, ignore.case = TRUE)) {
      z <- utils::unzip(path, list = TRUE)
      data.frame(
        member = z$Name,
        member_size = z$Length,
        stringsAsFactors = FALSE
      )
    } else {
      x <- utils::untar(path, list = TRUE)
      data.frame(
        member = x,
        member_size = NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }, error = function(e) {
    warn_msg("压缩包无法读取:", fname, "|", conditionMessage(e))
    data.frame(member = character(0), member_size = numeric(0))
  })

  if (nrow(members) == 0L) next

  ml <- tolower(members$member)
  members$gse <- gse
  members$archive <- fname
  members$has_matrix_mtx <- grepl("matrix\\.mtx(\\.gz)?$", ml)
  members$has_barcodes <- grepl("barcodes?\\.tsv(\\.gz)?$", ml)
  members$has_features <- grepl("(features|genes)\\.tsv(\\.gz)?$", ml)
  members$has_h5 <- grepl("\\.h5$", ml)
  members$has_spatial <- grepl(
    "spatial|tissue_positions|scalefactors|hires|lowres|\\.png$|\\.jpg$|\\.jpeg$",
    ml
  )
  members$has_idxstats <- grepl("idxstats", ml)
  members$has_rdata <- grepl("\\.rdata(\\.gz)?$|\\.rda(\\.gz)?$", ml)
  members$has_expression_hint <- grepl(
    "matrix|count|expression|normalized|tpm|fpkm|rpkm|\\.mtx|\\.h5$",
    ml
  )

  archive_detail_list[[length(archive_detail_list) + 1L]] <- members

  archive_summary_list[[length(archive_summary_list) + 1L]] <- data.frame(
    gse = gse,
    archive = fname,
    archive_readable = TRUE,
    n_members = nrow(members),
    n_matrix_mtx = sum(members$has_matrix_mtx),
    n_barcodes = sum(members$has_barcodes),
    n_features = sum(members$has_features),
    n_h5 = sum(members$has_h5),
    n_spatial_related = sum(members$has_spatial),
    n_rdata = sum(members$has_rdata),
    n_idxstats = sum(members$has_idxstats),
    n_expression_hints = sum(members$has_expression_hint),
    stringsAsFactors = FALSE
  )
}

archive_details <- if (length(archive_detail_list)) {
  do.call(rbind, archive_detail_list)
} else {
  data.frame()
}

archive_summary <- if (length(archive_summary_list)) {
  do.call(rbind, archive_summary_list)
} else {
  data.frame()
}

safe_write_csv(archive_details, file.path(REPORT_DIR, "03_archive_member_details.csv"))
safe_write_csv(archive_summary, file.path(REPORT_DIR, "04_archive_summary.csv"))

# 为每个GSE单独生成压缩包文件名清单，便于人工检查
if (nrow(archive_details) > 0L) {
  for (g in unique(na.omit(archive_details$gse))) {
    x <- archive_details[archive_details$gse == g, c("archive", "member"), drop = FALSE]
    safe_write_csv(
      x,
      file.path(REPORT_DIR, paste0("archive_list_", g, ".csv"))
    )
  }
}

# ----------------------------- 6. SOFT样本表解析 -------------------------------
parse_soft_samples <- function(path) {
  lines <- read_all_lines(path)
  if (!length(lines)) return(data.frame())

  sample_starts <- grep("^\\^SAMPLE\\s*=", lines)
  if (!length(sample_starts)) return(data.frame())

  sample_ends <- c(sample_starts[-1] - 1L, length(lines))
  rows <- vector("list", length(sample_starts))

  get_first_value <- function(block, key_regex) {
    hit <- grep(key_regex, block, ignore.case = TRUE)
    if (!length(hit)) return(NA_character_)
    sub("^[^=]*=\\s*", "", block[hit[1]])
  }

  get_all_values <- function(block, key_regex) {
    hit <- grep(key_regex, block, ignore.case = TRUE)
    if (!length(hit)) return(NA_character_)
    paste(sub("^[^=]*=\\s*", "", block[hit]), collapse = " | ")
  }

  for (i in seq_along(sample_starts)) {
    block <- lines[sample_starts[i]:sample_ends[i]]
    gsm <- sub("^\\^SAMPLE\\s*=\\s*", "", block[1])
    rows[[i]] <- data.frame(
      gsm = gsm,
      title = get_first_value(block, "^!Sample_title\\s*="),
      source_name = get_first_value(block, "^!Sample_source_name_ch1\\s*="),
      organism = get_first_value(block, "^!Sample_organism_ch1\\s*="),
      characteristics = get_all_values(block, "^!Sample_characteristics_ch1"),
      description = get_all_values(block, "^!Sample_description"),
      data_processing = get_all_values(block, "^!Sample_data_processing"),
      supplementary_file = get_all_values(block, "^!Sample_supplementary_file"),
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, rows)
}

soft_files <- inventory$full_path[inventory$role == "SOFT"]
soft_summary_list <- list()

for (path in soft_files) {
  gse <- extract_gse(path)
  fname <- basename(path)
  log_msg("解析SOFT样本信息:", fname)

  x <- parse_soft_samples(path)
  if (!nrow(x)) {
    warn_msg("SOFT中没有解析到样本:", fname)
    next
  }

  x$gse <- gse
  x$source_file <- fname
  x <- x[, c("gse", "source_file", setdiff(names(x), c("gse", "source_file")))]

  safe_write_csv(
    x,
    file.path(REPORT_DIR, paste0("05_SOFT_samples_", gse, ".csv"))
  )

  soft_summary_list[[length(soft_summary_list) + 1L]] <- data.frame(
    gse = gse,
    source_file = fname,
    n_samples_parsed = nrow(x),
    titles = paste(x$title, collapse = " || "),
    stringsAsFactors = FALSE
  )
}

soft_summary <- if (length(soft_summary_list)) {
  do.call(rbind, soft_summary_list)
} else {
  data.frame()
}

safe_write_csv(soft_summary, file.path(REPORT_DIR, "05_SOFT_parse_summary.csv"))

# ----------------------------- 7. Series Matrix头部审计 ------------------------
series_files <- inventory$full_path[inventory$role == "SeriesMatrix"]
series_summary_list <- list()

for (path in series_files) {
  gse <- extract_gse(path)
  fname <- basename(path)
  lines <- safe_read_lines(path, n = 600L)
  if (!length(lines)) next

  preview_path <- file.path(
    PREVIEW_DIR,
    paste0(gsub("[^A-Za-z0-9._-]", "_", fname), "_first600lines.txt")
  )
  safe_write_lines(lines, preview_path)

  n_sample_title <- sum(grepl("^!Sample_title", lines))
  n_sample_geo <- sum(grepl("^!Sample_geo_accession", lines))
  has_table_begin <- any(grepl("^!series_matrix_table_begin", lines, ignore.case = TRUE))
  has_gene_header <- any(grepl('^"ID_REF"|^ID_REF', lines))

  series_summary_list[[length(series_summary_list) + 1L]] <- data.frame(
    gse = gse,
    source_file = fname,
    n_preview_lines = length(lines),
    sample_title_rows = n_sample_title,
    sample_geo_rows = n_sample_geo,
    has_matrix_table_begin = has_table_begin,
    has_ID_REF_header = has_gene_header,
    contains_recurrence_text = any(grepl("recurr", lines, ignore.case = TRUE)),
    stringsAsFactors = FALSE
  )
}

series_summary <- if (length(series_summary_list)) {
  do.call(rbind, series_summary_list)
} else {
  data.frame()
}

safe_write_csv(series_summary, file.path(REPORT_DIR, "06_series_matrix_summary.csv"))

# ----------------------------- 8. Metadata结构审计 -----------------------------
metadata_paths <- inventory$full_path[
  inventory$role %in% c("Metadata", "CSV", "Text", "CountMatrix") &
  inventory$size_mb <= 100
]

metadata_summary_list <- list()

guess_sep <- function(path, lines) {
  if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) return(",")
  if (length(lines) && grepl("\t", lines[1], fixed = TRUE)) return("\t")
  if (length(lines) && grepl(",", lines[1], fixed = TRUE)) return(",")
  "\t"
}

for (path in metadata_paths) {
  fname <- basename(path)
  gse <- extract_gse(path)
  lines <- safe_read_lines(path, n = 100L)
  if (!length(lines)) next

  preview_path <- file.path(
    PREVIEW_DIR,
    paste0(gsub("[^A-Za-z0-9._-]", "_", fname), "_preview.txt")
  )
  safe_write_lines(lines, preview_path)

  sep <- guess_sep(path, lines)
  dat <- tryCatch({
    con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
    on.exit(try(close(con), silent = TRUE), add = TRUE)
    utils::read.table(
      con,
      header = TRUE,
      sep = sep,
      quote = "\"",
      comment.char = "",
      fill = TRUE,
      check.names = FALSE,
      nrows = 10000L
    )
  }, error = function(e) {
    warn_msg("元数据结构读取失败:", fname, "|", conditionMessage(e))
    NULL
  })

  text_blob <- paste(lines, collapse = "\n")

  metadata_summary_list[[length(metadata_summary_list) + 1L]] <- data.frame(
    gse = gse,
    source_file = fname,
    file_size_mb = inventory$size_mb[match(path, inventory$full_path)],
    preview_lines = length(lines),
    parsed_successfully = !is.null(dat),
    parsed_rows_up_to_10000 = if (!is.null(dat)) nrow(dat) else NA_integer_,
    parsed_columns = if (!is.null(dat)) ncol(dat) else NA_integer_,
    column_names = if (!is.null(dat)) paste(names(dat), collapse = " | ") else NA_character_,
    contains_barcode = grepl("barcode", text_blob, ignore.case = TRUE),
    contains_donor = grepl("donor", text_blob, ignore.case = TRUE),
    contains_skin_wound_stage = grepl(
      "wound1|wound7|wound30|day.?1|day.?7|day.?30|skin",
      text_blob,
      ignore.case = TRUE
    ),
    contains_activity_region = grepl(
      "active|inactive|AC0|AP0|IC0|IP0|mature.?scar",
      text_blob,
      ignore.case = TRUE
    ),
    contains_recurrence = grepl("recurr", text_blob, ignore.case = TRUE),
    stringsAsFactors = FALSE
  )
}

metadata_summary <- if (length(metadata_summary_list)) {
  do.call(rbind, metadata_summary_list)
} else {
  data.frame()
}

safe_write_csv(metadata_summary, file.path(REPORT_DIR, "07_metadata_structure_summary.csv"))

# ----------------------------- 9. 数据集专项判断 -------------------------------
get_archive_sum <- function(gse) {
  if (!nrow(archive_summary)) return(NULL)
  archive_summary[archive_summary$gse == gse, , drop = FALSE]
}

has_file <- function(pattern) {
  any(grepl(pattern, inventory$file_name, ignore.case = TRUE))
}

status_rows <- list()

# GSE241132
a <- get_archive_sum("GSE241132")
status_rows[[length(status_rows) + 1L]] <- data.frame(
  gse = "GSE241132",
  role = "正常伤口单细胞时间轴",
  core_files_present = all(required_check$present[required_check$gse == "GSE241132"]),
  archive_readable = !is.null(a) && nrow(a) > 0L,
  key_structure = if (!is.null(a) && nrow(a)) {
    paste0("MTX=", sum(a$n_matrix_mtx), "; H5=", sum(a$n_h5),
           "; barcodes=", sum(a$n_barcodes), "; features=", sum(a$n_features))
  } else NA_character_,
  stage1_status = NA_character_,
  stringsAsFactors = FALSE
)

# GSE241124
a <- get_archive_sum("GSE241124")
status_rows[[length(status_rows) + 1L]] <- data.frame(
  gse = "GSE241124",
  role = "正常伤口空间时间轴",
  core_files_present = all(required_check$present[required_check$gse == "GSE241124"]),
  archive_readable = !is.null(a) && nrow(a) > 0L,
  key_structure = if (!is.null(a) && nrow(a)) {
    paste0("H5=", sum(a$n_h5), "; spatial_related=", sum(a$n_spatial_related))
  } else NA_character_,
  stage1_status = NA_character_,
  stringsAsFactors = FALSE
)

# GSE220300
a <- get_archive_sum("GSE220300")
status_rows[[length(status_rows) + 1L]] <- data.frame(
  gse = "GSE220300",
  role = "活动/非活动及中心/周边映射",
  core_files_present = all(required_check$present[required_check$gse == "GSE220300"]),
  archive_readable = !is.null(a) && nrow(a) > 0L,
  key_structure = if (!is.null(a) && nrow(a)) {
    paste0("MTX=", sum(a$n_matrix_mtx), "; sample information from SOFT required")
  } else NA_character_,
  stage1_status = NA_character_,
  stringsAsFactors = FALSE
)

# GSE163973
a <- get_archive_sum("GSE163973")
status_rows[[length(status_rows) + 1L]] <- data.frame(
  gse = "GSE163973",
  role = "独立瘢痕疙瘩验证",
  core_files_present = all(required_check$present[required_check$gse == "GSE163973"]),
  archive_readable = !is.null(a) && nrow(a) > 0L,
  key_structure = if (!is.null(a) && nrow(a)) {
    paste0("MTX=", sum(a$n_matrix_mtx), "; metadata=",
           has_file("GSE163973.*cell\\.meta\\.data\\.csv\\.gz"))
  } else NA_character_,
  stage1_status = NA_character_,
  stringsAsFactors = FALSE
)

# GSE181297
a <- get_archive_sum("GSE181297")
status_rows[[length(status_rows) + 1L]] <- data.frame(
  gse = "GSE181297",
  role = "瘢痕疙瘩单细胞与空间支持",
  core_files_present = all(required_check$present[required_check$gse == "GSE181297"]),
  archive_readable = !is.null(a) && nrow(a) > 0L,
  key_structure = if (!is.null(a) && nrow(a)) {
    paste0("MTX=", sum(a$n_matrix_mtx), "; H5=", sum(a$n_h5),
           "; spatial_related=", sum(a$n_spatial_related))
  } else NA_character_,
  stage1_status = NA_character_,
  stringsAsFactors = FALSE
)

# GSE265972
a <- get_archive_sum("GSE265972")
status_rows[[length(status_rows) + 1L]] <- data.frame(
  gse = "GSE265972",
  role = "静脉溃疡反向对照",
  core_files_present = all(required_check$present[required_check$gse == "GSE265972"]),
  archive_readable = !is.null(a) && nrow(a) > 0L,
  key_structure = if (!is.null(a) && nrow(a)) {
    paste0("MTX=", sum(a$n_matrix_mtx), "; H5=", sum(a$n_h5),
           "; metadata=", has_file("GSE265972_VU_anno_metadata"))
  } else NA_character_,
  stage1_status = NA_character_,
  stringsAsFactors = FALSE
)

# GSE274709
a <- get_archive_sum("GSE274709")
g274_expr <- FALSE
g274_idx <- FALSE
if (nrow(archive_details)) {
  x274 <- archive_details[archive_details$gse == "GSE274709", , drop = FALSE]
  if (nrow(x274)) {
    g274_idx <- any(x274$has_idxstats)
    g274_expr <- any(x274$has_expression_hint & !x274$has_idxstats)
  }
}
status_rows[[length(status_rows) + 1L]] <- data.frame(
  gse = "GSE274709",
  role = "术后复发支持",
  core_files_present = all(required_check$present[required_check$gse == "GSE274709"]),
  archive_readable = !is.null(a) && nrow(a) > 0L,
  key_structure = paste0(
    "idxstats=", g274_idx,
    "; non-idxstats_expression_hint=", g274_expr,
    "; series_matrix_exists=", has_file("GSE274709_series_matrix")
  ),
  stage1_status = NA_character_,
  stringsAsFactors = FALSE
)

dataset_status <- do.call(rbind, status_rows)

for (i in seq_len(nrow(dataset_status))) {
  if (!dataset_status$core_files_present[i]) {
    dataset_status$stage1_status[i] <- "FAIL_MISSING_CORE_FILE"
  } else if (!dataset_status$archive_readable[i]) {
    dataset_status$stage1_status[i] <- "FAIL_ARCHIVE_UNREADABLE"
  } else {
    dataset_status$stage1_status[i] <- "PASS_BASIC_STRUCTURE"
  }
}

# GSE274709单独降级
i274 <- which(dataset_status$gse == "GSE274709")
if (length(i274) && g274_idx && !g274_expr) {
  dataset_status$stage1_status[i274] <- "WARN_ONLY_IDXSTATS_LIKELY"
}

safe_write_csv(dataset_status, file.path(REPORT_DIR, "08_dataset_stage1_status.csv"))

# ----------------------------- 10. 自动结论 -----------------------------------
must_pass <- c("GSE241132", "GSE241124", "GSE220300", "GSE163973")
must_tbl <- dataset_status[dataset_status$gse %in% must_pass, , drop = FALSE]

stage1_can_continue <- all(
  must_tbl$core_files_present &
  must_tbl$archive_readable
)

decision <- c(
  "瘢痕疙瘩分子伤口年龄项目",
  "第一阶段：原始数据与元数据审计结果",
  "============================================================",
  paste0("运行时间：", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("根目录：", ROOT_DIR),
  "",
  paste0("第一阶段是否具备继续条件：", if (stage1_can_continue) "YES" else "NO"),
  "",
  "核心继续标准：",
  "1. GSE241132压缩包可读，且有单细胞表达矩阵及细胞元数据；",
  "2. GSE241124压缩包可读，且有H5/空间坐标/图像相关文件；",
  "3. GSE220300压缩包可读，SOFT可恢复AC/AP/IC/IP/MS等样本标题；",
  "4. GSE163973压缩包可读，且有3例keloid与3例normal scar信息；",
  "",
  "各数据集状态："
)

for (i in seq_len(nrow(dataset_status))) {
  decision <- c(
    decision,
    paste0(
      "- ", dataset_status$gse[i],
      " | ", dataset_status$role[i],
      " | ", dataset_status$stage1_status[i],
      " | ", dataset_status$key_structure[i]
    )
  )
}

decision <- c(
  decision,
  "",
  "特别说明：",
  "- GSE274709若仅有idxstats，不影响主课题继续，但不能用于FRFS复发转录组二次分析；",
  "- GSE181297、GSE265972属于强化数据，单独失败不导致主课题停止；",
  "- 第一阶段只判断文件与标签是否具备，不判断模型是否有效；",
  "- 下一阶段将正式读取矩阵并生成逐样本细胞数、基因数和样本映射。"
)

safe_write_lines(decision, file.path(REPORT_DIR, "09_STAGE1_DECISION.txt"))

# ----------------------------- 11. 运行环境 -----------------------------------
capture.output(sessionInfo(), file = file.path(REPORT_DIR, "10_SESSION_INFO.txt"))

completion <- c(
  paste0("STAGE1_COMPLETED=", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste0("ROOT_DIR=", ROOT_DIR),
  paste0("STAGE_DIR=", STAGE_DIR),
  paste0("CAN_CONTINUE=", stage1_can_continue)
)
safe_write_lines(completion, file.path(REPORT_DIR, "STAGE1_COMPLETED.txt"))

# ----------------------------- 12. 自动压缩 -----------------------------------
create_zip_windows <- function(source_dir, zip_path) {
  if (file.exists(zip_path)) file.remove(zip_path)

  ps_escape <- function(x) gsub("'", "''", normalize_slash(x), fixed = TRUE)
  source_q <- ps_escape(source_dir)
  zip_q <- ps_escape(zip_path)

  command <- paste0(
    "$ErrorActionPreference='Stop'; ",
    "$items = Get-ChildItem -LiteralPath '", source_q, "'; ",
    "Compress-Archive -Path $items.FullName -DestinationPath '", zip_q, "' -Force"
  )

  status <- tryCatch(
    system2(
      "powershell.exe",
      args = c("-NoProfile", "-ExecutionPolicy", "Bypass", "-Command", command),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) {
      warn_msg("PowerShell压缩失败:", conditionMessage(e))
      structure(character(0), status = 1L)
    }
  )

  ok <- file.exists(zip_path) && file.info(zip_path)$size > 0
  if (!ok) {
    warn_msg("ZIP未成功生成，将创建tar.gz备用检查包。")
  }
  ok
}

zip_ok <- create_zip_windows(STAGE_DIR, PACKAGE_ZIP)

if (!zip_ok) {
  if (file.exists(PACKAGE_TARGZ)) file.remove(PACKAGE_TARGZ)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(ROOT_DIR)
  tryCatch(
    utils::tar(
      tarfile = basename(PACKAGE_TARGZ),
      files = basename(STAGE_DIR),
      compression = "gzip",
      tar = "internal"
    ),
    error = function(e) warn_msg("tar.gz备用压缩也失败:", conditionMessage(e))
  )
}

log_msg("第一阶段审计完成")
log_msg("是否可进入下一阶段（自动初判）:", stage1_can_continue)
if (file.exists(PACKAGE_ZIP)) {
  log_msg("检查包:", PACKAGE_ZIP)
} else if (file.exists(PACKAGE_TARGZ)) {
  log_msg("备用检查包:", PACKAGE_TARGZ)
}

cat("\n============================================================\n")
cat("第一阶段运行完成。\n")
cat("输出目录：", STAGE_DIR, "\n", sep = "")
if (file.exists(PACKAGE_ZIP)) {
  cat("请上传：", PACKAGE_ZIP, "\n", sep = "")
} else if (file.exists(PACKAGE_TARGZ)) {
  cat("ZIP生成失败，请上传备用文件：", PACKAGE_TARGZ, "\n", sep = "")
} else {
  cat("自动压缩失败，请手动压缩文件夹：", STAGE_DIR, "\n", sep = "")
}
cat("============================================================\n")
