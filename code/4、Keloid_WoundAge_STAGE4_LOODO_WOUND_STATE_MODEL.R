# =============================================================================
# 项目：瘢痕疙瘩“伤口状态持续与修复终止失败”
# 第四阶段：正常伤口状态模型、供者留一验证与固定参考签名
#
# 根目录：
#   C:/Users/33652/Desktop/LZZ/1、yssj
#
# 核心原则：
#   1. 统计单位是供者×阶段×细胞类型，而不是单个细胞；
#   2. 使用GSE241132的3名供者，每次留出1名供者验证；
#   3. 不随机拆分细胞，避免同一供者信息泄漏；
#   4. 先验证4阶段状态签名，再决定是否使用“分子伤口年龄”表述；
#   5. 使用样本内秩评分，便于后续跨队列映射；
#   6. 本阶段不分析任何瘢痕疙瘩数据。
#
# 输入：
#   03_STAGE3_REFERENCE_QC/objects/
#   GSE241132_mainCellType_pseudobulk_counts.rds
#
# 输出：
#   04_STAGE4_WOUND_STATE_MODEL/
#   第四阶段_正常伤口状态模型与供者留一验证检查包.zip
# =============================================================================

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(LANG = "en")

ROOT_DIR <- "C:/Users/33652/Desktop/LZZ/1、yssj"

INPUT_RDS <- file.path(
  ROOT_DIR,
  "03_STAGE3_REFERENCE_QC",
  "objects",
  "GSE241132_mainCellType_pseudobulk_counts.rds"
)

STAGE_DIR <- file.path(ROOT_DIR, "04_STAGE4_WOUND_STATE_MODEL")
REPORT_DIR <- file.path(STAGE_DIR, "report")
FIG_DIR <- file.path(STAGE_DIR, "figures")
OBJECT_DIR <- file.path(STAGE_DIR, "objects")

PACKAGE_ZIP <- file.path(
  ROOT_DIR,
  "第四阶段_正常伤口状态模型与供者留一验证检查包.zip"
)
PACKAGE_TARGZ <- file.path(
  ROOT_DIR,
  "第四阶段_正常伤口状态模型与供者留一验证检查包.tar.gz"
)

if (dir.exists(STAGE_DIR)) {
  unlink(STAGE_DIR, recursive = TRUE, force = TRUE)
}

dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OBJECT_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(REPORT_DIR, "00_stage4_log.txt")
WARNING_FILE <- file.path(REPORT_DIR, "00_stage4_warnings.txt")

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

safe_write_csv <- function(x, path) {
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

safe_write_lines <- function(x, path) {
  tryCatch(
    writeLines(x, path, useBytes = TRUE),
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
  normalizePath(x, winslash = "/", mustWork = FALSE)
}

softmax <- function(x) {
  if (!length(x) || all(!is.finite(x))) {
    return(rep(NA_real_, length(x)))
  }

  x[!is.finite(x)] <- min(x[is.finite(x)])
  z <- x - max(x)
  ez <- exp(z)
  ez / sum(ez)
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)

  if (sum(ok) < 3L) {
    return(NA_real_)
  }

  if (length(unique(x[ok])) < 2L ||
      length(unique(y[ok])) < 2L) {
    return(NA_real_)
  }

  suppressWarnings(
    stats::cor(
      x[ok],
      y[ok],
      method = "spearman"
    )
  )
}

collapse_duplicate_symbols <- function(
  counts,
  gene_symbol
) {
  gene_symbol <- trimws(as.character(gene_symbol))

  valid <- (
    !is.na(gene_symbol) &
      nzchar(gene_symbol)
  )

  counts <- counts[valid, , drop = FALSE]
  gene_symbol <- gene_symbol[valid]

  if (!is.matrix(counts)) {
    counts <- as.matrix(counts)
  }

  collapsed <- rowsum(
    counts,
    group = gene_symbol,
    reorder = FALSE
  )

  storage.mode(collapsed) <- "double"
  collapsed
}

log_cpm <- function(counts) {
  lib <- colSums(counts)

  if (any(lib <= 0)) {
    stop("发现文库总计数≤0的pseudobulk样本。")
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

get_sample_design <- function(
  counts,
  sample_info
) {
  if (is.null(colnames(counts))) {
    stop("pseudobulk计数矩阵缺少列名。")
  }

  sample_info$key <- paste(
    sample_info$sample_id,
    sample_info$main_cell_type,
    sep = "||"
  )

  if (anyDuplicated(sample_info$key)) {
    dup <- unique(
      sample_info$key[duplicated(sample_info$key)]
    )
    stop(
      "sample_info中出现重复键：",
      paste(dup, collapse = ", ")
    )
  }

  idx <- match(
    colnames(counts),
    sample_info$key
  )

  if (anyNA(idx)) {
    stop(
      "计数矩阵列名无法全部匹配sample_info：",
      paste(
        colnames(counts)[is.na(idx)],
        collapse = ", "
      )
    )
  }

  out <- sample_info[idx, , drop = FALSE]

  if (!identical(
    colnames(counts),
    out$key
  )) {
    stop("pseudobulk列顺序与sample_info匹配失败。")
  }

  out
}

stage_order <- c(
  "Skin",
  "Wound1",
  "Wound7",
  "Wound30"
)

stage_position <- c(
  Skin = 0,
  Wound1 = 1,
  Wound7 = 2,
  Wound30 = 3
)

stage_label_short <- c(
  Skin = "Skin",
  Wound1 = "D1",
  Wound7 = "D7",
  Wound30 = "D30"
)

exclude_gene <- function(gene) {
  g <- toupper(as.character(gene))

  grepl(
    "^MT-|^RPS[0-9]|^RPL[0-9]|^HB[ABDEGMQZ][0-9A-Z]*$",
    g
  )
}

filter_training_genes <- function(
  counts_train,
  min_cpm = 1,
  min_samples = 4L
) {
  lib <- colSums(counts_train)

  cpm <- sweep(
    counts_train,
    2,
    pmax(lib, 1),
    "/"
  ) * 1e6

  keep <- (
    rowSums(cpm >= min_cpm) >= min_samples &
      rowSums(counts_train >= 10) >= 3L &
      !exclude_gene(rownames(counts_train))
  )

  rownames(counts_train)[keep]
}

within_donor_stage_effect <- function(
  expr_train,
  design_train,
  gene_names
) {
  expr_train <- expr_train[
    gene_names,
    ,
    drop = FALSE
  ]

  donors <- unique(design_train$donor)

  effects <- vector(
    "list",
    length(stage_order)
  )

  names(effects) <- stage_order

  for (st in stage_order) {
    eff_by_donor <- lapply(
      donors,
      function(d) {
        ids <- which(
          design_train$donor == d
        )

        d_design <- design_train[
          ids,
          ,
          drop = FALSE
        ]

        d_expr <- expr_train[
          ,
          ids,
          drop = FALSE
        ]

        if (!all(stage_order %in% d_design$condition)) {
          stop(
            "训练供者",
            d,
            "缺少四个阶段。"
          )
        }

        target_col <- which(
          d_design$condition == st
        )

        other_cols <- which(
          d_design$condition != st
        )

        d_expr[, target_col] -
          rowMeans(
            d_expr[
              ,
              other_cols,
              drop = FALSE
            ]
          )
      }
    )

    effects[[st]] <- do.call(
      cbind,
      eff_by_donor
    )

    colnames(effects[[st]]) <- donors
  }

  effects
}

derive_stage_signatures <- function(
  counts_train,
  design_train,
  top_n = 75L,
  min_effect = 0.35,
  min_genes = 20L
) {
  train_genes <- filter_training_genes(
    counts_train,
    min_cpm = 1,
    min_samples = max(
      3L,
      floor(ncol(counts_train) / 2)
    )
  )

  if (length(train_genes) < 500L) {
    stop(
      "训练集中可用基因不足500个，实际为：",
      length(train_genes)
    )
  }

  expr_train <- log_cpm(
    counts_train[
      train_genes,
      ,
      drop = FALSE
    ]
  )

  effects <- within_donor_stage_effect(
    expr_train,
    design_train,
    train_genes
  )

  signature_list <- vector(
    "list",
    length(stage_order)
  )

  names(signature_list) <- stage_order

  signature_details <- list()

  for (st in stage_order) {
    eff <- effects[[st]]
    mean_eff <- rowMeans(eff)

    direction_consistent <- apply(
      eff,
      1,
      function(x) all(x > 0)
    )

    candidate <- (
      mean_eff >= min_effect &
        direction_consistent
    )

    selected <- names(
      sort(
        mean_eff[candidate],
        decreasing = TRUE
      )
    )

    if (length(selected) < min_genes) {
      candidate <- (
        mean_eff > 0 &
          direction_consistent
      )

      selected <- names(
        sort(
          mean_eff[candidate],
          decreasing = TRUE
        )
      )
    }

    selected <- head(
      selected,
      top_n
    )

    signature_list[[st]] <- selected

    signature_details[[st]] <- data.frame(
      stage = st,
      gene = names(mean_eff),
      mean_stage_specific_effect = as.numeric(mean_eff),
      direction_consistent = direction_consistent,
      selected = names(mean_eff) %in% selected,
      stringsAsFactors = FALSE
    )
  }

  list(
    signatures = signature_list,
    details = do.call(
      rbind,
      signature_details
    ),
    training_genes = train_genes
  )
}

rank_signature_scores <- function(
  expr_sample,
  signatures
) {
  if (is.matrix(expr_sample)) {
    if (ncol(expr_sample) != 1L) {
      stop("rank_signature_scores只接受单个样本。")
    }

    expr_sample <- expr_sample[, 1]
  }

  gene_names <- names(expr_sample)

  if (is.null(gene_names)) {
    stop("表达向量缺少基因名。")
  }

  ranks <- rank(
    expr_sample,
    ties.method = "average",
    na.last = "keep"
  )

  ranks <- ranks / sum(is.finite(ranks))
  names(ranks) <- gene_names

  vapply(
    stage_order,
    function(st) {
      genes <- intersect(
        signatures[[st]],
        names(ranks)
      )

      if (length(genes) < 5L) {
        return(NA_real_)
      }

      mean(
        ranks[genes],
        na.rm = TRUE
      )
    },
    numeric(1)
  )
}

evaluate_cell_type_loocv <- function(
  counts_ct,
  design_ct,
  cell_type
) {
  donors <- unique(design_ct$donor)

  if (length(donors) != 3L) {
    stop(
      cell_type,
      "不是3名供者，实际为：",
      length(donors)
    )
  }

  pred_rows <- list()
  signature_rows <- list()

  for (held_out in donors) {
    train_idx <- which(
      design_ct$donor != held_out
    )

    test_idx <- which(
      design_ct$donor == held_out
    )

    train_design <- design_ct[
      train_idx,
      ,
      drop = FALSE
    ]

    test_design <- design_ct[
      test_idx,
      ,
      drop = FALSE
    ]

    train_counts <- counts_ct[
      ,
      train_idx,
      drop = FALSE
    ]

    test_counts <- counts_ct[
      ,
      test_idx,
      drop = FALSE
    ]

    model <- derive_stage_signatures(
      train_counts,
      train_design,
      top_n = 75L,
      min_effect = 0.35,
      min_genes = 20L
    )

    fold_sig <- do.call(
      rbind,
      lapply(
        stage_order,
        function(st) {
          data.frame(
            cell_type = cell_type,
            held_out_donor = held_out,
            stage = st,
            gene = model$signatures[[st]],
            stringsAsFactors = FALSE
          )
        }
      )
    )

    signature_rows[[length(signature_rows) + 1L]] <- fold_sig

    test_genes <- intersect(
      model$training_genes,
      rownames(test_counts)
    )

    test_expr <- log_cpm(
      test_counts[
        test_genes,
        ,
        drop = FALSE
      ]
    )

    for (j in seq_len(ncol(test_expr))) {
      scores <- rank_signature_scores(
        test_expr[, j],
        model$signatures
      )

      pred_stage <- names(
        which.max(scores)
      )

      score_sd <- stats::sd(scores)

      if (!is.finite(score_sd) ||
          score_sd == 0) {
        probs <- rep(
          1 / length(stage_order),
          length(stage_order)
        )
        names(probs) <- stage_order
      } else {
        z <- (
          scores - mean(scores)
        ) / score_sd

        probs <- softmax(z * 2)
        names(probs) <- stage_order
      }

      predicted_position <- sum(
        probs * stage_position[
          names(probs)
        ]
      )

      true_stage <- as.character(
        test_design$condition[j]
      )

      pred_rows[[length(pred_rows) + 1L]] <- data.frame(
        cell_type = cell_type,
        held_out_donor = held_out,
        sample_id = test_design$sample_id[j],
        true_stage = true_stage,
        true_position = unname(
          stage_position[true_stage]
        ),
        predicted_stage = pred_stage,
        predicted_position = predicted_position,
        score_Skin = scores["Skin"],
        score_Wound1 = scores["Wound1"],
        score_Wound7 = scores["Wound7"],
        score_Wound30 = scores["Wound30"],
        probability_Skin = probs["Skin"],
        probability_Wound1 = probs["Wound1"],
        probability_Wound7 = probs["Wound7"],
        probability_Wound30 = probs["Wound30"],
        exact_correct = pred_stage == true_stage,
        ordinal_error = abs(
          unname(stage_position[pred_stage]) -
            unname(stage_position[true_stage])
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  predictions <- do.call(
    rbind,
    pred_rows
  )

  donor_metrics <- do.call(
    rbind,
    lapply(
      split(
        predictions,
        predictions$held_out_donor
      ),
      function(x) {
        data.frame(
          cell_type = cell_type,
          held_out_donor = x$held_out_donor[1],
          exact_accuracy = mean(x$exact_correct),
          mean_ordinal_error = mean(x$ordinal_error),
          spearman_true_vs_predicted_position =
            safe_spearman(
              x$true_position,
              x$predicted_position
            ),
          stringsAsFactors = FALSE
        )
      }
    )
  )

  overall_metrics <- data.frame(
    cell_type = cell_type,
    n_predictions = nrow(predictions),
    exact_accuracy = mean(
      predictions$exact_correct
    ),
    adjacent_accuracy = mean(
      predictions$ordinal_error <= 1
    ),
    mean_ordinal_error = mean(
      predictions$ordinal_error
    ),
    pooled_spearman = safe_spearman(
      predictions$true_position,
      predictions$predicted_position
    ),
    donors_with_positive_rho = sum(
      donor_metrics$spearman_true_vs_predicted_position > 0,
      na.rm = TRUE
    ),
    donors_with_rho_ge_0_5 = sum(
      donor_metrics$spearman_true_vs_predicted_position >= 0.5,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )

  list(
    predictions = predictions,
    donor_metrics = donor_metrics,
    overall_metrics = overall_metrics,
    fold_signatures = do.call(
      rbind,
      signature_rows
    )
  )
}

derive_final_reference <- function(
  counts_ct,
  design_ct,
  cell_type
) {
  model <- derive_stage_signatures(
    counts_ct,
    design_ct,
    top_n = 100L,
    min_effect = 0.35,
    min_genes = 20L
  )

  final_signatures <- do.call(
    rbind,
    lapply(
      stage_order,
      function(st) {
        details <- model$details[
          model$details$stage == st &
            model$details$selected,
          ,
          drop = FALSE
        ]

        details$cell_type <- cell_type

        details[
          ,
          c(
            "cell_type",
            "stage",
            "gene",
            "mean_stage_specific_effect",
            "direction_consistent"
          )
        ]
      }
    )
  )

  list(
    signatures = model$signatures,
    signature_table = final_signatures,
    training_genes = model$training_genes
  )
}

plot_loocv_predictions <- function(
  predictions,
  path
) {
  cell_types <- unique(
    predictions$cell_type
  )

  pdf(
    path,
    width = 10,
    height = 7,
    onefile = TRUE
  )

  for (ct in cell_types) {
    x <- predictions[
      predictions$cell_type == ct,
      ,
      drop = FALSE
    ]

    x <- x[
      order(
        x$held_out_donor,
        x$true_position
      ),
      ,
      drop = FALSE
    ]

    plot(
      jitter(
        x$true_position,
        amount = 0.04
      ),
      x$predicted_position,
      pch = 16,
      xlim = c(-0.2, 3.2),
      ylim = c(-0.2, 3.2),
      xaxt = "n",
      xlab = "True wound-state position",
      ylab = "Predicted ordinal wound-state position",
      main = paste0(
        ct,
        ": leave-one-donor-out validation"
      )
    )

    axis(
      1,
      at = 0:3,
      labels = stage_label_short[
        stage_order
      ]
    )

    abline(
      a = 0,
      b = 1,
      lty = 2
    )

    text(
      jitter(
        x$true_position,
        amount = 0.04
      ),
      x$predicted_position,
      labels = x$held_out_donor,
      pos = 3,
      cex = 0.7
    )
  }

  dev.off()
}

plot_confusion_matrices <- function(
  predictions,
  path
) {
  cell_types <- unique(
    predictions$cell_type
  )

  pdf(
    path,
    width = 7,
    height = 6,
    onefile = TRUE
  )

  for (ct in cell_types) {
    x <- predictions[
      predictions$cell_type == ct,
      ,
      drop = FALSE
    ]

    tab <- table(
      factor(
        x$true_stage,
        levels = stage_order
      ),
      factor(
        x$predicted_stage,
        levels = stage_order
      )
    )

    image(
      x = seq_len(ncol(tab)),
      y = seq_len(nrow(tab)),
      z = t(tab),
      axes = FALSE,
      xlab = "Predicted stage",
      ylab = "True stage",
      main = paste0(
        ct,
        ": LOODO confusion matrix"
      )
    )

    axis(
      1,
      at = seq_along(stage_order),
      labels = stage_label_short[
        stage_order
      ]
    )

    axis(
      2,
      at = seq_along(stage_order),
      labels = stage_label_short[
        stage_order
      ]
    )

    for (i in seq_len(nrow(tab))) {
      for (j in seq_len(ncol(tab))) {
        text(
          j,
          i,
          labels = tab[i, j]
        )
      }
    }

    box()
  }

  dev.off()
}

plot_fibroblast_signature_heatmap <- function(
  counts_ct,
  design_ct,
  final_signature_table,
  path
) {
  sig_genes <- unique(
    final_signature_table$gene
  )

  sig_genes <- intersect(
    sig_genes,
    rownames(counts_ct)
  )

  if (length(sig_genes) < 10L) {
    warn_msg(
      "成纤维细胞最终签名基因不足10个，跳过热图。"
    )
    return(invisible(NULL))
  }

  expr <- log_cpm(
    counts_ct[
      sig_genes,
      ,
      drop = FALSE
    ]
  )

  # 供者内中心化，展示时间变化而非供者基线差异。
  for (d in unique(design_ct$donor)) {
    idx <- which(
      design_ct$donor == d
    )

    expr[, idx] <- sweep(
      expr[, idx, drop = FALSE],
      1,
      rowMeans(
        expr[, idx, drop = FALSE]
      ),
      "-"
    )
  }

  gene_sd <- apply(
    expr,
    1,
    stats::sd
  )

  keep <- is.finite(gene_sd) &
    gene_sd > 0

  expr <- expr[keep, , drop = FALSE]

  if (nrow(expr) > 120L) {
    expr <- expr[
      order(
        apply(expr, 1, stats::sd),
        decreasing = TRUE
      )[seq_len(120L)],
      ,
      drop = FALSE
    ]
  }

  expr_z <- t(
    scale(
      t(expr)
    )
  )

  sample_ord <- order(
    design_ct$donor,
    match(
      design_ct$condition,
      stage_order
    )
  )

  expr_z <- expr_z[
    ,
    sample_ord,
    drop = FALSE
  ]

  sample_labels <- paste(
    design_ct$donor[sample_ord],
    stage_label_short[
      design_ct$condition[sample_ord]
    ],
    sep = "_"
  )

  pdf(
    path,
    width = 10,
    height = 12
  )

  par(
    mar = c(9, 7, 3, 2)
  )

  image(
    x = seq_len(ncol(expr_z)),
    y = seq_len(nrow(expr_z)),
    z = t(expr_z),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = "Fibroblast fixed wound-state signatures"
  )

  axis(
    1,
    at = seq_len(ncol(expr_z)),
    labels = sample_labels,
    las = 2,
    cex.axis = 0.75
  )

  axis(
    2,
    at = seq_len(nrow(expr_z)),
    labels = rownames(expr_z),
    las = 2,
    cex.axis = 0.42
  )

  box()
  dev.off()
}

# ----------------------------- 读取输入 ----------------------------------------
log_msg("第四阶段开始")

if (!file.exists(INPUT_RDS)) {
  stop(
    "未找到第三阶段pseudobulk对象：",
    INPUT_RDS
  )
}

obj <- readRDS(INPUT_RDS)

required_obj <- c(
  "counts",
  "sample_info",
  "gene_symbol"
)

missing_obj <- setdiff(
  required_obj,
  names(obj)
)

if (length(missing_obj)) {
  stop(
    "第三阶段RDS缺少：",
    paste(
      missing_obj,
      collapse = ", "
    )
  )
}

counts_raw <- obj$counts
sample_info <- obj$sample_info
gene_symbol <- obj$gene_symbol

if (nrow(counts_raw) != length(gene_symbol)) {
  stop(
    "gene_symbol长度与计数矩阵行数不一致。"
  )
}

design <- get_sample_design(
  counts_raw,
  sample_info
)

counts <- collapse_duplicate_symbols(
  counts_raw,
  gene_symbol
)

rm(counts_raw)
gc(verbose = FALSE)

design$donor <- as.character(
  design$donor
)

design$condition <- as.character(
  design$condition
)

design$main_cell_type <- as.character(
  design$main_cell_type
)

if (!all(
  design$condition %in% stage_order
)) {
  stop(
    "发现未知阶段：",
    paste(
      setdiff(
        unique(design$condition),
        stage_order
      ),
      collapse = ", "
    )
  )
}

# 只分析第三阶段确认全部12组均至少50细胞的5种主细胞类型
primary_cell_types <- c(
  "Fibroblast",
  "Keratinocyte",
  "Myeloid",
  "Endothelial",
  "Lymphoid"
)

available_cell_types <- intersect(
  primary_cell_types,
  unique(design$main_cell_type)
)

if (!"Fibroblast" %in% available_cell_types) {
  stop("pseudobulk中不存在Fibroblast。")
}

log_msg(
  "进入LOODO的细胞类型：",
  paste(
    available_cell_types,
    collapse = ", "
  )
)

# ----------------------------- 模型验证 ----------------------------------------
prediction_list <- list()
donor_metric_list <- list()
overall_metric_list <- list()
fold_signature_list <- list()
final_signature_tables <- list()
final_reference_objects <- list()

for (ct in available_cell_types) {
  log_msg("供者留一验证：", ct)

  idx <- which(
    design$main_cell_type == ct
  )

  counts_ct <- counts[
    ,
    idx,
    drop = FALSE
  ]

  design_ct <- design[
    idx,
    ,
    drop = FALSE
  ]

  ord <- order(
    design_ct$donor,
    match(
      design_ct$condition,
      stage_order
    )
  )

  counts_ct <- counts_ct[
    ,
    ord,
    drop = FALSE
  ]

  design_ct <- design_ct[
    ord,
    ,
    drop = FALSE
  ]

  per_donor_stage <- table(
    design_ct$donor,
    factor(
      design_ct$condition,
      levels = stage_order
    )
  )

  if (!all(per_donor_stage == 1L)) {
    stop(
      ct,
      "不是每名供者每个阶段恰好1个pseudobulk样本。"
    )
  }

  cv <- evaluate_cell_type_loocv(
    counts_ct,
    design_ct,
    ct
  )

  prediction_list[[ct]] <- cv$predictions
  donor_metric_list[[ct]] <- cv$donor_metrics
  overall_metric_list[[ct]] <- cv$overall_metrics
  fold_signature_list[[ct]] <- cv$fold_signatures

  final_ref <- derive_final_reference(
    counts_ct,
    design_ct,
    ct
  )

  final_signature_tables[[ct]] <-
    final_ref$signature_table

  final_reference_objects[[ct]] <- list(
    signatures = final_ref$signatures,
    training_genes = final_ref$training_genes,
    stage_order = stage_order,
    stage_position = stage_position
  )

  if (ct == "Fibroblast") {
    plot_fibroblast_signature_heatmap(
      counts_ct,
      design_ct,
      final_ref$signature_table,
      file.path(
        FIG_DIR,
        "Figure_fibroblast_fixed_signature_heatmap.pdf"
      )
    )
  }
}

predictions <- do.call(
  rbind,
  prediction_list
)

donor_metrics <- do.call(
  rbind,
  donor_metric_list
)

overall_metrics <- do.call(
  rbind,
  overall_metric_list
)

fold_signatures <- do.call(
  rbind,
  fold_signature_list
)

final_signatures <- do.call(
  rbind,
  final_signature_tables
)

rownames(predictions) <- NULL
rownames(donor_metrics) <- NULL
rownames(overall_metrics) <- NULL
rownames(fold_signatures) <- NULL
rownames(final_signatures) <- NULL

safe_write_csv(
  predictions,
  file.path(
    REPORT_DIR,
    "01_LOODO_predictions.csv"
  )
)

safe_write_csv(
  donor_metrics,
  file.path(
    REPORT_DIR,
    "02_LOODO_donor_metrics.csv"
  )
)

safe_write_csv(
  overall_metrics,
  file.path(
    REPORT_DIR,
    "03_LOODO_overall_metrics.csv"
  )
)

safe_write_csv(
  fold_signatures,
  file.path(
    REPORT_DIR,
    "04_fold_specific_signatures.csv"
  )
)

safe_write_csv(
  final_signatures,
  file.path(
    REPORT_DIR,
    "05_fixed_final_stage_signatures.csv"
  )
)

# 签名稳定性：每个基因在3次LOODO中被选择几次
signature_stability <- aggregate(
  x = list(
    n_folds_selected = rep(
      1L,
      nrow(fold_signatures)
    )
  ),
  by = fold_signatures[
    ,
    c(
      "cell_type",
      "stage",
      "gene"
    )
  ],
  FUN = sum
)

signature_stability$stable_all_3_folds <-
  signature_stability$n_folds_selected == 3L

safe_write_csv(
  signature_stability,
  file.path(
    REPORT_DIR,
    "06_signature_stability_across_folds.csv"
  )
)

# 混淆矩阵
confusion_rows <- list()

for (ct in unique(predictions$cell_type)) {
  x <- predictions[
    predictions$cell_type == ct,
    ,
    drop = FALSE
  ]

  tab <- as.data.frame(
    table(
      true_stage = factor(
        x$true_stage,
        levels = stage_order
      ),
      predicted_stage = factor(
        x$predicted_stage,
        levels = stage_order
      )
    ),
    stringsAsFactors = FALSE
  )

  tab$cell_type <- ct

  confusion_rows[[ct]] <- tab[
    ,
    c(
      "cell_type",
      "true_stage",
      "predicted_stage",
      "Freq"
    )
  ]
}

confusion_table <- do.call(
  rbind,
  confusion_rows
)

rownames(confusion_table) <- NULL

safe_write_csv(
  confusion_table,
  file.path(
    REPORT_DIR,
    "07_LOODO_confusion_matrices.csv"
  )
)

# ----------------------------- 固定模型对象 ------------------------------------
saveRDS(
  list(
    cell_type_models = final_reference_objects,
    fixed_signature_table = final_signatures,
    signature_stability = signature_stability,
    validation_predictions = predictions,
    validation_metrics = overall_metrics,
    stage_order = stage_order,
    stage_position = stage_position,
    scoring_method = paste(
      "Within-sample percentile-rank mean score;",
      "four stage-specific positive gene sets;",
      "LOODO by donor."
    )
  ),
  file = file.path(
    OBJECT_DIR,
    "GSE241132_fixed_wound_state_reference.rds"
  ),
  compress = "gzip"
)

# ----------------------------- 图形 --------------------------------------------
plot_loocv_predictions(
  predictions,
  file.path(
    FIG_DIR,
    "Figure_LOODO_ordinal_predictions.pdf"
  )
)

plot_confusion_matrices(
  predictions,
  file.path(
    FIG_DIR,
    "Figure_LOODO_confusion_matrices.pdf"
  )
)

# 性能汇总图
pdf(
  file.path(
    FIG_DIR,
    "Figure_model_performance_summary.pdf"
  ),
  width = 10,
  height = 6
)

op <- par(no.readonly = TRUE)
par(
  mfrow = c(1, 2),
  mar = c(8, 4, 3, 1)
)

barplot(
  overall_metrics$exact_accuracy,
  names.arg = overall_metrics$cell_type,
  las = 2,
  ylim = c(0, 1),
  ylab = "Exact stage accuracy",
  main = "Leave-one-donor-out exact accuracy"
)

abline(
  h = c(0.50, 0.67),
  lty = c(3, 2)
)

barplot(
  overall_metrics$pooled_spearman,
  names.arg = overall_metrics$cell_type,
  las = 2,
  ylim = c(-1, 1),
  ylab = "Spearman correlation",
  main = "True vs predicted ordinal position"
)

abline(
  h = c(0, 0.70),
  lty = c(3, 2)
)

par(op)
dev.off()

# ----------------------------- 决策闸门 ----------------------------------------
fib_metric <- overall_metrics[
  overall_metrics$cell_type == "Fibroblast",
  ,
  drop = FALSE
]

if (nrow(fib_metric) != 1L) {
  stop("无法唯一获得Fibroblast验证指标。")
}

fib_strong_pass <- (
  fib_metric$exact_accuracy >= 0.67 &
    fib_metric$pooled_spearman >= 0.70 &
    fib_metric$adjacent_accuracy >= 0.90 &
    fib_metric$donors_with_positive_rho >= 2L
)

fib_moderate_pass <- (
  fib_metric$exact_accuracy >= 0.50 &
    fib_metric$pooled_spearman >= 0.50 &
    fib_metric$mean_ordinal_error <= 0.75
)

secondary_pass <- with(
  overall_metrics[
    overall_metrics$cell_type != "Fibroblast",
    ,
    drop = FALSE
  ],
  exact_accuracy >= 0.50 &
    pooled_spearman >= 0.50 &
    adjacent_accuracy >= 0.80
)

n_secondary_pass <- sum(
  secondary_pass,
  na.rm = TRUE
)

if (fib_strong_pass &&
    n_secondary_pass >= 2L) {
  project_gate <- "PASS_STRONG"
  primary_term <- "ordinal molecular wound-state position"
  next_step <- paste(
    "进入瘢痕疙瘩映射阶段；",
    "可保留“分子伤口状态/相对伤口年龄”概念，",
    "但仍不声称精确日龄。"
  )
} else if (fib_moderate_pass) {
  project_gate <- "PASS_STAGE_SIGNATURE_ONLY"
  primary_term <- "wound-state persistence and failure of resolution"
  next_step <- paste(
    "进入瘢痕疙瘩映射阶段；",
    "不使用精确“分子伤口年龄”，",
    "仅使用Skin/D1/D7/D30-like状态及早期持续评分。"
  )
} else {
  project_gate <- "FAIL_TIME_REFERENCE"
  primary_term <- "stage-associated wound signatures"
  next_step <- paste(
    "停止分子伤口年龄主线；",
    "改为较保守的正常伤口阶段签名与瘢痕疙瘩偏轨分析。"
  )
}

stable_signature_summary <- aggregate(
  x = list(
    n_stable_genes = signature_stability$stable_all_3_folds
  ),
  by = signature_stability[
    ,
    c(
      "cell_type",
      "stage"
    )
  ],
  FUN = sum
)

safe_write_csv(
  stable_signature_summary,
  file.path(
    REPORT_DIR,
    "08_stable_signature_counts.csv"
  )
)

decision_lines <- c(
  "第四阶段：正常伤口状态模型与供者留一验证结论",
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
  paste0(
    "建议正文主术语：",
    primary_term
  ),
  "",
  "Fibroblast验证结果：",
  paste0(
    "- exact accuracy = ",
    round(
      fib_metric$exact_accuracy,
      3
    )
  ),
  paste0(
    "- adjacent accuracy = ",
    round(
      fib_metric$adjacent_accuracy,
      3
    )
  ),
  paste0(
    "- mean ordinal error = ",
    round(
      fib_metric$mean_ordinal_error,
      3
    )
  ),
  paste0(
    "- pooled Spearman = ",
    round(
      fib_metric$pooled_spearman,
      3
    )
  ),
  paste0(
    "- donors with positive rho = ",
    fib_metric$donors_with_positive_rho
  ),
  "",
  paste0(
    "通过预设标准的次要细胞类型数量：",
    n_secondary_pass
  ),
  "",
  "下一步：",
  next_step,
  "",
  "解释边界：",
  "- 只有3名正常伤口供者，因此不能宣称通用临床预测模型。",
  "- 四个离散时间点不能支持精确的连续日龄估计。",
  "- 即使PASS_STRONG，也优先写作ordinal wound-state position。",
  "- 活动度、复发和空间结果必须在后续独立队列中验证。"
)

safe_write_lines(
  decision_lines,
  file.path(
    REPORT_DIR,
    "09_STAGE4_DECISION.txt"
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
      "STAGE4_COMPLETED=",
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
      "FIBROBLAST_STRONG_PASS=",
      fib_strong_pass
    ),
    paste0(
      "FIBROBLAST_MODERATE_PASS=",
      fib_moderate_pass
    ),
    paste0(
      "N_SECONDARY_PASS=",
      n_secondary_pass
    )
  ),
  file.path(
    REPORT_DIR,
    "STAGE4_COMPLETED.txt"
  )
)

# ----------------------------- 自动打包 ----------------------------------------
package_dir <- file.path(
  STAGE_DIR,
  "_package_for_review"
)

if (dir.exists(package_dir)) {
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
  if (file.exists(zip_path)) {
    file.remove(zip_path)
  }

  esc <- function(x) {
    gsub(
      "'",
      "''",
      normalize_slash(x),
      fixed = TRUE
    )
  }

  cmd <- paste0(
    "$ErrorActionPreference='Stop'; ",
    "$items=Get-ChildItem -LiteralPath '",
    esc(source_dir),
    "'; ",
    "Compress-Archive -Path $items.FullName ",
    "-DestinationPath '",
    esc(zip_path),
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
        cmd
      ),
      stdout = TRUE,
      stderr = TRUE
    ),
    silent = TRUE
  )

  file.exists(zip_path) &&
    file.info(zip_path)$size > 0
}

zip_ok <- create_zip_windows(
  package_dir,
  PACKAGE_ZIP
)

if (!zip_ok) {
  if (file.exists(PACKAGE_TARGZ)) {
    file.remove(PACKAGE_TARGZ)
  }

  old_wd <- getwd()
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

  setwd(old_wd)
}

unlink(
  package_dir,
  recursive = TRUE,
  force = TRUE
)

log_msg("第四阶段完成")

cat("\n============================================================\n")
cat("第四阶段运行完成。\n")
cat(
  "固定参考模型保存在：",
  file.path(
    OBJECT_DIR,
    "GSE241132_fixed_wound_state_reference.rds"
  ),
  "\n",
  sep = ""
)

if (file.exists(PACKAGE_ZIP)) {
  cat(
    "请上传：",
    PACKAGE_ZIP,
    "\n",
    sep = ""
  )
} else if (file.exists(PACKAGE_TARGZ)) {
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

cat("============================================================\n")
