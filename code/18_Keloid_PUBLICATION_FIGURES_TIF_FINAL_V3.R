# ============================================================================
# Keloid study - FINAL ALL-IN-ONE FIGURE PIPELINE
# Version: FINAL FIXED | 2026-06-19
# ============================================================================
# Fixed issue:
#   The previous all-in-one file stopped after individual-panel generation because
#   the panels-only section incorrectly required composite TIFFs before the
#   composite stitching section had run. This version allows Step 1 to finish and
#   then executes Step 2.
# ============================================================================


# ============================================================================
# PART 1/2: INDIVIDUAL PANEL GENERATION
# ============================================================================

# ============================================================================
# Keloid persistent late-remodeling study
# INDIVIDUAL-PANEL publication figure pipeline
# Version: FINAL SINGLE PANELS + FIGURE2B DISEASE PROJECTION | FIXED COMPOSITE HANDOFF | 2026-06-19
#
# PURPOSE
#   1) Generate the complete publication figure set as TIFF only.
#   2) Assemble related panels into non-uniform, story-driven composite figures.
#   3) Export panel-level source data for every figure.
#   4) Export main and supplementary tables used by the manuscript.
#
# EXPECTED INPUT
#   The script uses the fixed Windows project root below and also performs
#   a validated fallback search when the folder has been moved.
#   Required stage outputs must be extracted; ZIP files alone are not used.
#
# OUTPUT
#   15_PUBLICATION_FIGURES_FINAL/
#     main_tif/          [empty unless SAVE_COMPOSITE_TIFFS <- TRUE]
#     supplementary_tif/ [empty unless SAVE_COMPOSITE_TIFFS <- TRUE]Supplementary_Figure_S1.tif ... S8.tif
#     individual_panel_tif/All main and supplementary panels as TIFF
#     source_data/       Panel-level CSV files + index
#     tables/            Main and supplementary table CSV files
#     logs/              Input manifest, sessionInfo, completion flag
#
# IMPORTANT
#   - Figures are TIFF only, white background, LZW compressed, 600 dpi.
#   - Composite widths are set to approximately 7.2 inches (full journal width),
#     with figure-specific heights to avoid oversized files and unreadably small text.
#   - No cell-level P values are used as patient-level evidence.
#   - Figure 6 uses association language, not causal mechanism language.
#   - GSE181297 is kept supplementary because its cross-modal evidence is mixed.
#   - No workflow/flowchart panel is included in the main figures.
#   - Main figures emphasize the prespecified positive evidence chain; competing
#     hypotheses and modality-dependent findings remain transparently supplementary.
#   - Version PANELS-ONLY 1.1 fixes single-panel problems directly:
#     clean confusion matrices, donor-split spatial sections, and resource-split
#     non-redundant GSEA plots. Composite export remains disabled by default.
#   - Version 1.2 additionally redraws Figure5_F as a legend-free, unclipped
#     wound-state failure plane with expanded axes and manual group labels.
# ============================================================================

options(stringsAsFactors = FALSE, scipen = 999)
set.seed(20260618)

# ------------------------------- USER SETTINGS -------------------------------
# Fixed project path used throughout Stages 1-12.
USER_PROJECT_ROOT <- "C:/Users/33652/Desktop/LZZ/1、yssj"

DPI <- 600
BASE_FONT <- if (.Platform$OS.type == "windows") "Arial" else "sans"
AUTO_INSTALL_PACKAGES <- TRUE
SAVE_PANEL_TIFFS <- TRUE
SAVE_COMPOSITE_TIFFS <- FALSE
PANEL_BASE_WIDTH <- 6.0
PANEL_BASE_HEIGHT <- 4.8
OUTPUT_FOLDER_NAME <- "__TEMP_KELOID_PANEL_ENGINE__"

# The script first validates USER_PROJECT_ROOT. If the project was moved, it
# evaluates the script directory, current working directory and their parents.
get_current_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg)) {
    return(sub("^--file=", "", file_arg[[1]]))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    context_path <- tryCatch(
      rstudioapi::getActiveDocumentContext()$path,
      error = function(e) ""
    )
    if (nzchar(context_path)) {
      return(context_path)
    }
  }

  ""
}

root_anchor_files <- c(
  "01_LOODO_predictions.csv",
  "09_GSE163973_patient_wound_state_scores.csv",
  "04_patient_level_scores_after_K3_merge.csv",
  "09_patient_lesion_and_scar_scores.csv",
  "09_spatial_section_wound_state_scores.csv",
  "10_GSE265972_wound_state_scores.csv",
  "02_core_patient_level_data.csv",
  "08_consensus_persistent_late_remodeling_genes.csv"
)

# Fast project-root validation.
#
# V2.0 recursively scanned multiple parent directories before accepting the
# already specified fixed project path. On Windows this could reach C:/Users or
# even the drive root and appear to freeze at vapply(). V2.1 never recursively
# scans parent directories. The fixed project path is accepted immediately when
# it exists. A limited fallback is used only if the fixed path has moved.

stage_directory_pattern <- paste0(
  "(?i)(",
  "^0[1-9][A-Z]?_STAGE|",
  "^1[0-2]_STAGE|",
  "STAGE[0-9]|",
  "阶段|",
  "stage[0-9]+_|",
  "_review|",
  "_inspect|",
  "_reinspect|",
  "_audit",
  ")"
)

fast_anchor_score <- function(path) {
  if (
    is.na(path) ||
      !nzchar(path) ||
      !dir.exists(path)
  ) {
    return(-Inf)
  }

  path <- normalizePath(
    path,
    winslash = "/",
    mustWork = TRUE
  )

  root_files <- list.files(
    path,
    recursive = FALSE,
    full.names = TRUE,
    include.dirs = FALSE,
    all.files = FALSE
  )

  first_level_dirs <- list.dirs(
    path,
    recursive = FALSE,
    full.names = TRUE
  )

  stage_dirs <- first_level_dirs[
    grepl(
      stage_directory_pattern,
      basename(first_level_dirs),
      perl = TRUE
    )
  ]

  stage_files <- unlist(
    lapply(
      stage_dirs,
      function(d) {
        list.files(
          d,
          recursive = TRUE,
          full.names = TRUE,
          include.dirs = FALSE,
          all.files = FALSE
        )
      }
    ),
    use.names = FALSE
  )

  files <- unique(
    c(
      root_files,
      stage_files
    )
  )

  sum(
    tolower(
      basename(files)
    ) %in%
      tolower(
        root_anchor_files
      )
  )
}

script_path <- get_current_script_path()
script_dir <- if (
  nzchar(script_path)
) {
  dirname(script_path)
} else {
  ""
}

desktop_candidate <- file.path(
  Sys.getenv("USERPROFILE"),
  "Desktop",
  "LZZ",
  "1、yssj"
)

if (dir.exists(USER_PROJECT_ROOT)) {
  PROJECT_ROOT <- normalizePath(
    USER_PROJECT_ROOT,
    winslash = "/",
    mustWork = TRUE
  )

  message(
    "Using fixed PROJECT_ROOT: ",
    PROJECT_ROOT
  )
} else {
  # Deliberately limited candidates: no multi-level ancestor traversal.
  fallback_candidates <- unique(
    c(
      desktop_candidate,
      script_dir,
      getwd(),
      if (nzchar(script_dir)) dirname(script_dir) else "",
      dirname(getwd())
    )
  )

  fallback_candidates <- fallback_candidates[
    !is.na(fallback_candidates) &
      nzchar(fallback_candidates) &
      dir.exists(fallback_candidates)
  ]

  message(
    "Fixed PROJECT_ROOT was not found. Checking ",
    length(fallback_candidates),
    " limited fallback location(s)..."
  )

  fallback_scores <- vapply(
    fallback_candidates,
    function(path) {
      message(
        "  Checking: ",
        normalizePath(
          path,
          winslash = "/",
          mustWork = FALSE
        )
      )

      fast_anchor_score(path)
    },
    numeric(1)
  )

  if (
    length(fallback_scores) &&
      max(fallback_scores) >= 6L
  ) {
    PROJECT_ROOT <- normalizePath(
      fallback_candidates[
        which.max(fallback_scores)
      ],
      winslash = "/",
      mustWork = TRUE
    )

    message(
      "Validated fallback PROJECT_ROOT: ",
      PROJECT_ROOT
    )
  } else {
    score_text <- if (
      length(fallback_scores)
    ) {
      paste0(
        "  - ",
        normalizePath(
          fallback_candidates,
          winslash = "/",
          mustWork = FALSE
        ),
        " | anchor files found: ",
        fallback_scores,
        collapse = "\n"
      )
    } else {
      "  No existing limited fallback directory was found."
    }

    stop(
      "Unable to locate the extracted Stage 1-12 project folder.\n",
      "Expected path: ",
      USER_PROJECT_ROOT,
      "\nCandidate assessment:\n",
      score_text,
      "\nUpdate USER_PROJECT_ROOT to the exact project directory.\n",
      "The script no longer scans C:/Users or the disk root automatically."
    )
  }
}

OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  OUTPUT_FOLDER_NAME
)

message("Validated PROJECT_ROOT: ", PROJECT_ROOT)

# ------------------------------- PACKAGES ------------------------------------
required_packages <- c(
  "ggplot2", "dplyr", "tidyr", "readr", "stringr", "forcats",
  "patchwork", "ggrepel", "scales", "purrr", "tibble", "ragg"
)

install_and_load <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    if (!AUTO_INSTALL_PACKAGES) {
      stop("Missing packages: ", paste(missing, collapse = ", "))
    }
    install.packages(missing, repos = "https://cloud.r-project.org", dependencies = TRUE)
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}
install_and_load(required_packages)

# ------------------------------ OUTPUT FOLDERS -------------------------------
# Always rebuild from a clean directory so stale figures or source-data files
# from an earlier failed run cannot enter the final package.
if (dir.exists(OUTPUT_DIR)) {
  unlink(
    OUTPUT_DIR,
    recursive = TRUE,
    force = TRUE
  )
}

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
DIR_MAIN   <- file.path(OUTPUT_DIR, "main_tif")
DIR_SUPP   <- file.path(OUTPUT_DIR, "supplementary_tif")
DIR_PANELS <- file.path(OUTPUT_DIR, "individual_panel_tif")
DIR_SOURCE <- file.path(OUTPUT_DIR, "source_data")
DIR_TABLES <- file.path(OUTPUT_DIR, "tables")
DIR_LOGS   <- file.path(OUTPUT_DIR, "logs")
for (d in c(DIR_MAIN, DIR_SUPP, DIR_PANELS, DIR_SOURCE, DIR_TABLES, DIR_LOGS)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------- INPUT FINDER --------------------------------
# Scan only Stage/output folders first. Raw GEO matrices and unrelated folders
# are not needed for publication plotting and can contain very large file trees.
message(
  "Indexing extracted Stage 1-12 result folders..."
)

first_level_dirs <- list.dirs(
  PROJECT_ROOT,
  recursive = FALSE,
  full.names = TRUE
)

stage_scan_dirs <- first_level_dirs[
  grepl(
    stage_directory_pattern,
    basename(first_level_dirs),
    perl = TRUE
  )
]

root_level_files <- list.files(
  PROJECT_ROOT,
  recursive = FALSE,
  full.names = TRUE,
  include.dirs = FALSE,
  all.files = FALSE
)

stage_level_files <- unlist(
  lapply(
    stage_scan_dirs,
    function(stage_dir) {
      message(
        "  Indexing: ",
        basename(stage_dir)
      )

      list.files(
        stage_dir,
        recursive = TRUE,
        full.names = TRUE,
        include.dirs = FALSE,
        all.files = FALSE
      )
    }
  ),
  use.names = FALSE
)

all_project_files <- unique(
  c(
    root_level_files,
    stage_level_files
  )
)

# Safety fallback: only when the targeted scan clearly found too few report
# files. This scans the project root itself, never its parent directories.
n_target_csv <- sum(
  grepl(
    "\\.csv$",
    all_project_files,
    ignore.case = TRUE
  )
)

if (
  length(stage_scan_dirs) == 0L ||
    n_target_csv < 20L
) {
  message(
    "Targeted Stage-folder scan found too few CSV files; ",
    "performing one recursive scan inside PROJECT_ROOT only..."
  )

  all_project_files <- list.files(
    PROJECT_ROOT,
    recursive = TRUE,
    full.names = TRUE,
    include.dirs = FALSE,
    all.files = FALSE
  )
}

all_project_files <- normalizePath(
  all_project_files,
  winslash = "/",
  mustWork = FALSE
)

output_norm <- normalizePath(
  OUTPUT_DIR,
  winslash = "/",
  mustWork = FALSE
)

all_project_files <- all_project_files[
  !startsWith(
    all_project_files,
    output_norm
  )
]

message(
  "Input index completed: ",
  length(all_project_files),
  " file(s) indexed."
)

input_manifest <- tibble::tibble(
  logical_name = character(),
  stage_tag = character(),
  file_name = character(),
  candidate_count = integer(),
  resolved_path = character(),
  status = character(),
  selection_note = character()
)

stage_alias_map <- list(
  "第一阶段" = c("第一阶段", "01_STAGE1_"),
  "第二阶段" = c("第二阶段", "02_STAGE2_"),
  "第三阶段" = c("第三阶段", "03_STAGE3_"),
  "第四阶段" = c("第四阶段", "04_STAGE4_"),
  "第五阶段" = c("第五阶段", "05_STAGE5_"),
  "第六阶段A" = c("第六阶段A", "06A_STAGE6A_"),
  "第六阶段B" = c("第六阶段B", "06B_STAGE6B_"),
  "第七阶段" = c("第七阶段", "07_STAGE7_"),
  "第八阶段" = c("第八阶段", "08_STAGE8_"),
  "第九阶段" = c("第九阶段", "09_STAGE9_"),
  "第十阶段" = c("第十阶段", "10_STAGE10_"),
  "第十一阶段" = c("第十一阶段", "11_STAGE11_"),
  "第十二阶段" = c("第十二阶段", "12_STAGE12_")
)

path_preference_score <- function(path, stage_tag = NULL) {
  score <- 0

  if (!is.null(stage_tag)) {
    aliases <- unique(c(
      stage_tag,
      stage_alias_map[[stage_tag]]
    ))

    if (any(vapply(
      aliases,
      function(alias) grepl(alias, path, fixed = TRUE),
      logical(1)
    ))) {
      score <- score + 100
    }
  }

  if (grepl("/report/", path, fixed = TRUE)) {
    score <- score + 20
  }

  if (grepl("STAGE12_GENE_PROGRAM_INTERPRETATION", path, fixed = TRUE)) {
    score <- score + 5
  }

  if (grepl(
    "backup|old|deprecated|archive|临时|旧版",
    path,
    ignore.case = TRUE
  )) {
    score <- score - 20
  }

  if (grepl("\\(1\\)|\\(2\\)|copy", path, ignore.case = TRUE)) {
    score <- score - 2
  }

  score
}

find_input <- function(
  file_name,
  stage_tag = NULL,
  required = TRUE,
  logical_name = file_name
) {
  hits <- all_project_files[
    tolower(basename(all_project_files)) ==
      tolower(file_name)
  ]

  if (!is.null(stage_tag) && length(hits)) {
    aliases <- unique(c(
      stage_tag,
      stage_alias_map[[stage_tag]]
    ))

    tagged <- hits[
      vapply(
        hits,
        function(path) {
          any(vapply(
            aliases,
            function(alias) grepl(alias, path, fixed = TRUE),
            logical(1)
          ))
        },
        logical(1)
      )
    ]

    if (length(tagged)) {
      hits <- tagged
    }
  }

  candidate_count <- length(hits)
  selection_note <- ""

  if (candidate_count > 1L) {
    pref_score <- vapply(
      hits,
      path_preference_score,
      numeric(1),
      stage_tag = stage_tag
    )

    mtime <- file.info(hits)$mtime
    ord <- order(
      pref_score,
      mtime,
      decreasing = TRUE,
      na.last = TRUE
    )

    hits <- hits[ord]
    selection_note <- paste0(
      "Multiple candidates found; selected highest-ranked/latest file. Other candidates: ",
      paste(hits[-1], collapse = " | ")
    )
  }

  path <- if (length(hits)) hits[[1]] else NA_character_
  status <- ifelse(is.na(path), "missing", "found")

  input_manifest <<- dplyr::bind_rows(
    input_manifest,
    tibble::tibble(
      logical_name = logical_name,
      stage_tag = ifelse(is.null(stage_tag), "", stage_tag),
      file_name = file_name,
      candidate_count = candidate_count,
      resolved_path = ifelse(is.na(path), "", path),
      status = status,
      selection_note = selection_note
    )
  )

  if (required && is.na(path)) {
    stop(
      "Required input not found: ", file_name,
      if (!is.null(stage_tag)) paste0(" [stage: ", stage_tag, "]") else "",
      "\nValidated PROJECT_ROOT: ", PROJECT_ROOT,
      "\nExtract the relevant stage check package or verify the stage output folder."
    )
  }

  path
}

read_csv_input <- function(
  file_name,
  stage_tag = NULL,
  required = TRUE,
  logical_name = file_name
) {
  path <- find_input(
    file_name,
    stage_tag,
    required,
    logical_name
  )

  if (is.na(path)) {
    return(NULL)
  }

  readr::read_csv(
    path,
    show_col_types = FALSE,
    progress = FALSE
  )
}

assert_columns <- function(data, required_columns, object_name) {
  if (is.null(data)) {
    stop(object_name, " is NULL.")
  }

  missing_columns <- setdiff(
    required_columns,
    names(data)
  )

  if (length(missing_columns)) {
    stop(
      object_name,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }

  invisible(TRUE)
}

# ------------------------------- CORE INPUTS ---------------------------------
# Stage 2-4: reference design and donor-held-out model
master_dictionary <- read_csv_input("09_master_sample_dictionary.csv", "第二阶段", FALSE, "Master sample dictionary")
ref_qc            <- read_csv_input("01_sample_QC_summary.csv", "第三阶段", TRUE, "GSE241132 sample QC")
ref_cell_counts   <- read_csv_input("03_main_celltype_counts_after_QC.csv", "第三阶段", TRUE, "GSE241132 cell counts")
lodo_pred         <- read_csv_input("01_LOODO_predictions.csv", "第四阶段", TRUE, "LODO predictions")
lodo_donor        <- read_csv_input("02_LOODO_donor_metrics.csv", "第四阶段", TRUE, "LODO donor metrics")
lodo_overall      <- read_csv_input("03_LOODO_overall_metrics.csv", "第四阶段", TRUE, "LODO overall metrics")
fixed_signatures  <- read_csv_input("05_fixed_final_stage_signatures.csv", "第四阶段", TRUE, "Frozen stage signatures")
sig_stability     <- read_csv_input("06_signature_stability_across_folds.csv", "第四阶段", TRUE, "Signature stability")
confusion         <- read_csv_input("07_LOODO_confusion_matrices.csv", "第四阶段", TRUE, "LODO confusion matrices")
stable_counts     <- read_csv_input("08_stable_signature_counts.csv", "第四阶段", TRUE, "Stable signature counts")

# Stage 5: discovery cohort
s5_qc             <- read_csv_input("03_sample_QC_summary.csv", "第五阶段", FALSE, "GSE163973 sample QC")
s5_cell_counts    <- read_csv_input("05_harmonized_celltype_counts.csv", "第五阶段", FALSE, "GSE163973 harmonized cell counts")
s5_scores         <- read_csv_input("09_GSE163973_patient_wound_state_scores.csv", "第五阶段", TRUE, "GSE163973 patient scores")
ref_calibration   <- read_csv_input("10_reference_sample_calibration_scores.csv", "第五阶段", TRUE, "Reference calibration scores")
s5_effects        <- read_csv_input("11_patient_level_group_comparisons.csv", "第五阶段", TRUE, "GSE163973 patient effects")

# Stage 6: independent replication
s6a_qc            <- read_csv_input("03_cell_calling_QC_fibroblast_summary.csv", "第六阶段A", TRUE, "GSE181316 cell calling QC")
s6_gene_coverage <- read_csv_input("02_signature_gene_coverage.csv", "第六阶段B", FALSE, "GSE181316 signature coverage")
s6_scores         <- read_csv_input("04_patient_level_scores_after_K3_merge.csv", "第六阶段B", TRUE, "GSE181316 patient scores")
s6_effects        <- read_csv_input("05_patient_level_group_comparisons.csv", "第六阶段B", TRUE, "GSE181316 patient effects")
k3_lr             <- read_csv_input("06_K3_left_right_concordance.csv", "第六阶段B", TRUE, "K3 left-right concordance")
replication_table <- read_csv_input("07_discovery_validation_replication_table.csv", "第六阶段B", TRUE, "Discovery-validation replication")

# Stage 7: third cohort support
s7_qc             <- read_csv_input("04_sample_QC_fibroblast_summary.csv", "第七阶段", TRUE, "GSE220300 sample QC")
s7_scores         <- read_csv_input("08_sample_level_wound_state_scores.csv", "第七阶段", TRUE, "GSE220300 sample scores")
s7_patient        <- read_csv_input("09_patient_lesion_and_scar_scores.csv", "第七阶段", TRUE, "GSE220300 patient lesion/scar scores")
s7_center_periph  <- read_csv_input("10_center_periphery_paired_differences.csv", "第七阶段", TRUE, "GSE220300 center-periphery pairs")
s7_inactive_scar  <- read_csv_input("11_inactive_lesion_vs_paired_mature_scar.csv", "第七阶段", TRUE, "GSE220300 inactive-scar pairs")
s7_activity       <- read_csv_input("12_active_vs_inactive_descriptive_comparisons.csv", "第七阶段", TRUE, "GSE220300 activity comparisons")
s7_lesion_scar    <- read_csv_input("13_keloid_lesion_vs_scar_reference_descriptive.csv", "第七阶段", TRUE, "GSE220300 lesion-scar effects")
s7_cross          <- read_csv_input("14_cross_cohort_direction_summary.csv", "第七阶段", TRUE, "Cross-cohort direction summary")

# Stage 8: normal wound spatial validation
s8_spot_preview   <- read_csv_input("07_spot_score_preview_max2000_per_section.csv", "第八阶段", TRUE, "Spatial spot source data")
s8_sections       <- read_csv_input("09_spatial_section_wound_state_scores.csv", "第八阶段", TRUE, "Spatial section scores")
s8_overall        <- read_csv_input("10_spatial_validation_overall_metrics.csv", "第八阶段", TRUE, "Spatial overall metrics")
s8_donor          <- read_csv_input("11_spatial_validation_donor_metrics.csv", "第八阶段", TRUE, "Spatial donor metrics")
s8_stage          <- read_csv_input("12_stage_level_spatial_summary.csv", "第八阶段", TRUE, "Spatial stage summary")

# Stage 9: disease contrast
s9_qc             <- read_csv_input("05_sample_QC_summary.csv", "第九阶段", TRUE, "GSE265972 sample QC")
s9_scores         <- read_csv_input("10_GSE265972_wound_state_scores.csv", "第九阶段", TRUE, "Venous ulcer patient scores")
s9_effects        <- read_csv_input("11_GSE265972_group_comparisons.csv", "第九阶段", TRUE, "Venous ulcer effects")
s9_cross          <- read_csv_input("12_cross_disease_failure_mode_summary.csv", "第九阶段", TRUE, "Cross-disease failure modes")

# Stage 10: exploratory mixed cross-modal evidence
s10_scores        <- read_csv_input("10_GSE181297_cross_modal_wound_state_scores.csv", "第十阶段", TRUE, "GSE181297 cross-modal scores")
s10_direction     <- read_csv_input("11_cross_modal_direction_summary.csv", "第十阶段", TRUE, "GSE181297 direction summary")
s10_spot_summary  <- read_csv_input("08_Visium_spot_distribution_summary.csv", "第十阶段", TRUE, "GSE181297 spot summary")
s10_qc            <- read_csv_input("06_Visium_QC_fibroblast_spot_summary.csv", "第十阶段", TRUE, "GSE181297 Visium QC")

# Stage 11: final evidence lock
core_patient      <- read_csv_input("02_core_patient_level_data.csv", "第十一阶段", TRUE, "Core patient-level data")
core_effects      <- read_csv_input("03_core_cohort_specific_effects.csv", "第十一阶段", TRUE, "Core cohort effects")
core_integrated   <- read_csv_input("04_stratified_exact_core_integration.csv", "第十一阶段", TRUE, "Stratified exact integration")
perm_distrib      <- read_csv_input("05_exact_joint_permutation_distributions.csv", "第十一阶段", TRUE, "Exact permutation distributions")
evidence_registry<- read_csv_input("06_supportive_evidence_registry.csv", "第十一阶段", FALSE, "Supportive evidence registry")
failure_plane     <- read_csv_input("10_failure_mode_contrast_source_data.csv", "第十一阶段", TRUE, "Failure-mode contrast")

# Stage 12: gene programs and pathways
all_gene_effects  <- read_csv_input("05_all_gene_core_effects_and_exact_integration.csv", "第十二阶段", TRUE, "All-gene core effects")
third_gene_effects<- read_csv_input("06_GSE220300_supportive_gene_effects.csv", "第十二阶段", TRUE, "Third-cohort gene effects")
gene_contrib      <- read_csv_input("07_fixed_signature_gene_contribution_table.csv", "第十二阶段", TRUE, "Fixed-signature gene contributions")
late_genes        <- read_csv_input("08_consensus_persistent_late_remodeling_genes.csv", "第十二阶段", TRUE, "Persistent late-remodeling genes")
skin_loss_genes   <- read_csv_input("09_consensus_loss_of_skin_homeostasis_genes.csv", "第十二阶段", TRUE, "Skin-homeostasis loss genes")
stage_competition <- read_csv_input("10_stage_signature_competition_summary.csv", "第十二阶段", TRUE, "Stage signature competition")
state_panels      <- read_csv_input("11_predefined_fibroblast_state_panels.csv", "第十二阶段", FALSE, "Predefined fibroblast panels")
panel_effects     <- read_csv_input("13_fibroblast_state_panel_patient_effects.csv", "第十二阶段", TRUE, "Fibroblast state panel effects")
panel_scores      <- read_csv_input("14_fibroblast_state_panel_sample_scores.csv", "第十二阶段", TRUE, "Fibroblast panel sample scores")
stage_gsea        <- read_csv_input("15_frozen_stage_signature_GSEA.csv", "第十二阶段", TRUE, "Frozen stage signature GSEA")
hallmark_gsea     <- read_csv_input("17_Hallmark_GSEA.csv", "第十二阶段", TRUE, "Hallmark GSEA")
gobp_gsea         <- read_csv_input("18_GO_BP_GSEA.csv", "第十二阶段", TRUE, "GO BP GSEA")
reactome_gsea     <- read_csv_input("19_Reactome_GSEA.csv", "第十二阶段", TRUE, "Reactome GSEA")
late_hallmark_ora <- read_csv_input("20_persistent_late_Hallmark_ORA.csv", "第十二阶段", TRUE, "Persistent late Hallmark ORA")
late_gobp_ora     <- read_csv_input("21_persistent_late_GO_BP_ORA.csv", "第十二阶段", TRUE, "Persistent late GO BP ORA")
late_reactome_ora <- read_csv_input("22_persistent_late_Reactome_ORA.csv", "第十二阶段", TRUE, "Persistent late Reactome ORA")
skin_hallmark_ora <- read_csv_input("23_skin_homeostasis_loss_Hallmark_ORA.csv", "第十二阶段", TRUE, "Skin loss Hallmark ORA")
skin_gobp_ora     <- read_csv_input("24_skin_homeostasis_loss_GO_BP_ORA.csv", "第十二阶段", TRUE, "Skin loss GO BP ORA")
skin_reactome_ora <- read_csv_input("25_skin_homeostasis_loss_Reactome_ORA.csv", "第十二阶段", TRUE, "Skin loss Reactome ORA")

# Fail early with precise, human-readable structure checks.
assert_columns(ref_qc, c("sample_id", "donor", "condition", "qc_retention_rate"), "ref_qc")
assert_columns(ref_cell_counts, c("sample_id", "donor", "condition", "main_cell_type", "n_cells"), "ref_cell_counts")
assert_columns(lodo_pred, c("cell_type", "held_out_donor", "true_stage", "predicted_stage", "true_position", "predicted_position"), "lodo_pred")
assert_columns(lodo_donor, c("cell_type", "held_out_donor", "exact_accuracy", "mean_ordinal_error", "spearman_true_vs_predicted_position"), "lodo_donor")
assert_columns(confusion, c("cell_type", "true_stage", "predicted_stage", "Freq"), "confusion")
assert_columns(stable_counts, c("cell_type", "stage", "n_stable_genes"), "stable_counts")
assert_columns(ref_calibration, c("cell_type", "sample_id", "donor", "true_stage", "z_Skin", "z_Wound1", "z_Wound7", "z_Wound30"), "ref_calibration")
assert_columns(s5_scores, c("cell_type", "sample_id", "group", "remodeling_completion", "ordinal_position", "z_Skin", "z_Wound30"), "s5_scores")
assert_columns(s6_scores, c("analysis_set", "patient_id", "group", "late_remodeling_state", "ordinal_wound_state_position", "z_Skin", "z_Wound30"), "s6_scores")
assert_columns(s7_patient, c("analysis_set", "sample_key", "patient_id", "aggregate_state", "late_remodeling_state", "z_Skin", "z_Wound30"), "s7_patient")
assert_columns(s8_spot_preview, c("donor", "condition", "pxl_row_in_fullres", "pxl_col_in_fullres", "high_fibroblast_enriched", "spot_late_remodeling"), "s8_spot_preview")
assert_columns(s8_sections, c("analysis_set", "donor", "condition", "stage_code", "nearest_stage", "ordinal_wound_state_position"), "s8_sections")
assert_columns(s9_scores, c("analysis_set", "sample_key", "group", "early_state_persistence", "late_remodeling_state", "late_vs_early_balance"), "s9_scores")
assert_columns(core_patient, c("cohort", "patient_id", "group", "late_remodeling_state", "early_state_persistence"), "core_patient")
assert_columns(core_integrated, c("metric", "equal_weight_mean_difference", "exact_stratified_two_sided_p", "endpoint_role"), "core_integrated")
assert_columns(all_gene_effects, c("gene", "mean_difference_GSE163973", "mean_difference_GSE181316", "integrated_rank_score"), "all_gene_effects")
assert_columns(late_genes, c("gene", "mean_difference_GSE163973", "mean_difference_GSE181316", "minimum_oriented_core_effect"), "late_genes")
assert_columns(skin_loss_genes, c("gene", "mean_difference_GSE163973", "mean_difference_GSE181316", "minimum_oriented_core_effect"), "skin_loss_genes")

if (nrow(late_genes) != 51L || nrow(skin_loss_genes) != 29L) {
  stop(
    "Stage 12 is not the locked V1.2 result. Expected 51 persistent-late genes and ",
    "29 Skin-homeostasis-loss genes; found ",
    nrow(late_genes), " and ", nrow(skin_loss_genes), "."
  )
}

readr::write_csv(
  input_manifest,
  file.path(DIR_LOGS, "resolved_input_manifest.csv")
)

# ------------------------------- CONSTANTS -----------------------------------
stage_levels <- c("Skin", "Wound1", "Wound7", "Wound30")
stage_labels <- c(Skin = "Skin", Wound1 = "Day 1", Wound7 = "Day 7", Wound30 = "Day 30")
stage_colors <- c(
  Skin = "#4E5D66", Wound1 = "#D95F02", Wound7 = "#E6AB02", Wound30 = "#1F78B4"
)
group_colors <- c(
  normal_scar = "#5A6268", healthy_skin = "#5A6268", normal_skin = "#5A6268",
  adjacent_normal = "#7F878D", keloid = "#B2182B", venous_ulcer = "#7B3294",
  scar_reference = "#5A6268", keloid_lesion = "#B2182B"
)
cohort_colors <- c(GSE163973 = "#1B9E77", GSE181316 = "#7570B3", GSE220300 = "#D95F02")

metric_labels <- c(
  late_remodeling_state = "Late-remodeling state",
  remodeling_completion = "Late-remodeling state",
  ordinal_wound_state_position = "Ordinal wound-state position",
  ordinal_position = "Ordinal wound-state position",
  z_Wound30 = "Day-30-like program",
  wound_activation = "Overall wound activation",
  z_Skin = "Skin-like homeostasis",
  skin_return_state = "Skin-return state",
  early_state_persistence = "Early-state persistence",
  off_trajectory_ratio = "Off-trajectory ratio",
  late_vs_early_balance = "Late-minus-early balance"
)

pretty_metric <- function(x) {
  out <- unname(metric_labels[x])
  out[is.na(out)] <- stringr::str_replace_all(x[is.na(out)], "_", " ")
  out
}

pretty_group <- function(x) {
  dplyr::recode(
    x,
    normal_scar = "Normal scar", healthy_skin = "Normal scar",
    normal_skin = "Normal skin", adjacent_normal = "Adjacent normal",
    keloid = "Keloid", venous_ulcer = "Venous ulcer", scar_reference = "Scar reference",
    keloid_lesion = "Keloid lesion",
    .default = stringr::str_to_sentence(stringr::str_replace_all(x, "_", " "))
  )
}

clean_pathway <- function(x) {
  x |>
    stringr::str_replace("^HALLMARK_", "") |>
    stringr::str_replace("^GOBP_", "") |>
    stringr::str_replace("^REACTOME_", "") |>
    stringr::str_replace_all("_", " ") |>
    stringr::str_to_sentence()
}


clean_pathway_short <- function(x, width = 34) {
  y <- clean_pathway(x)
  y <- stringr::str_replace_all(y, "\\bextracellular matrix\\b", "ECM")
  y <- stringr::str_replace_all(y, "\\bbiological process\\b", "")
  y <- stringr::str_replace_all(y, "\\borganization\\b", "organisation")
  y <- stringr::str_replace_all(y, "\\bmodifying enzymes\\b", "modifying enzymes")
  y <- stringr::str_replace_all(y, "\\band\\b", "&")
  y <- stringr::str_replace_all(y, "\\s+", " ")
  y <- stringr::str_squish(y)
  stringr::str_wrap(y, width = width)
}

safe_file_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "_+", "_")
  x <- stringr::str_replace_all(x, "^_|_$", "")
  ifelse(nzchar(x), x, "NA")
}

confusion_matrix_clean <- function(
  data,
  x_col,
  y_col,
  count_col,
  fill_col,
  title = "Confusion matrix",
  x_lab = "Predicted stage",
  y_lab = "True stage",
  x_labels = stage_labels,
  y_labels = stage_labels
) {
  d <- data |>
    dplyr::mutate(
      .x = .data[[x_col]],
      .y = .data[[y_col]],
      .n = as.numeric(.data[[count_col]]),
      .fill = as.numeric(.data[[fill_col]]),
      .label = ifelse(
        !is.na(.n) & .n > 0,
        paste0(
          .n,
          "\n",
          scales::percent(.fill, accuracy = 1)
        ),
        ""
      )
    )

  ggplot2::ggplot(
    d,
    ggplot2::aes(.x, .y, fill = .fill)
  ) +
    ggplot2::geom_tile(
      colour = "white",
      linewidth = 0.55
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = .label),
      size = 3.1,
      lineheight = 0.88,
      family = BASE_FONT,
      colour = "black"
    ) +
    ggplot2::scale_fill_gradient(
      low = "white",
      high = "#2166AC",
      limits = c(0, 1),
      labels = scales::percent,
      name = "Row %",
      guide = ggplot2::guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barheight = grid::unit(42, "pt"),
        barwidth = grid::unit(7, "pt"),
        ticks.linewidth = 0.35
      )
    ) +
    ggplot2::scale_x_discrete(labels = x_labels) +
    ggplot2::scale_y_discrete(labels = y_labels) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      title = title,
      x = x_lab,
      y = y_lab
    ) +
    theme_pub(base_size = 9.2, title_size = 9.4) +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x = ggplot2::element_text(
        angle = 0,
        hjust = 0.5,
        vjust = 0.5
      ),
      axis.text.y = ggplot2::element_text(
        angle = 0,
        hjust = 1
      ),
      plot.margin = ggplot2::margin(8, 12, 8, 10),
      legend.margin = ggplot2::margin(0, 0, 0, 4)
    )
}

pathway_resource_plot <- function(
  data,
  resource_name,
  title,
  n_terms = 12,
  wrap_width = 34
) {
  d <- data |>
    dplyr::filter(
      resource == resource_name,
      !is.na(NES),
      !is.na(p.adjust),
      p.adjust < 0.05
    ) |>
    dplyr::arrange(
      p.adjust,
      dplyr::desc(abs(NES))
    ) |>
    dplyr::mutate(
      pathway_clean = clean_pathway_short(Description, width = wrap_width)
    ) |>
    dplyr::distinct(
      pathway_clean,
      .keep_all = TRUE
    ) |>
    dplyr::slice_head(n = n_terms) |>
    dplyr::mutate(
      pathway_clean = forcats::fct_reorder(pathway_clean, NES),
      logFDR = -log10(
        pmax(p.adjust, .Machine$double.xmin)
      )
    )

  if (!nrow(d)) {
    return(
      panel_placeholder(
        title,
        "No FDR-significant term"
      )
    )
  }

  ggplot2::ggplot(
    d,
    ggplot2::aes(
      NES,
      pathway_clean,
      size = logFDR,
      colour = NES
    )
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.35,
      colour = "#8A8A8A"
    ) +
    ggplot2::geom_point(alpha = 0.95) +
    ggplot2::scale_colour_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      name = "NES"
    ) +
    ggplot2::scale_size_continuous(
      range = c(2.8, 6.0),
      name = "-log10 FDR"
    ) +
    ggplot2::labs(
      title = title,
      x = "Normalized enrichment score",
      y = NULL
    ) +
    theme_pub(base_size = 8.4, title_size = 9.0) +
    ggplot2::theme(
      legend.position = "right",
      axis.text.y = ggplot2::element_text(size = 7.8, lineheight = 0.88),
      plot.margin = ggplot2::margin(8, 14, 8, 12)
    )
}

# ------------------------------- THEMES --------------------------------------
wrap_title <- function(x, width = 28) {
  if (is.null(x) || is.na(x) || !nzchar(x)) return(x)
  stringr::str_wrap(x, width = width)
}

theme_pub <- function(base_size = 8.6, title_size = 8.8) {
  ggplot2::theme_classic(base_size = base_size, base_family = BASE_FONT) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.42, colour = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.36, colour = "black"),
      axis.ticks.length = grid::unit(1.8, "pt"),
      axis.title = ggplot2::element_text(
        size = base_size + 0.15, colour = "black",
        margin = ggplot2::margin(3, 3, 3, 3)
      ),
      axis.text = ggplot2::element_text(size = base_size - 0.35, colour = "black"),
      plot.title = ggplot2::element_text(
        size = title_size, face = "bold", hjust = 0,
        margin = ggplot2::margin(0, 0, 3, 0), lineheight = 0.95
      ),
      plot.subtitle = ggplot2::element_text(size = base_size - 0.55, colour = "#4D4D4D", hjust = 0),
      plot.caption = ggplot2::element_text(size = base_size - 1.0, colour = "#606060", hjust = 0),
      plot.margin = ggplot2::margin(8, 10, 8, 10),
      legend.position = "top",
      legend.box = "horizontal",
      legend.justification = "left",
      legend.title = ggplot2::element_text(size = base_size - 0.35, face = "bold"),
      legend.text = ggplot2::element_text(size = base_size - 0.75),
      legend.key.height = grid::unit(8, "pt"),
      legend.key.width = grid::unit(10, "pt"),
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.box.margin = ggplot2::margin(0, 0, 0, 0),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        size = base_size - 0.55, face = "bold", colour = "black",
        margin = ggplot2::margin(1, 1, 1, 1)
      )
    )
}

theme_heatmap <- function(base_size = 8.0, title_size = 8.6) {
  ggplot2::theme_minimal(base_size = base_size, base_family = BASE_FONT) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_text(colour = "black", size = base_size - 0.2),
      axis.text.x = ggplot2::element_text(angle = 40, hjust = 1, vjust = 1),
      plot.title = ggplot2::element_text(size = title_size, face = "bold", hjust = 0, margin = ggplot2::margin(0,0,3,0)),
      legend.title = ggplot2::element_text(face = "bold", size = base_size - 0.2),
      legend.text = ggplot2::element_text(size = base_size - 0.6),
      legend.key.height = grid::unit(8, "pt"),
      legend.key.width = grid::unit(10, "pt"),
      plot.margin = ggplot2::margin(8, 10, 8, 10)
    )
}

theme_patch_tags <- ggplot2::theme(
  plot.tag = ggplot2::element_text(
    family = BASE_FONT, face = "bold", size = 11.2, colour = "black",
    margin = ggplot2::margin(0, 4, 2, 0)
  ),
  plot.tag.position = c(0, 1)
)

panel_placeholder <- function(title, text = "Input unavailable") {
  ggplot2::ggplot() +
    ggplot2::annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
                      fill = "#F7F7F7", colour = "#B0B0B0") +
    ggplot2::annotate("text", x = 0.5, y = 0.55, label = title,
                      family = BASE_FONT, fontface = "bold", size = 3.2) +
    ggplot2::annotate("text", x = 0.5, y = 0.40, label = text,
                      family = BASE_FONT, size = 2.8, colour = "#666666") +
    ggplot2::coord_cartesian(xlim = c(0,1), ylim = c(0,1), expand = FALSE) +
    ggplot2::theme_void()
}

# ---------------------------- SOURCE DATA EXPORT -----------------------------
source_index <- tibble::tibble(
  figure = character(), panel = character(), file = character(),
  description = character(), n_rows = integer(), n_columns = integer()
)

write_source <- function(data, figure, panel, description) {
  stopifnot(is.data.frame(data))
  safe_fig <- stringr::str_replace_all(figure, "[^A-Za-z0-9_]+", "_")
  safe_pan <- stringr::str_replace_all(panel, "[^A-Za-z0-9_]+", "_")
  out <- file.path(DIR_SOURCE, paste0(safe_fig, "_", safe_pan, ".csv"))
  readr::write_csv(data, out, na = "")
  source_index <<- dplyr::bind_rows(
    source_index,
    tibble::tibble(
      figure = figure, panel = panel, file = basename(out), description = description,
      n_rows = nrow(data), n_columns = ncol(data)
    )
  )
  invisible(out)
}

# ------------------------------- TIFF OUTPUT ---------------------------------
save_tiff <- function(
  plot,
  filename,
  width,
  height,
  dpi = DPI,
  dir = DIR_MAIN
) {
  path <- file.path(dir, filename)
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)

  device_open <- FALSE

  if (requireNamespace("ragg", quietly = TRUE)) {
    device_result <- try(
      {
        ragg::agg_tiff(
          filename = path,
          width = width,
          height = height,
          units = "in",
          res = dpi,
          background = "white",
          compression = "lzw"
        )
        device_open <- TRUE
      },
      silent = TRUE
    )

    if (inherits(device_result, "try-error")) {
      device_open <- FALSE
    }
  }

  if (!device_open) {
    cairo_result <- try(
      {
        grDevices::tiff(
          filename = path,
          width = width,
          height = height,
          units = "in",
          res = dpi,
          compression = "lzw",
          type = "cairo",
          bg = "white"
        )
        device_open <- TRUE
      },
      silent = TRUE
    )

    if (inherits(cairo_result, "try-error")) {
      device_open <- FALSE
    }
  }

  if (!device_open) {
    grDevices::tiff(
      filename = path,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      compression = "lzw",
      bg = "white"
    )
  }

  print(plot)
  grDevices::dev.off()

  if (!file.exists(path) || file.info(path)$size <= 0) {
    stop("TIFF output failed: ", path)
  }

  invisible(path)
}

save_panel_set <- function(
  plots,
  figure_name,
  width = PANEL_BASE_WIDTH,
  height = PANEL_BASE_HEIGHT
) {
  if (!SAVE_PANEL_TIFFS) {
    return(invisible(NULL))
  }

  panel_names <- names(plots)

  for (i in seq_along(plots)) {
    panel <- panel_names[[i]]
    panel_width <- if (length(width) == 1L) width else width[[i]]
    panel_height <- if (length(height) == 1L) height else height[[i]]

    save_tiff(
      plots[[i]],
      paste0(figure_name, "_", panel, ".tif"),
      width = panel_width,
      height = panel_height,
      dir = DIR_PANELS
    )
  }

  invisible(NULL)
}

# ---------------------------- GENERAL PLOT HELPERS ---------------------------
patient_dot_plot <- function(data, metric, title, ylab = NULL, group_col = "group",
                             id_col = "patient_id", order = c("normal_scar", "keloid")) {
  d <- data |>
    dplyr::filter(!is.na(.data[[metric]])) |>
    dplyr::mutate(
      group_plot = factor(.data[[group_col]], levels = order),
      group_label = factor(
        pretty_group(as.character(group_plot)),
        levels = pretty_group(order)
      )
    )
  ggplot2::ggplot(d, ggplot2::aes(x = group_label, y = .data[[metric]], fill = group_plot)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, linetype = "dashed", colour = "#9A9A9A") +
    ggplot2::geom_point(
      shape = 21, size = 3.2, stroke = 0.55, colour = "black",
      position = ggplot2::position_jitter(width = 0.07, height = 0, seed = 20260618)
    ) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 95, size = 9, colour = "black") +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = .data[[id_col]]), size = 2.55, family = BASE_FONT,
      box.padding = 0.25, point.padding = 0.25, min.segment.length = 0,
      segment.size = 0.35, max.overlaps = Inf, show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = group_colors, drop = FALSE) +
    ggplot2::labs(title = title, x = NULL, y = ifelse(is.null(ylab), pretty_metric(metric), ylab)) +
    theme_pub() +
    ggplot2::theme(legend.position = "none")
}

heatmap_stage_scores <- function(data, sample_col, group_col = NULL, title = NULL,
                                 score_cols = c("z_Skin", "z_Wound1", "z_Wound7", "z_Wound30"),
                                 sample_order = NULL) {
  d <- data |>
    dplyr::select(dplyr::all_of(c(sample_col, group_col, score_cols))) |>
    tidyr::pivot_longer(dplyr::all_of(score_cols), names_to = "stage", values_to = "score") |>
    dplyr::mutate(
      stage = stringr::str_remove(stage, "^z_"),
      stage = factor(stage, levels = stage_levels),
      sample = .data[[sample_col]]
    )
  if (!is.null(sample_order)) d$sample <- factor(d$sample, levels = sample_order)
  else d$sample <- factor(d$sample, levels = unique(d$sample))
  ggplot2::ggplot(d, ggplot2::aes(x = sample, y = stage, fill = score)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.45) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", score)), size = 2.2, family = BASE_FONT) +
    ggplot2::scale_y_discrete(labels = stage_labels) +
    ggplot2::scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
      name = "Standardized\nscore"
    ) +
    ggplot2::labs(title = title) +
    theme_heatmap() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

stage_probability_bars <- function(data, sample_col, group_col = NULL, title = NULL,
                                   weight_prefix = "weight_") {
  cols <- paste0(weight_prefix, stage_levels)
  d <- data |>
    dplyr::select(dplyr::all_of(c(sample_col, group_col, cols))) |>
    tidyr::pivot_longer(dplyr::all_of(cols), names_to = "stage", values_to = "weight") |>
    dplyr::mutate(
      stage = stringr::str_remove(stage, paste0("^", weight_prefix)),
      stage = factor(stage, levels = stage_levels),
      sample = factor(.data[[sample_col]], levels = unique(.data[[sample_col]]))
    )
  ggplot2::ggplot(d, ggplot2::aes(x = sample, y = weight, fill = stage)) +
    ggplot2::geom_col(width = 0.78, colour = "white", linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = stage_colors, labels = stage_labels) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = c(0,0)) +
    ggplot2::labs(title = title, x = NULL, y = "Stage weight", fill = "Reference stage") +
    theme_pub() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "top")
}

forest_plot <- function(data, estimate_col, label_col, title, xlab,
                        cohort_col = NULL, point_size = 2.8, zero_line = TRUE,
                        colour_values = NULL) {
  d <- data |>
    dplyr::mutate(label_plot = forcats::fct_rev(factor(.data[[label_col]], levels = unique(.data[[label_col]]))))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data[[estimate_col]], y = label_plot))
  if (zero_line) p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "#8C8C8C")
  if (!is.null(cohort_col)) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(colour = .data[[cohort_col]]), size = point_size) +
      ggplot2::scale_colour_manual(values = colour_values %||% cohort_colors)
  } else {
    p <- p + ggplot2::geom_point(size = point_size, colour = "#2B2B2B")
  }
  p +
    ggplot2::labs(title = title, x = xlab, y = NULL, colour = NULL) +
    theme_pub() +
    ggplot2::theme(legend.position = "top")
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# -------------------------- DATA STANDARDIZATION -----------------------------
ref_fib <- ref_calibration |>
  dplyr::filter(cell_type == "Fibroblast") |>
  dplyr::mutate(
    true_stage = factor(true_stage, levels = stage_levels),
    donor = factor(donor),
    stage_code = as.numeric(true_stage) - 1
  )

lodo_fib <- lodo_pred |>
  dplyr::filter(cell_type == "Fibroblast") |>
  dplyr::mutate(
    true_stage = factor(true_stage, levels = stage_levels),
    predicted_stage = factor(predicted_stage, levels = stage_levels)
  )

s5_fib <- s5_scores |>
  dplyr::filter(cell_type == "Fibroblast") |>
  dplyr::mutate(
    late_remodeling_state = remodeling_completion,
    ordinal_wound_state_position = ordinal_position,
    patient_id = sample_id,
    group = dplyr::recode(group, healthy_skin = "normal_scar", .default = group),
    group = factor(group, levels = c("normal_scar", "keloid"))
  )

s6_primary <- s6_scores |>
  dplyr::filter(analysis_set == "high_specificity") |>
  dplyr::mutate(
    group = dplyr::recode(group, healthy_skin = "normal_scar", normal_skin = "normal_scar", .default = group),
    group = factor(group, levels = c("normal_scar", "keloid"))
  )

core_patient <- core_patient |>
  dplyr::mutate(group = factor(group, levels = c("normal_scar", "keloid")))

s8_primary <- s8_sections |>
  dplyr::filter(analysis_set == "high_specificity_spots") |>
  dplyr::mutate(
    condition = factor(condition, levels = stage_levels),
    nearest_stage = factor(nearest_stage, levels = stage_levels)
  )

s9_primary <- s9_scores |>
  dplyr::filter(analysis_set == "FB1_FB2_combined_primary") |>
  dplyr::mutate(group = factor(group, levels = c("normal_skin", "venous_ulcer")))

# =============================== FIGURE 1 ====================================
message("Building Figure 1...")

# 1A Reference fibroblast yield across donors and stages
ref_fib_counts <- ref_cell_counts |>
  dplyr::filter(main_cell_type == "Fibroblast") |>
  dplyr::mutate(
    condition = factor(condition, levels = stage_levels),
    donor = factor(donor)
  )

write_source(
  ref_fib_counts,
  "Figure1",
  "A",
  "Fibroblast cell counts across all reference donors and wound stages"
)

p1a <- ggplot2::ggplot(
  ref_fib_counts,
  ggplot2::aes(condition, n_cells, colour = donor, group = donor)
) +
  ggplot2::geom_line(linewidth = 0.65) +
  ggplot2::geom_point(size = 3.0, stroke = 0.55) +
  ggplot2::scale_x_discrete(labels = stage_labels) +
  ggplot2::scale_y_log10(labels = scales::comma) +
  ggplot2::labs(
    title = "Reference fibroblast yield",
    x = NULL,
    y = "Fibroblasts per donor-stage sample",
    colour = "Donor"
  ) +
  theme_pub() +
  ggplot2::theme(
    legend.position = "top",
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1)
  )

# 1B PCA of frozen stage scores
pca_ref <- stats::prcomp(
  ref_fib |> dplyr::select(z_Skin, z_Wound1, z_Wound7, z_Wound30),
  center = TRUE, scale. = TRUE
)
pca_ref_df <- dplyr::bind_cols(
  ref_fib |> dplyr::select(sample_id, donor, true_stage),
  as.data.frame(pca_ref$x[, 1:2, drop = FALSE])
)
var_exp <- round(100 * summary(pca_ref)$importance[2, 1:2], 1)
write_source(pca_ref_df, "Figure1", "B", "PCA coordinates from frozen four-stage fibroblast scores")
p1b <- ggplot2::ggplot(pca_ref_df, ggplot2::aes(PC1, PC2, colour = true_stage, shape = donor, group = donor)) +
  ggplot2::geom_path(linewidth = 0.55, colour = "#7A7A7A", alpha = 0.75, arrow = grid::arrow(length = grid::unit(0.08, "inches"))) +
  ggplot2::geom_point(size = 3.2, stroke = 0.75) +
  ggplot2::scale_colour_manual(values = stage_colors, labels = stage_labels) +
  ggplot2::labs(
    title = "Reference projection",
    x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)"),
    colour = "Stage", shape = "Donor"
  ) + theme_pub() + ggplot2::theme(legend.position = "top")

# 1C Reference heatmap
ref_order <- ref_fib |>
  dplyr::arrange(donor, true_stage) |>
  dplyr::pull(sample_id)
p1c <- heatmap_stage_scores(ref_fib, "sample_id", title = "Reference stage scores", sample_order = ref_order)
write_source(
  ref_fib |> dplyr::select(sample_id, donor, true_stage, z_Skin, z_Wound1, z_Wound7, z_Wound30),
  "Figure1", "C", "Frozen stage scores for all 12 reference pseudobulk samples"
)

# 1D Confusion matrix
conf_fib <- confusion |>
  dplyr::filter(cell_type == "Fibroblast") |>
  dplyr::mutate(
    true_stage = factor(true_stage, levels = stage_levels),
    predicted_stage = factor(predicted_stage, levels = stage_levels)
  ) |>
  dplyr::group_by(true_stage) |>
  dplyr::mutate(row_percent = Freq / sum(Freq)) |>
  dplyr::ungroup()
write_source(conf_fib, "Figure1", "D", "Donor-held-out confusion matrix")
p1d <- confusion_matrix_clean(
  data = conf_fib,
  x_col = "predicted_stage",
  y_col = "true_stage",
  count_col = "Freq",
  fill_col = "row_percent",
  title = "LODO confusion matrix",
  x_lab = "Predicted stage",
  y_lab = "True stage"
)

# 1E Donor metrics
lodo_donor_fib <- lodo_donor |>
  dplyr::filter(cell_type == "Fibroblast") |>
  dplyr::select(held_out_donor, exact_accuracy, mean_ordinal_error, spearman_true_vs_predicted_position) |>
  tidyr::pivot_longer(-held_out_donor, names_to = "metric", values_to = "value") |>
  dplyr::mutate(
    metric = factor(metric,
      levels = c("exact_accuracy", "mean_ordinal_error", "spearman_true_vs_predicted_position"),
      labels = c("Exact accuracy", "Mean ordinal error", "Spearman rho")
    )
  )
write_source(lodo_donor_fib, "Figure1", "E", "Per-donor held-out performance metrics")
p1e <- ggplot2::ggplot(lodo_donor_fib, ggplot2::aes(held_out_donor, value, group = 1)) +
  ggplot2::geom_hline(yintercept = c(0, 0.5, 1), linewidth = 0.30, colour = "#E0E0E0") +
  ggplot2::geom_line(linewidth = 0.55, colour = "#4D4D4D") +
  ggplot2::geom_point(size = 2.9, shape = 21, fill = "white", stroke = 0.65) +
  ggplot2::facet_wrap(~ metric, nrow = 1, scales = "free_y") +
  ggplot2::labs(title = "Unseen-donor performance", x = "Held-out donor", y = NULL) +
  theme_pub() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

# 1F Ordinal prediction
write_source(lodo_fib, "Figure1", "F", "True and predicted ordinal wound-state positions")
p1f <- ggplot2::ggplot(lodo_fib, ggplot2::aes(true_position, predicted_position, colour = true_stage, group = held_out_donor)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.55, colour = "#777777") +
  ggplot2::geom_line(linewidth = 0.48, colour = "#8C8C8C", alpha = 0.7) +
  ggplot2::geom_point(size = 3.0, stroke = 0.5) +
  ggplot2::scale_colour_manual(values = stage_colors, labels = stage_labels) +
  ggplot2::scale_x_continuous(breaks = 0:3, labels = stage_labels, limits = c(-0.1, 3.1)) +
  ggplot2::scale_y_continuous(breaks = 0:3, labels = stage_labels, limits = c(-0.1, 3.1)) +
  ggplot2::coord_equal() +
  ggplot2::labs(title = "Ordinal prediction", x = "True state", y = "Predicted position", colour = "Stage") +
  theme_pub() + ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

# 1G Stable genes
stable_fib <- stable_counts |>
  dplyr::filter(cell_type == "Fibroblast") |>
  dplyr::mutate(stage = factor(stage, levels = stage_levels))
write_source(stable_fib, "Figure1", "G", "Number of cross-fold stable genes in each frozen fibroblast signature")
p1g <- ggplot2::ggplot(stable_fib, ggplot2::aes(stage, n_stable_genes, fill = stage)) +
  ggplot2::geom_col(width = 0.68, colour = "black", linewidth = 0.45) +
  ggplot2::geom_text(ggplot2::aes(label = n_stable_genes), vjust = -0.35, size = 3, family = BASE_FONT) +
  ggplot2::scale_fill_manual(values = stage_colors, guide = "none") +
  ggplot2::scale_x_discrete(labels = stage_labels) +
  ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
  ggplot2::labs(title = "Stable gene counts", x = NULL, y = "Stable genes") +
  theme_pub() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

fig1_design <- "
AABBB
CCCCC
DDEEE
FFGGG
"
figure1 <- p1a + p1b + p1c + p1d + p1e + p1f + p1g +
  patchwork::plot_layout(design = fig1_design, guides = "keep") +
  patchwork::plot_annotation(tag_levels = "A", theme = theme_patch_tags)
if (SAVE_COMPOSITE_TIFFS) save_tiff(figure1, "Figure1.tif", width = 7.2, height = 9.4)
save_panel_set(
  list(A = p1a, B = p1b, C = p1c, D = p1d, E = p1e, F = p1f, G = p1g),
  "Figure1",
  width = c(6.2, 6.2, 6.4, 5.9, 5.8, 6.0, 5.4),
  height = c(5.0, 4.8, 5.0, 5.0, 4.4, 4.6, 4.2)
)

# =============================== FIGURE 2 ====================================
message("Building Figure 2...")

s5_order <- s5_fib |>
  dplyr::arrange(group, dplyr::desc(late_remodeling_state)) |>
  dplyr::pull(patient_id)

# 2A stage-score heatmap
p2a <- heatmap_stage_scores(s5_fib, "patient_id", "group", "Discovery stage scores", sample_order = s5_order)
write_source(s5_fib |> dplyr::select(patient_id, group, z_Skin, z_Wound1, z_Wound7, z_Wound30), "Figure2", "A", "Discovery cohort frozen stage scores")

# 2B disease projection into frozen reference score space
ref_matrix <- ref_fib |> dplyr::select(z_Skin, z_Wound1, z_Wound7, z_Wound30)
pca_projection <- stats::prcomp(ref_matrix, center = TRUE, scale. = TRUE)
ref_proj <- dplyr::bind_cols(
  ref_fib |> dplyr::select(sample_id, donor, true_stage),
  as.data.frame(pca_projection$x[,1:2, drop = FALSE])
) |>
  dplyr::mutate(type = "Reference", group = as.character(true_stage), id = sample_id)

dis_proj_coords <- predict(
  pca_projection,
  newdata = s5_fib |> dplyr::select(z_Skin, z_Wound1, z_Wound7, z_Wound30)
)[,1:2, drop = FALSE]

dis_proj <- dplyr::bind_cols(
  s5_fib |> dplyr::select(patient_id, group),
  as.data.frame(dis_proj_coords)
) |>
  dplyr::mutate(type = "Disease cohort", id = patient_id)

ref_stage_centers <- ref_proj |>
  dplyr::group_by(true_stage) |>
  dplyr::summarise(
    PC1 = mean(PC1, na.rm = TRUE),
    PC2 = mean(PC2, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(stage_label = dplyr::recode(as.character(true_stage), !!!stage_labels))

write_source(ref_proj, "Figure2", "B_reference", "Reference PCA coordinates used as the background manifold for disease projection")
write_source(dis_proj, "Figure2", "B_disease", "Discovery patient coordinates projected into the frozen reference score space")
write_source(ref_stage_centers, "Figure2", "B_reference_centers", "Reference stage centroids used to orient the disease projection")

p2b <- ggplot2::ggplot() +
  ggplot2::geom_path(
    data = ref_proj,
    ggplot2::aes(PC1, PC2, group = donor),
    linewidth = 0.55,
    colour = "#B8B8B8",
    alpha = 0.80,
    arrow = grid::arrow(length = grid::unit(0.06, "inches"))
  ) +
  ggplot2::geom_point(
    data = ref_proj,
    ggplot2::aes(PC1, PC2),
    size = 1.8,
    shape = 16,
    colour = "#9AA3AD",
    alpha = 0.65
  ) +
  ggrepel::geom_text_repel(
    data = ref_stage_centers,
    ggplot2::aes(PC1, PC2, label = stage_label),
    size = 2.6,
    family = BASE_FONT,
    colour = "#6F7780",
    fontface = "bold",
    segment.size = 0.25,
    box.padding = 0.20,
    point.padding = 0.15,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123
  ) +
  ggplot2::geom_point(
    data = dis_proj,
    ggplot2::aes(PC1, PC2, fill = group),
    shape = 21,
    size = 5.0,
    stroke = 0.95,
    colour = "black"
  ) +
  ggrepel::geom_text_repel(
    data = dis_proj,
    ggplot2::aes(PC1, PC2, label = patient_id),
    size = 3.0,
    family = BASE_FONT,
    fontface = "bold",
    segment.size = 0.35,
    box.padding = 0.28,
    point.padding = 0.22,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 123
  ) +
  ggplot2::scale_fill_manual(
    values = group_colors[c("normal_scar", "keloid")],
    breaks = c("normal_scar", "keloid"),
    labels = c("Normal scar", "Keloid")
  ) +
  ggplot2::labs(
    title = "Disease projection onto wound reference",
    x = "Reference PC1",
    y = "Reference PC2",
    fill = NULL
  ) +
  theme_pub() +
  ggplot2::theme(
    legend.position = "top",
    panel.grid.minor = ggplot2::element_blank()
  )

# 2C stage weights
s5_bar <- s5_fib |> dplyr::arrange(factor(patient_id, levels = s5_order))
p2c <- stage_probability_bars(s5_bar, "patient_id", "group", "Stage weights by patient")
write_source(s5_bar |> dplyr::select(patient_id, group, dplyr::starts_with("weight_")), "Figure2", "C", "Discovery patient reference-stage weights")

# 2D primary endpoint
p2d <- patient_dot_plot(s5_fib, "late_remodeling_state", "Persistent late remodeling", group_col = "group", id_col = "patient_id")
write_source(s5_fib |> dplyr::select(patient_id, group, late_remodeling_state), "Figure2", "D", "Discovery patient late-remodeling state")

# 2E Skin-Wound30 plane
stage_centroids <- ref_fib |>
  dplyr::group_by(true_stage) |>
  dplyr::summarise(z_Skin = mean(z_Skin), z_Wound30 = mean(z_Wound30), .groups = "drop")
write_source(stage_centroids, "Figure2", "E_reference", "Reference-stage centroids for the Skin-Day30 plane")
write_source(s5_fib |> dplyr::select(patient_id, group, z_Skin, z_Wound30), "Figure2", "E_disease", "Discovery patients in the Skin-Day30 plane")
p2e <- ggplot2::ggplot() +
  ggplot2::geom_path(
    data = stage_centroids, ggplot2::aes(z_Skin, z_Wound30, group = 1),
    linewidth = 0.65, colour = "#777777", arrow = grid::arrow(length = grid::unit(0.08, "inches"))
  ) +
  ggplot2::geom_point(
    data = stage_centroids, ggplot2::aes(z_Skin, z_Wound30, colour = true_stage), size = 3.0
  ) +
  ggrepel::geom_text_repel(
    data = stage_centroids, ggplot2::aes(z_Skin, z_Wound30, label = stage_labels[as.character(true_stage)]),
    size = 2.7, family = BASE_FONT, segment.size = 0.3
  ) +
  ggplot2::geom_point(
    data = s5_fib, ggplot2::aes(z_Skin, z_Wound30, fill = group),
    shape = 21, size = 4.0, stroke = 0.75, colour = "black"
  ) +
  ggrepel::geom_text_repel(
    data = s5_fib, ggplot2::aes(z_Skin, z_Wound30, label = patient_id),
    size = 2.55, family = BASE_FONT, segment.size = 0.3, max.overlaps = Inf
  ) +
  ggplot2::scale_colour_manual(values = stage_colors, guide = "none") +
  ggplot2::scale_fill_manual(values = group_colors, labels = c("Normal scar", "Keloid")) +
  ggplot2::labs(title = "Skin loss with Day 30 persistence", x = "Skin-like score", y = "Day-30-like score", fill = NULL) +
  theme_pub() + ggplot2::theme(legend.position = "top")

# 2F effect summary
selected_s5_metrics <- c("remodeling_completion", "ordinal_position", "z_Wound30", "wound_activation", "z_Skin")
s5_forest <- s5_effects |>
  dplyr::filter(cell_type == "Fibroblast", metric %in% selected_s5_metrics) |>
  dplyr::mutate(
    metric_label = pretty_metric(metric),
    oriented_g = ifelse(metric == "z_Skin", -hedges_g, hedges_g),
    metric_label = factor(metric_label, levels = rev(pretty_metric(selected_s5_metrics)))
  )
write_source(s5_forest, "Figure2", "F", "Discovery cohort patient-level effect estimates")
p2f <- ggplot2::ggplot(s5_forest, ggplot2::aes(oriented_g, metric_label)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "#8A8A8A") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = oriented_g, yend = metric_label), linewidth = 0.65, colour = "#A7A7A7") +
  ggplot2::geom_point(size = 3.1, colour = "#B2182B") +
  ggplot2::labs(title = "Patient-level effect sizes", x = "Oriented Hedges g (positive supports keloid persistence)", y = NULL) +
  theme_pub()

fig2_design <- "
AABBB
CCCDD
EEEFF
"
figure2 <- p2a + p2b + p2c + p2d + p2e + p2f +
  patchwork::plot_layout(design = fig2_design, guides = "keep") +
  patchwork::plot_annotation(tag_levels = "A", theme = theme_patch_tags)
if (SAVE_COMPOSITE_TIFFS) save_tiff(figure2, "Figure2.tif", width = 7.2, height = 8.8)
save_panel_set(
  list(A = p2a, B = p2b, C = p2c, D = p2d, E = p2e, F = p2f),
  "Figure2",
  width = c(6.4, 6.4, 5.8, 6.0, 6.0, 5.8),
  height = c(5.0, 5.0, 4.4, 4.8, 4.8, 4.6)
)

# =============================== FIGURE 3 ====================================
message("Building Figure 3...")

s6_order <- s6_primary |>
  dplyr::arrange(group, dplyr::desc(late_remodeling_state)) |>
  dplyr::pull(patient_id)

p3a <- heatmap_stage_scores(s6_primary, "patient_id", "group", "Frozen stage scores in the independent cohort", sample_order = s6_order)
write_source(s6_primary |> dplyr::select(patient_id, group, z_Skin, z_Wound1, z_Wound7, z_Wound30), "Figure3", "A", "Independent replication cohort frozen stage scores")

p3b <- stage_probability_bars(s6_primary |> dplyr::arrange(factor(patient_id, levels = s6_order)), "patient_id", "group", "Reference-stage weights after frozen mapping")
write_source(s6_primary |> dplyr::select(patient_id, group, dplyr::starts_with("weight_")), "Figure3", "B", "Independent replication reference-stage weights")

p3c <- patient_dot_plot(s6_primary, "late_remodeling_state", "Independent replication of the primary state", group_col = "group", id_col = "patient_id")
write_source(s6_primary |> dplyr::select(patient_id, group, late_remodeling_state), "Figure3", "C", "Independent replication patient late-remodeling values")

core_late <- core_patient |>
  dplyr::select(cohort, patient_id, group, late_remodeling_state) |>
  dplyr::mutate(
    group_label = factor(
      pretty_group(as.character(group)),
      levels = c("Normal scar", "Keloid")
    )
  )
write_source(core_late, "Figure3", "D", "Patient-level primary endpoint in both core cohorts")
p3d <- ggplot2::ggplot(core_late, ggplot2::aes(group_label, late_remodeling_state, fill = group)) +
  ggplot2::geom_point(shape = 21, size = 3.0, stroke = 0.55, colour = "black",
                      position = ggplot2::position_jitter(width = 0.07, seed = 20260618)) +
  ggplot2::stat_summary(fun = mean, geom = "point", shape = 95, size = 8, colour = "black") +
  ggplot2::facet_wrap(~cohort, nrow = 1) +
  ggplot2::scale_fill_manual(values = group_colors, guide = "none") +
  ggplot2::labs(title = "Patient-level replication across core cohorts", x = NULL, y = "Late-remodeling state") +
  theme_pub()

rep_forest <- replication_table |>
  dplyr::filter(validation_analysis_set == "high_specificity") |>
  dplyr::filter(metric %in% c("late_remodeling_state", "ordinal_wound_state_position", "z_Wound30", "wound_activation", "z_Skin")) |>
  dplyr::select(metric, discovery_hedges_g, validation_hedges_g) |>
  tidyr::pivot_longer(c(discovery_hedges_g, validation_hedges_g), names_to = "cohort", values_to = "hedges_g") |>
  dplyr::mutate(
    cohort = dplyr::recode(cohort, discovery_hedges_g = "GSE163973", validation_hedges_g = "GSE181316"),
    oriented_g = ifelse(metric == "z_Skin", -hedges_g, hedges_g),
    metric_label = pretty_metric(metric),
    metric_label = factor(metric_label, levels = rev(unique(metric_label)))
  )
write_source(rep_forest, "Figure3", "E", "Discovery and validation Hedges g for the frozen endpoints")
p3e <- ggplot2::ggplot(rep_forest, ggplot2::aes(oriented_g, metric_label, colour = cohort)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "#888888") +
  ggplot2::geom_line(ggplot2::aes(group = metric_label), colour = "#B5B5B5", linewidth = 0.55) +
  ggplot2::geom_point(size = 3.0) +
  ggplot2::scale_colour_manual(values = cohort_colors) +
  ggplot2::labs(title = "Replication effect sizes", x = "Oriented Hedges g", y = NULL, colour = NULL) +
  theme_pub() + ggplot2::theme(legend.position = "top")

perm_late <- perm_distrib |> dplyr::filter(metric == "late_remodeling_state")
late_integrated_row <- core_integrated |>
  dplyr::filter(metric == "late_remodeling_state")
obs_late <- late_integrated_row$equal_weight_mean_difference[[1]]
obs_late_p <- late_integrated_row$exact_stratified_two_sided_p[[1]]
write_source(perm_late, "Figure3", "F", "Exact stratified null distribution for the primary endpoint")
p3f <- ggplot2::ggplot(perm_late, ggplot2::aes(null_equal_weight_mean_difference)) +
  ggplot2::geom_histogram(bins = 28, fill = "#D9D9D9", colour = "white", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = obs_late, colour = "#B2182B", linewidth = 0.9) +
  ggplot2::annotate(
    "text",
    x = obs_late,
    y = Inf,
    label = paste0(
      "Observed = ",
      sprintf("%.3f", obs_late),
      "\nExact P = ",
      format(obs_late_p, digits = 3)
    ),
    vjust = 1.3,
    hjust = 1.05,
    colour = "#B2182B",
    family = BASE_FONT,
    size = 2.8
  ) +
  ggplot2::labs(title = "Exact integration", x = "Null equal-weight mean difference", y = "Permutation count") +
  theme_pub()

integrated_plot <- core_integrated |>
  dplyr::filter(
    metric %in% c(
      "late_remodeling_state",
      "ordinal_wound_state_position",
      "z_Wound30",
      "wound_activation",
      "z_Skin"
    )
  ) |>
  dplyr::mutate(
    metric_label = pretty_metric(metric),
    oriented_effect = ifelse(
      expected_direction < 0,
      -equal_weight_mean_difference,
      equal_weight_mean_difference
    ),
    endpoint_role = factor(
      endpoint_role,
      levels = c("primary", "key_secondary", "supporting")
    ),
    metric_label = factor(
      metric_label,
      levels = rev(pretty_metric(metric))
    )
  )
write_source(integrated_plot, "Figure3", "G", "Stratified exact integration of core endpoints")
p3g <- ggplot2::ggplot(integrated_plot, ggplot2::aes(oriented_effect, metric_label)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "#888888") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = oriented_effect, yend = metric_label), linewidth = 0.7, colour = "#A8A8A8") +
  ggplot2::geom_point(ggplot2::aes(fill = endpoint_role), shape = 21, size = 3.3, stroke = 0.5, colour = "black") +
  ggplot2::geom_text(ggplot2::aes(label = paste0("P=", format(exact_stratified_two_sided_p, trim = TRUE))), hjust = -0.1, size = 2.45, family = BASE_FONT) +
  ggplot2::scale_fill_manual(
    values = c(
      primary = "#B2182B",
      key_secondary = "#7B3294",
      supporting = "#2166AC"
    ),
    guide = "none"
  ) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.35))) +
  ggplot2::labs(title = "Integrated core endpoints", x = "Oriented equal-weight mean difference", y = NULL) +
  theme_pub()

k3_long <- k3_lr |>
  dplyr::filter(analysis_set == "high_specificity") |>
  tidyr::pivot_longer(c(late_remodeling_state, ordinal_wound_state_position, z_Wound30, z_Skin), names_to = "metric", values_to = "value") |>
  dplyr::mutate(
    side = dplyr::case_when(stringr::str_detect(sample_key, "L$") ~ "Left", stringr::str_detect(sample_key, "R$") ~ "Right", TRUE ~ sample_key),
    metric_label = pretty_metric(metric)
  )
write_source(k3_long, "Figure3", "H", "Within-patient left-right concordance for K3")
p3h <- ggplot2::ggplot(k3_long, ggplot2::aes(side, value, group = metric_label, colour = metric_label)) +
  ggplot2::geom_line(linewidth = 0.65) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::facet_wrap(~metric_label, scales = "free_y", ncol = 2) +
  ggplot2::labs(title = "K3 left-right consistency", x = NULL, y = NULL, colour = NULL) +
  theme_pub() + ggplot2::theme(legend.position = "none")

fig3_design <- "
AAABB
CCDDD
EEEFF
GGGHH
"
figure3 <- p3a + p3b + p3c + p3d + p3e + p3f + p3g + p3h +
  patchwork::plot_layout(design = fig3_design, guides = "keep") +
  patchwork::plot_annotation(tag_levels = "A", theme = theme_patch_tags)
if (SAVE_COMPOSITE_TIFFS) save_tiff(figure3, "Figure3.tif", width = 7.2, height = 9.8)
save_panel_set(
  list(A = p3a, B = p3b, C = p3c, D = p3d, E = p3e, F = p3f, G = p3g, H = p3h),
  "Figure3",
  width = c(6.4, 5.8, 5.8, 6.0, 5.8, 6.0, 6.2, 5.6),
  height = c(5.0, 4.4, 4.6, 4.8, 4.6, 4.8, 5.0, 4.4)
)

# =============================== FIGURE 4 ====================================
message("Building Figure 4...")

s7_high <- s7_patient |> dplyr::filter(analysis_set == "high_specificity")
s7_order <- s7_high |> dplyr::arrange(aggregate_state, patient_id) |> dplyr::pull(sample_key)
p4a <- heatmap_stage_scores(s7_high, "sample_key", title = "Frozen stage scores in the external keloid cohort", sample_order = s7_order)
write_source(s7_high |> dplyr::select(sample_key, patient_id, aggregate_state, z_Skin, z_Wound1, z_Wound7, z_Wound30), "Figure4", "A", "GSE220300 high-specificity patient-level stage scores")

s7_support_group <- s7_high |>
  dplyr::mutate(
    support_group = dplyr::case_when(
      aggregate_state %in% c("active_lesion", "inactive_lesion") ~ "keloid_lesion",
      aggregate_state %in% c("mature_scar", "normal_scar") ~ "scar_reference",
      TRUE ~ NA_character_
    ),
    support_group = factor(
      support_group,
      levels = c("scar_reference", "keloid_lesion")
    ),
    support_label = factor(
      dplyr::recode(
        as.character(support_group),
        scar_reference = "Scar reference",
        keloid_lesion = "Keloid lesion"
      ),
      levels = c("Scar reference", "Keloid lesion")
    )
  ) |>
  dplyr::filter(!is.na(support_group))

write_source(
  s7_support_group,
  "Figure4",
  "B",
  "Patient-level lesion-versus-scar support in GSE220300"
)

p4b <- ggplot2::ggplot(
  s7_support_group,
  ggplot2::aes(support_label, late_remodeling_state, fill = support_group)
) +
  ggplot2::geom_point(
    shape = 21,
    size = 3.2,
    stroke = 0.55,
    colour = "black",
    position = ggplot2::position_jitter(
      width = 0.07,
      seed = 20260618
    )
  ) +
  ggplot2::stat_summary(
    fun = mean,
    geom = "point",
    shape = 95,
    size = 8,
    colour = "black"
  ) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = sample_key),
    size = 2.4,
    family = BASE_FONT,
    segment.size = 0.3,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      scar_reference = "#5A6268",
      keloid_lesion = "#B2182B"
    ),
    guide = "none"
  ) +
  ggplot2::labs(
    title = "Patient-level lesion-versus-scar support",
    x = NULL,
    y = "Late-remodeling state"
  ) +
  theme_pub()

cross_support <- s7_cross |>
  dplyr::filter(metric %in% c("late_remodeling_state", "ordinal_wound_state_position", "z_Wound30", "z_Skin")) |>
  dplyr::filter(analysis_set %in% c("author_annotated_fibroblast", "high_specificity", "high_specificity_primary")) |>
  dplyr::mutate(
    oriented_g = ifelse(metric == "z_Skin", -hedges_g, hedges_g),
    metric_label = pretty_metric(metric),
    cohort = factor(cohort, levels = c("GSE163973", "GSE181316", "GSE220300"))
  )
write_source(cross_support, "Figure4", "C", "Cross-cohort direction summary including GSE220300")
p4c <- ggplot2::ggplot(cross_support, ggplot2::aes(oriented_g, metric_label, colour = cohort)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "#8A8A8A") +
  ggplot2::geom_point(size = 3.0, position = ggplot2::position_dodge(width = 0.45)) +
  ggplot2::scale_colour_manual(values = cohort_colors, na.translate = FALSE) +
  ggplot2::labs(title = "Direction persists in a third disease cohort", x = "Oriented Hedges g", y = NULL, colour = NULL) +
  theme_pub() + ggplot2::theme(legend.position = "top")

sp_conf <- s8_primary |>
  dplyr::count(condition, nearest_stage, name = "n") |>
  dplyr::group_by(condition) |>
  dplyr::mutate(row_percent = n / sum(n)) |>
  dplyr::ungroup()
write_source(sp_conf, "Figure4", "D", "Spatial section-level confusion matrix")
p4d <- confusion_matrix_clean(
  data = sp_conf,
  x_col = "nearest_stage",
  y_col = "condition",
  count_col = "n",
  fill_col = "row_percent",
  title = "Spatial stage classification",
  x_lab = "Predicted stage",
  y_lab = "Actual section"
)

write_source(s8_primary, "Figure4", "E", "Spatial section wound-state positions by donor")
p4e <- ggplot2::ggplot(s8_primary, ggplot2::aes(stage_code, ordinal_wound_state_position, group = donor, colour = donor)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.45, colour = "#999999") +
  ggplot2::geom_line(linewidth = 0.65) +
  ggplot2::geom_point(size = 2.8) +
  ggplot2::scale_x_continuous(breaks = 0:3, labels = stage_labels) +
  ggplot2::labs(title = "Ordinal trajectories in four new donors", x = "Actual stage", y = "Predicted ordinal position", colour = "Donor") +
  theme_pub() + ggplot2::theme(legend.position = "top", axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

# Representative coordinate maps: choose first donor with all 4 stages
rep_donors <- s8_spot_preview |> dplyr::distinct(donor, condition) |> dplyr::count(donor) |> dplyr::filter(n >= 4) |> dplyr::pull(donor)
rep_donor <- ifelse(length(rep_donors) > 0, rep_donors[[1]], unique(s8_spot_preview$donor)[[1]])
spot_rep <- s8_spot_preview |>
  dplyr::filter(donor == rep_donor, high_fibroblast_enriched) |>
  dplyr::mutate(condition = factor(condition, levels = stage_levels))
write_source(spot_rep, "Figure4", "F", paste0("Downsampled representative fibroblast-enriched spatial spots from ", rep_donor))
p4f <- ggplot2::ggplot(spot_rep, ggplot2::aes(pxl_col_in_fullres, -pxl_row_in_fullres, colour = spot_late_remodeling)) +
  ggplot2::geom_point(size = 0.55, alpha = 0.9) +
  ggplot2::facet_wrap(~condition, nrow = 1, labeller = ggplot2::as_labeller(stage_labels)) +
  ggplot2::scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, name = "Late-remodeling\nspot score") +
  ggplot2::coord_equal() +
  ggplot2::labs(title = paste0("Representative spatial maps (", rep_donor, "; display subset)"), x = NULL, y = NULL) +
  theme_pub() +
  ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.line = ggplot2::element_blank(), legend.position = "right")

stage_sp <- s8_stage |>
  dplyr::filter(analysis_set == "high_specificity_spots") |>
  dplyr::select(condition, mean_late_remodeling, mean_skin_return) |>
  tidyr::pivot_longer(-condition, names_to = "metric", values_to = "mean_score") |>
  dplyr::mutate(
    condition = factor(condition, levels = stage_levels),
    metric = dplyr::recode(metric, mean_late_remodeling = "Late-remodeling state", mean_skin_return = "Skin-return state")
  )
write_source(stage_sp, "Figure4", "G", "Stage-level spatial summary")
p4g <- ggplot2::ggplot(stage_sp, ggplot2::aes(condition, mean_score, colour = metric, group = metric)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#A0A0A0") +
  ggplot2::geom_line(linewidth = 0.75) +
  ggplot2::geom_point(size = 3.0) +
  ggplot2::scale_x_discrete(labels = stage_labels) +
  ggplot2::scale_colour_manual(values = c("Late-remodeling state" = "#B2182B", "Skin-return state" = "#2166AC")) +
  ggplot2::labs(title = "Late-remodeling versus skin return", x = NULL, y = "Mean section score", colour = NULL) +
  theme_pub() + ggplot2::theme(legend.position = "top", axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

s8_donor_primary <- s8_donor |> dplyr::filter(analysis_set == "high_specificity_spots") |>
  dplyr::select(donor, exact_accuracy, adjacent_accuracy, mean_ordinal_error, spearman) |>
  tidyr::pivot_longer(-donor, names_to = "metric", values_to = "value") |>
  dplyr::mutate(metric = dplyr::recode(metric,
    exact_accuracy = "Exact accuracy", adjacent_accuracy = "Adjacent accuracy",
    mean_ordinal_error = "Mean ordinal error", spearman = "Spearman rho"
  ))
write_source(s8_donor_primary, "Figure4", "H", "Spatial donor-level validation metrics")
p4h <- ggplot2::ggplot(s8_donor_primary, ggplot2::aes(donor, value, group = 1)) +
  ggplot2::geom_line(linewidth = 0.5, colour = "#666666") +
  ggplot2::geom_point(size = 2.7, shape = 21, fill = "white", stroke = 0.6) +
  ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 2) +
  ggplot2::labs(title = "Spatial donor metrics", x = NULL, y = NULL) +
  theme_pub() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

fig4_design <- "
AAABB
CCDDD
EEEEE
FFFFF
GGGHH
"
figure4 <- p4a + p4b + p4c + p4d + p4e + p4f + p4g + p4h +
  patchwork::plot_layout(design = fig4_design, guides = "keep") +
  patchwork::plot_annotation(tag_levels = "A", theme = theme_patch_tags)
if (SAVE_COMPOSITE_TIFFS) save_tiff(figure4, "Figure4.tif", width = 7.2, height = 10.4)
save_panel_set(
  list(A = p4a, B = p4b, C = p4c, D = p4d, E = p4e, F = p4f, G = p4g, H = p4h),
  "Figure4",
  width = c(6.4, 5.8, 5.8, 5.9, 6.0, 6.6, 6.0, 5.8),
  height = c(5.0, 4.6, 4.6, 5.0, 5.0, 5.8, 4.8, 4.4)
)

# =============================== FIGURE 5 ====================================
message("Building Figure 5...")

s9_order <- s9_primary |> dplyr::arrange(group, dplyr::desc(early_state_persistence)) |> dplyr::pull(sample_key)
p5a <- heatmap_stage_scores(s9_primary, "sample_key", "group", "Frozen wound-state scores in venous ulcer", sample_order = s9_order)
write_source(s9_primary |> dplyr::select(sample_key, group, z_Skin, z_Wound1, z_Wound7, z_Wound30), "Figure5", "A", "Venous ulcer primary-analysis frozen stage scores")

p5b <- stage_probability_bars(s9_primary |> dplyr::arrange(factor(sample_key, levels = s9_order)), "sample_key", "group", "Reference-stage weights in venous ulcer")
write_source(s9_primary |> dplyr::select(sample_key, group, dplyr::starts_with("weight_")), "Figure5", "B", "Venous ulcer reference-stage weights")

vu_metrics <- s9_primary |>
  dplyr::select(sample_key, group, early_state_persistence, late_remodeling_state, late_vs_early_balance) |>
  tidyr::pivot_longer(c(early_state_persistence, late_remodeling_state, late_vs_early_balance), names_to = "metric", values_to = "value") |>
  dplyr::mutate(
    metric_label = pretty_metric(metric),
    group_label = factor(
      pretty_group(as.character(group)),
      levels = c("Normal skin", "Venous ulcer")
    )
  )
write_source(vu_metrics, "Figure5", "C", "Patient-level venous ulcer state metrics")
p5c <- ggplot2::ggplot(vu_metrics, ggplot2::aes(group_label, value, fill = group)) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#999999") +
  ggplot2::geom_point(shape = 21, size = 2.9, stroke = 0.5, colour = "black",
                      position = ggplot2::position_jitter(width = 0.07, seed = 20260618)) +
  ggplot2::stat_summary(fun = mean, geom = "point", shape = 95, size = 7, colour = "black") +
  ggplot2::facet_wrap(~metric_label, scales = "free_y", nrow = 1) +
  ggplot2::scale_fill_manual(values = group_colors, guide = "none") +
  ggplot2::labs(title = "Venous ulcer preferentially retains early-state activity", x = NULL, y = "State score") +
  theme_pub()

failure_plane_plot <- failure_plane |>
  dplyr::mutate(
    disease_label = dplyr::recode(disease,
      Keloid_core_integrated = "Keloid",
      Venous_ulcer = "Venous ulcer",
      Venous_ulcer_primary = "Venous ulcer",
      venous_ulcer = "Venous ulcer",
      .default = disease
    )
  )
write_source(failure_plane_plot, "Figure5", "D", "Aggregated cross-disease early versus late failure-mode effects")
p5d <- ggplot2::ggplot(failure_plane_plot, ggplot2::aes(early_state_effect, late_remodeling_effect, colour = disease)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.4, colour = "#999999") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.4, colour = "#999999") +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dotted", linewidth = 0.55, colour = "#777777") +
  ggplot2::geom_point(size = 4.2) +
  ggrepel::geom_text_repel(
    ggplot2::aes(label = disease_label),
    size = 3.0, family = BASE_FONT, segment.size = 0.35
  ) +
  ggplot2::scale_colour_manual(values = c(Keloid_core_integrated = "#B2182B", Venous_ulcer = "#7B3294", Venous_ulcer_primary = "#7B3294", venous_ulcer = "#7B3294"), guide = "none") +
  ggplot2::annotate(
    "text",
    x = min(failure_plane_plot$early_state_effect),
    y = max(failure_plane_plot$late_remodeling_effect),
    label = "Late-remodeling dominant",
    hjust = 0,
    vjust = 1.1,
    family = BASE_FONT,
    size = 2.6,
    colour = "#666666"
  ) +
  ggplot2::annotate(
    "text",
    x = max(failure_plane_plot$early_state_effect),
    y = min(failure_plane_plot$late_remodeling_effect),
    label = "Early-state dominant",
    hjust = 1,
    vjust = -0.1,
    family = BASE_FONT,
    size = 2.6,
    colour = "#666666"
  ) +
  ggplot2::labs(title = "Distinct wound-resolution failure modes", x = "Early-state persistence effect", y = "Late-remodeling effect") +
  theme_pub()

cross_disease_effects <- dplyr::bind_rows(
  core_integrated |>
    dplyr::select(metric, mean_difference = equal_weight_mean_difference) |>
    dplyr::mutate(disease = "Keloid"),
  s9_effects |>
    dplyr::filter(analysis_set == "FB1_FB2_combined_primary") |>
    dplyr::select(metric, mean_difference) |>
    dplyr::mutate(disease = "Venous ulcer")
) |>
  dplyr::filter(metric %in% c("late_remodeling_state", "early_state_persistence", "ordinal_wound_state_position", "z_Wound30", "z_Skin", "wound_activation")) |>
  dplyr::mutate(
    oriented_effect = ifelse(metric == "z_Skin", -mean_difference, mean_difference),
    metric_label = pretty_metric(metric),
    disease = factor(disease, levels = c("Keloid", "Venous ulcer"))
  )
write_source(cross_disease_effects, "Figure5", "E", "Cross-disease comparison of oriented wound-state effects")
p5e <- ggplot2::ggplot(cross_disease_effects, ggplot2::aes(oriented_effect, metric_label, colour = disease)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "#888888") +
  ggplot2::geom_line(ggplot2::aes(group = metric_label), colour = "#C0C0C0", linewidth = 0.55) +
  ggplot2::geom_point(size = 3.1) +
  ggplot2::scale_colour_manual(values = c("Keloid" = "#B2182B", "Venous ulcer" = "#7B3294")) +
  ggplot2::labs(title = "The contrast extends beyond a single composite score", x = "Oriented disease-control mean difference", y = NULL, colour = NULL) +
  theme_pub() + ggplot2::theme(legend.position = "top")

patient_failure_plane <- dplyr::bind_rows(
  core_patient |>
    dplyr::transmute(
      sample_id = patient_id,
      cohort = cohort,
      state_group = dplyr::if_else(
        as.character(group) == "keloid",
        "Keloid",
        "Normal scar"
      ),
      early_state_persistence = early_state_persistence,
      late_remodeling_state = late_remodeling_state
    ),
  s9_primary |>
    dplyr::transmute(
      sample_id = sample_key,
      cohort = "GSE265972",
      state_group = dplyr::if_else(
        as.character(group) == "venous_ulcer",
        "Venous ulcer",
        "Normal skin"
      ),
      early_state_persistence = early_state_persistence,
      late_remodeling_state = late_remodeling_state
    )
) |>
  dplyr::mutate(
    state_group = factor(
      state_group,
      levels = c(
        "Normal scar",
        "Keloid",
        "Normal skin",
        "Venous ulcer"
      )
    )
  )

write_source(
  patient_failure_plane,
  "Figure5",
  "F",
  "Patient-level early-state and late-remodeling positions across keloid and venous-ulcer cohorts"
)

# Figure 5F is a single-panel state plane. It should not use a heavy
# top legend or dataset-specific shape coding, because that caused clipping
# and over-encoding. Disease groups are emphasized, controls are subdued, and
# group labels are placed inside an expanded coordinate window.
p5f_range_x <- range(
  patient_failure_plane$early_state_persistence,
  na.rm = TRUE
)
p5f_range_y <- range(
  patient_failure_plane$late_remodeling_state,
  na.rm = TRUE
)

p5f_span_x <- diff(
  p5f_range_x
)
p5f_span_y <- diff(
  p5f_range_y
)

if (
  !is.finite(p5f_span_x) ||
    p5f_span_x == 0
) {
  p5f_span_x <- 1
}

if (
  !is.finite(p5f_span_y) ||
    p5f_span_y == 0
) {
  p5f_span_y <- 1
}

p5f_xlim <- c(
  p5f_range_x[1] - 0.22 * p5f_span_x,
  p5f_range_x[2] + 0.34 * p5f_span_x
)

p5f_ylim <- c(
  p5f_range_y[1] - 0.22 * p5f_span_y,
  p5f_range_y[2] + 0.28 * p5f_span_y
)

p5f_label_df <- patient_failure_plane |>
  dplyr::group_by(
    state_group
  ) |>
  dplyr::summarise(
    x = mean(
      early_state_persistence,
      na.rm = TRUE
    ),
    y = mean(
      late_remodeling_state,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    state_group_chr = as.character(
      state_group
    ),
    dx = dplyr::case_when(
      state_group_chr == "Keloid" ~ 0.10 * p5f_span_x,
      state_group_chr == "Venous ulcer" ~ 0.10 * p5f_span_x,
      state_group_chr == "Normal scar" ~ -0.08 * p5f_span_x,
      state_group_chr == "Normal skin" ~ -0.08 * p5f_span_x,
      TRUE ~ 0
    ),
    dy = dplyr::case_when(
      state_group_chr == "Keloid" ~ 0.10 * p5f_span_y,
      state_group_chr == "Venous ulcer" ~ -0.10 * p5f_span_y,
      state_group_chr == "Normal scar" ~ 0.08 * p5f_span_y,
      state_group_chr == "Normal skin" ~ -0.08 * p5f_span_y,
      TRUE ~ 0
    ),
    label_x = pmin(
      pmax(
        x + dx,
        p5f_xlim[1] + 0.06 * diff(p5f_xlim)
      ),
      p5f_xlim[2] - 0.06 * diff(p5f_xlim)
    ),
    label_y = pmin(
      pmax(
        y + dy,
        p5f_ylim[1] + 0.06 * diff(p5f_ylim)
      ),
      p5f_ylim[2] - 0.06 * diff(p5f_ylim)
    ),
    label_colour = dplyr::case_when(
      state_group_chr == "Keloid" ~ "#B2182B",
      state_group_chr == "Venous ulcer" ~ "#7B3294",
      TRUE ~ "#4F565C"
    )
  )

p5f_controls <- patient_failure_plane |>
  dplyr::filter(
    state_group %in% c(
      "Normal scar",
      "Normal skin"
    )
  )

p5f_diseases <- patient_failure_plane |>
  dplyr::filter(
    state_group %in% c(
      "Keloid",
      "Venous ulcer"
    )
  )

p5f_centroids <- patient_failure_plane |>
  dplyr::group_by(
    state_group
  ) |>
  dplyr::summarise(
    early_state_persistence = mean(
      early_state_persistence,
      na.rm = TRUE
    ),
    late_remodeling_state = mean(
      late_remodeling_state,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

p5f <- ggplot2::ggplot() +
  ggplot2::annotate(
    "text",
    x = p5f_xlim[1] + 0.22 * diff(p5f_xlim),
    y = p5f_ylim[2] - 0.08 * diff(p5f_ylim),
    label = "Late-remodeling\ndominant",
    family = BASE_FONT,
    size = 3.0,
    lineheight = 0.92,
    colour = "#777777"
  ) +
  ggplot2::annotate(
    "text",
    x = p5f_xlim[2] - 0.22 * diff(p5f_xlim),
    y = p5f_ylim[1] + 0.10 * diff(p5f_ylim),
    label = "Early-state\ndominant",
    family = BASE_FONT,
    size = 3.0,
    lineheight = 0.92,
    colour = "#777777"
  ) +
  ggplot2::geom_hline(
    yintercept = 0,
    linewidth = 0.35,
    colour = "#B0B0B0"
  ) +
  ggplot2::geom_vline(
    xintercept = 0,
    linewidth = 0.35,
    colour = "#B0B0B0"
  ) +
  ggplot2::geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dotted",
    linewidth = 0.45,
    colour = "#8C8C8C"
  ) +
  ggplot2::geom_point(
    data = p5f_controls,
    ggplot2::aes(
      early_state_persistence,
      late_remodeling_state,
      fill = state_group
    ),
    shape = 21,
    size = 3.2,
    stroke = 0.45,
    colour = "#4D4D4D",
    alpha = 0.68
  ) +
  ggplot2::geom_point(
    data = p5f_diseases,
    ggplot2::aes(
      early_state_persistence,
      late_remodeling_state,
      fill = state_group
    ),
    shape = 21,
    size = 4.1,
    stroke = 0.70,
    colour = "black",
    alpha = 0.95
  ) +
  ggplot2::geom_point(
    data = p5f_centroids,
    ggplot2::aes(
      early_state_persistence,
      late_remodeling_state,
      fill = state_group
    ),
    shape = 23,
    size = 5.0,
    stroke = 0.85,
    colour = "black"
  ) +
  ggplot2::geom_label(
    data = p5f_label_df,
    ggplot2::aes(
      label_x,
      label_y,
      label = state_group_chr,
      colour = state_group_chr
    ),
    family = BASE_FONT,
    fontface = "bold",
    size = 3.0,
    label.size = 0.20,
    label.r = grid::unit(2.0, "pt"),
    label.padding = grid::unit(2.0, "pt"),
    fill = "white"
  ) +
  ggplot2::scale_colour_manual(
    values = c(
      "Normal scar" = "#4F565C",
      "Keloid" = "#B2182B",
      "Normal skin" = "#4F565C",
      "Venous ulcer" = "#7B3294"
    ),
    guide = "none"
  ) +
  ggplot2::scale_fill_manual(
    values = c(
      "Normal scar" = "#9EA5AA",
      "Keloid" = "#B2182B",
      "Normal skin" = "#C4C4C4",
      "Venous ulcer" = "#7B3294"
    ),
    guide = "none"
  ) +
  ggplot2::coord_cartesian(
    xlim = p5f_xlim,
    ylim = p5f_ylim,
    clip = "off"
  ) +
  ggplot2::labs(
    title = "Wound-state failure plane",
    x = "Early-state score",
    y = "Late-remodeling score"
  ) +
  theme_pub(base_size = 9.2, title_size = 9.5) +
  ggplot2::theme(
    legend.position = "none",
    plot.margin = ggplot2::margin(12, 28, 14, 16),
    axis.title = ggplot2::element_text(size = 9.4),
    axis.text = ggplot2::element_text(size = 8.4)
  )

fig5_design <- "
AABBB
CCDDE
FFFFF
"
figure5 <- p5a + p5b + p5c + p5d + p5e + p5f +
  patchwork::plot_layout(design = fig5_design, guides = "keep") +
  patchwork::plot_annotation(tag_levels = "A", theme = theme_patch_tags)
if (SAVE_COMPOSITE_TIFFS) save_tiff(figure5, "Figure5.tif", width = 7.2, height = 9.0)
save_panel_set(
  list(A = p5a, B = p5b, C = p5c, D = p5d, E = p5e, F = p5f),
  "Figure5",
  width = c(5.8, 6.0, 5.8, 5.8, 6.0, 7.4),
  height = c(4.4, 4.8, 4.6, 4.6, 4.8, 5.6)
)

# =============================== FIGURE 6 ====================================
message("Building Figure 6...")

stage_comp <- stage_competition |>
  dplyr::mutate(stage = factor(stage, levels = stage_levels))
write_source(stage_comp, "Figure6", "A", "Cross-cohort direction consistency of frozen stage signatures")
p6a <- ggplot2::ggplot(stage_comp, ggplot2::aes(stage, proportion_expected_direction_in_both, fill = stage)) +
  ggplot2::geom_col(width = 0.68, colour = "black", linewidth = 0.45) +
  ggplot2::geom_text(ggplot2::aes(label = paste0(n_expected_direction_in_both, "/", n_measured_in_both_core_cohorts)), vjust = -0.35, size = 2.8, family = BASE_FONT) +
  ggplot2::scale_fill_manual(values = stage_colors, guide = "none") +
  ggplot2::scale_x_discrete(labels = stage_labels) +
  ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1.08), expand = c(0,0)) +
  ggplot2::labs(title = "Frozen stage genes show selective directional consistency", x = NULL, y = "Expected direction in both core cohorts") +
  theme_pub() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

stage_gsea_plot <- stage_gsea |>
  dplyr::mutate(
    stage = dplyr::case_when(
      stringr::str_detect(ID, "SKIN") ~ "Skin",
      stringr::str_detect(ID, "WOUND1") ~ "Wound1",
      stringr::str_detect(ID, "WOUND7") ~ "Wound7",
      stringr::str_detect(ID, "WOUND30") ~ "Wound30",
      TRUE ~ Description
    ),
    stage = factor(stage, levels = stage_levels),
    logFDR = -log10(pmax(p.adjust, .Machine$double.xmin))
  )
write_source(stage_gsea_plot, "Figure6", "B", "GSEA of the frozen stage signatures")
p6b <- ggplot2::ggplot(stage_gsea_plot, ggplot2::aes(NES, stage, size = logFDR, colour = stage)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "#888888") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = stage), linewidth = 0.7, colour = "#A5A5A5") +
  ggplot2::geom_point() +
  ggplot2::scale_colour_manual(values = stage_colors, guide = "none") +
  ggplot2::scale_y_discrete(labels = stage_labels) +
  ggplot2::scale_size_continuous(range = c(3, 7), name = "-log10 FDR") +
  ggplot2::labs(title = "Day-30-like enrichment and Skin-like depletion", x = "Normalized enrichment score", y = NULL) +
  theme_pub() + ggplot2::theme(legend.position = "top")

# Gene concordance scatter
late_gene_set <- unique(late_genes$gene)
skin_gene_set <- unique(skin_loss_genes$gene)
concord <- all_gene_effects |>
  dplyr::mutate(
    class = dplyr::case_when(gene %in% late_gene_set ~ "Persistent late-remodeling", gene %in% skin_gene_set ~ "Loss of Skin homeostasis", TRUE ~ "Other genes"),
    class = factor(class, levels = c("Other genes", "Persistent late-remodeling", "Loss of Skin homeostasis"))
  )
write_source(concord, "Figure6", "C", "Cross-cohort all-gene effect concordance")
label_genes <- dplyr::bind_rows(
  concord |>
    dplyr::filter(class == "Persistent late-remodeling") |>
    dplyr::arrange(dplyr::desc(abs(integrated_rank_score))) |>
    dplyr::slice_head(n = 7),
  concord |>
    dplyr::filter(class == "Loss of Skin homeostasis") |>
    dplyr::arrange(dplyr::desc(abs(integrated_rank_score))) |>
    dplyr::slice_head(n = 7)
)
p6c <- ggplot2::ggplot(concord, ggplot2::aes(mean_difference_GSE163973, mean_difference_GSE181316, colour = class)) +
  ggplot2::geom_hline(yintercept = 0, linewidth = 0.35, colour = "#A0A0A0") +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.35, colour = "#A0A0A0") +
  ggplot2::geom_point(data = dplyr::filter(concord, class == "Other genes"), size = 0.65, alpha = 0.22) +
  ggplot2::geom_point(data = dplyr::filter(concord, class != "Other genes"), size = 2.1, alpha = 0.9) +
  ggrepel::geom_text_repel(data = label_genes, ggplot2::aes(label = gene), size = 2.35, family = BASE_FONT, segment.size = 0.3, max.overlaps = Inf) +
  ggplot2::scale_colour_manual(values = c("Other genes" = "#BDBDBD", "Persistent late-remodeling" = "#B2182B", "Loss of Skin homeostasis" = "#2166AC")) +
  ggplot2::labs(title = "Gene-level effects are concordant across core cohorts", x = "GSE163973 mean difference", y = "GSE181316 mean difference", colour = NULL) +
  theme_pub() + ggplot2::theme(legend.position = "top")

# Top consensus gene dumbbell
late_top <- late_genes |>
  dplyr::arrange(dplyr::desc(minimum_oriented_core_effect)) |>
  dplyr::slice_head(n = 12) |>
  dplyr::mutate(class = "Persistent late-remodeling")

skin_top <- skin_loss_genes |>
  dplyr::arrange(dplyr::desc(minimum_oriented_core_effect)) |>
  dplyr::slice_head(n = 10) |>
  dplyr::mutate(class = "Loss of Skin homeostasis")
top_genes <- dplyr::bind_rows(late_top, skin_top) |>
  dplyr::select(gene, class, mean_difference_GSE163973, mean_difference_GSE181316, equal_weight_mean_difference) |>
  dplyr::mutate(gene = factor(gene, levels = rev(gene)))
write_source(top_genes, "Figure6", "D", "Top cross-cohort consensus genes")
p6d <- ggplot2::ggplot(top_genes, ggplot2::aes(y = gene)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "#888888") +
  ggplot2::geom_segment(ggplot2::aes(x = mean_difference_GSE163973, xend = mean_difference_GSE181316, yend = gene), linewidth = 0.75, colour = "#B8B8B8") +
  ggplot2::geom_point(ggplot2::aes(x = mean_difference_GSE163973, colour = "GSE163973"), size = 2.4) +
  ggplot2::geom_point(ggplot2::aes(x = mean_difference_GSE181316, colour = "GSE181316"), size = 2.4) +
  ggplot2::facet_grid(class ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_colour_manual(values = cohort_colors[c("GSE163973", "GSE181316")]) +
  ggplot2::labs(title = "Consensus genes retain their direction in both patient cohorts", x = "Keloid-control mean difference", y = NULL, colour = NULL) +
  theme_pub() + ggplot2::theme(legend.position = "top")

# 3-cohort gene heatmap
third_cols <- names(third_gene_effects)
third_gene_col <- intersect(c("gene", "Gene", "symbol"), third_cols)[1]
third_eff_col <- intersect(c("mean_difference", "mean_difference_GSE220300", "effect", "mean_diff"), third_cols)[1]
if (!is.na(third_gene_col) && !is.na(third_eff_col)) {
  third_small <- third_gene_effects |> dplyr::transmute(gene = .data[[third_gene_col]], GSE220300 = .data[[third_eff_col]])
} else {
  third_small <- gene_contrib |> dplyr::select(gene, GSE220300 = mean_difference_GSE220300)
}
three_gene <- dplyr::bind_rows(late_top |> dplyr::slice_head(n = 10), skin_top |> dplyr::slice_head(n = 8)) |>
  dplyr::select(gene, stage, GSE163973 = mean_difference_GSE163973, GSE181316 = mean_difference_GSE181316) |>
  dplyr::left_join(third_small, by = "gene") |>
  tidyr::pivot_longer(c(GSE163973, GSE181316, GSE220300), names_to = "cohort", values_to = "effect") |>
  dplyr::mutate(gene = factor(gene, levels = rev(unique(gene))))
write_source(three_gene, "Figure6", "E", "Direction of selected consensus genes in two core and one supportive cohort")
p6e <- ggplot2::ggplot(three_gene, ggplot2::aes(cohort, gene, fill = effect)) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.45) +
  ggplot2::geom_text(ggplot2::aes(label = ifelse(is.na(effect), "NA", sprintf("%.1f", effect))), size = 2.15, family = BASE_FONT) +
  ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "Mean\ndifference", na.value = "#E6E6E6") +
  ggplot2::labs(title = "Selected programs retain support across datasets", x = NULL, y = NULL) +
  theme_heatmap() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))

# Hallmark top GSEA
hall_top <- hallmark_gsea |>
  dplyr::filter(
    !is.na(NES),
    !is.na(p.adjust),
    p.adjust < 0.05
  ) |>
  dplyr::arrange(p.adjust, dplyr::desc(abs(NES))) |>
  dplyr::slice_head(n = 12) |>
  dplyr::mutate(
    pathway = stringr::str_wrap(clean_pathway(Description), width = 30),
    pathway = forcats::fct_reorder(pathway, NES),
    logFDR = -log10(pmax(p.adjust, .Machine$double.xmin))
  )
write_source(hall_top, "Figure6", "F", "Top Hallmark GSEA pathways")
p6f <- ggplot2::ggplot(hall_top, ggplot2::aes(NES, pathway, size = logFDR, colour = NES)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "#888888") +
  ggplot2::geom_point() +
  ggplot2::scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, name = "NES") +
  ggplot2::scale_size_continuous(range = c(2.6, 6), name = "-log10 FDR") +
  ggplot2::labs(title = "Hallmark programs associated with persistent remodeling", x = "Normalized enrichment score", y = NULL) +
  theme_pub() + ggplot2::theme(legend.position = "top")

# GO and Reactome top GSEA
select_top_pathways <- function(df, resource, n = 9) {
  df |>
    dplyr::filter(!is.na(NES), !is.na(p.adjust), p.adjust < 0.05) |>
    dplyr::arrange(p.adjust, dplyr::desc(abs(NES))) |>
    dplyr::slice_head(n = n) |>
    dplyr::mutate(
      resource = resource,
      pathway = stringr::str_wrap(
        clean_pathway(Description),
        width = 34
      ),
      logFDR = -log10(
        pmax(p.adjust, .Machine$double.xmin)
      )
    )
}
path_top <- dplyr::bind_rows(
  select_top_pathways(gobp_gsea, "GO biological process", 9),
  select_top_pathways(reactome_gsea, "Reactome", 9)
) |>
  dplyr::group_by(resource) |>
  dplyr::mutate(pathway = forcats::fct_reorder(pathway, NES)) |>
  dplyr::ungroup()
write_source(path_top, "Figure6", "G", "Top GO BP and Reactome GSEA pathways")
p6g <- ggplot2::ggplot(path_top, ggplot2::aes(NES, pathway, size = logFDR, colour = NES)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "#888888") +
  ggplot2::geom_point() +
  ggplot2::facet_grid(resource ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_colour_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, name = "NES") +
  ggplot2::scale_size_continuous(range = c(2.4, 5.8), name = "-log10 FDR") +
  ggplot2::labs(title = "ECM, collagen and matrix-remodeling pathways dominate", x = "Normalized enrichment score", y = NULL) +
  theme_pub() + ggplot2::theme(legend.position = "top", axis.text.y = ggplot2::element_text(size = 7.2))

panel_eff <- panel_effects |>
  dplyr::mutate(
    panel_label = stringr::str_replace_all(panel, "_", " ") |>
      stringr::str_to_sentence() |>
      stringr::str_wrap(width = 30),
    panel_label = forcats::fct_reorder(panel_label, integrated_effect)
  )
write_source(panel_eff, "Figure6", "H", "Predefined fibroblast state-panel effects")
p6h <- ggplot2::ggplot(panel_eff, ggplot2::aes(integrated_effect, panel_label)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.45, colour = "#888888") +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = integrated_effect, yend = panel_label), linewidth = 0.65, colour = "#B0B0B0") +
  ggplot2::geom_point(ggplot2::aes(size = -log10(exact_stratified_two_sided_p), fill = integrated_effect > 0), shape = 21, stroke = 0.5, colour = "black") +
  ggplot2::scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#2166AC"), guide = "none") +
  ggplot2::scale_size_continuous(range = c(2.8, 5.5), name = "-log10 exact P") +
  ggplot2::labs(title = "Predefined fibroblast states provide biological context", x = "Integrated patient-level effect", y = NULL) +
  theme_pub() + ggplot2::theme(legend.position = "top")

fig6_design <- "
AABBB
CCCCC
DDEEE
FFFFF
GGGGH
"
figure6 <- p6a + p6b + p6c + p6d + p6e + p6f + p6g + p6h +
  patchwork::plot_layout(design = fig6_design, guides = "keep") +
  patchwork::plot_annotation(tag_levels = "A", theme = theme_patch_tags)
if (SAVE_COMPOSITE_TIFFS) save_tiff(figure6, "Figure6.tif", width = 7.2, height = 10.8)
save_panel_set(
  list(A = p6a, B = p6b, C = p6c, D = p6d, E = p6e, F = p6f, G = p6g, H = p6h),
  "Figure6",
  width = c(6.0, 6.2, 6.2, 6.0, 6.0, 6.2, 6.2, 5.8),
  height = c(4.8, 5.0, 5.0, 4.8, 4.8, 5.0, 5.0, 4.6)
)

# ========================= SUPPLEMENTARY FIGURE S1 ===========================
message("Building Supplementary Figure S1...")

s1a_data <- ref_qc |> dplyr::mutate(condition = factor(condition, levels = stage_levels))
write_source(s1a_data, "Supplementary_Figure_S1", "A", "GSE241132 sample-level QC summaries")
s1a <- ggplot2::ggplot(s1a_data, ggplot2::aes(condition, qc_retention_rate, colour = donor, group = donor)) +
  ggplot2::geom_line(linewidth = 0.55) + ggplot2::geom_point(size = 2.7) +
  ggplot2::scale_x_discrete(labels = stage_labels) + ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0.8, 1.01)) +
  ggplot2::labs(title = "Reference QC retention", x = NULL, y = "Retention", colour = "Donor") + theme_pub()

s1b_data <- ref_cell_counts |> dplyr::mutate(condition = factor(condition, levels = stage_levels))
write_source(s1b_data, "Supplementary_Figure_S1", "B", "Reference cell-type counts after QC")
s1b <- ggplot2::ggplot(s1b_data, ggplot2::aes(condition, main_cell_type, fill = log10(n_cells + 1))) +
  ggplot2::geom_tile(colour = "white", linewidth = 0.3) +
  ggplot2::facet_wrap(~donor, nrow = 1) + ggplot2::scale_x_discrete(labels = stage_labels) +
  ggplot2::scale_fill_gradient(low = "white", high = "#2166AC", name = "log10 cells") +
  ggplot2::labs(title = "Cell-type coverage across donors and stages", x = NULL, y = NULL) + theme_heatmap()

s1c_data <- s6a_qc |> dplyr::select(sample_id, group, n_called_union, n_QC_keep, n_fibroblast_high, n_fibroblast_broad) |>
  tidyr::pivot_longer(-c(sample_id, group), names_to = "step", values_to = "n")
write_source(s1c_data, "Supplementary_Figure_S1", "C", "GSE181316 cell-calling and fibroblast-yield audit")
s1c <- ggplot2::ggplot(s1c_data, ggplot2::aes(sample_id, n, fill = step)) +
  ggplot2::geom_col(position = "dodge", width = 0.78) + ggplot2::scale_y_log10(labels = scales::comma) +
  ggplot2::labs(title = "Cell calling and fibroblast extraction in GSE181316", x = NULL, y = "Count (log scale)", fill = NULL) +
  theme_pub() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "top")

qc_combined <- dplyr::bind_rows(
  s6a_qc |> dplyr::transmute(cohort = "GSE181316", sample = sample_id, group, retention = QC_retention_from_called, nFeature = nFeature_median, pctMT = median_percent_mt),
  s7_qc |> dplyr::transmute(cohort = "GSE220300", sample = sample_id, group = tissue_state, retention = QC_retention_rate, nFeature = median_nFeature, pctMT = median_percent_mt),
  s9_qc |> dplyr::transmute(cohort = "GSE265972", sample = sample_id, group, retention = audit_QC_pass_fraction, nFeature = median_nFeature, pctMT = median_percent_mt)
)
write_source(qc_combined, "Supplementary_Figure_S1", "D", "Cross-cohort sample-level QC summary")
s1d <- ggplot2::ggplot(qc_combined, ggplot2::aes(cohort, retention, fill = cohort)) +
  ggplot2::geom_point(shape = 21, size = 2.8, stroke = 0.45, colour = "black", position = ggplot2::position_jitter(width = 0.09, seed = 1)) +
  ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0.5,1.02)) + ggplot2::scale_fill_brewer(palette = "Set2", guide = "none") +
  ggplot2::labs(title = "QC retention across external datasets", x = NULL, y = "Retention") + theme_pub()

s1e <- ggplot2::ggplot(qc_combined, ggplot2::aes(cohort, nFeature, fill = cohort)) +
  ggplot2::geom_point(shape = 21, size = 2.8, stroke = 0.45, colour = "black", position = ggplot2::position_jitter(width = 0.09, seed = 2)) +
  ggplot2::scale_fill_brewer(palette = "Set2", guide = "none") +
  ggplot2::labs(title = "Median detected features", x = NULL, y = "Genes per cell/spot") + theme_pub()

s1f <- ggplot2::ggplot(qc_combined, ggplot2::aes(cohort, pctMT, fill = cohort)) +
  ggplot2::geom_point(shape = 21, size = 2.8, stroke = 0.45, colour = "black", position = ggplot2::position_jitter(width = 0.09, seed = 3)) +
  ggplot2::scale_fill_brewer(palette = "Set2", guide = "none") +
  ggplot2::labs(title = "Median mitochondrial proportion", x = NULL, y = "% mitochondrial RNA") + theme_pub()

supp1 <- s1a + s1b + s1c + s1d + s1e + s1f + patchwork::plot_layout(ncol = 2, heights = c(1,1,1)) +
  patchwork::plot_annotation(tag_levels = "A", theme = theme_patch_tags)
if (SAVE_COMPOSITE_TIFFS) save_tiff(supp1, "Supplementary_Figure_S1.tif", 7.2, 9.0, dir = DIR_SUPP)
save_panel_set(
  list(A = s1a, B = s1b, C = s1c, D = s1d, E = s1e, F = s1f),
  "Supplementary_Figure_S1",
  width = c(6.2, 5.8, 5.8, 5.8, 5.6, 5.6),
  height = c(4.8, 4.6, 4.6, 4.6, 4.4, 4.4)
)

# ========================= SUPPLEMENTARY FIGURE S2 ===========================
message("Building Supplementary Figure S2...")
all_metrics_long <- lodo_overall |>
  dplyr::select(cell_type, exact_accuracy, adjacent_accuracy, mean_ordinal_error, pooled_spearman) |>
  tidyr::pivot_longer(-cell_type, names_to = "metric", values_to = "value") |>
  dplyr::mutate(metric = dplyr::recode(metric,
    exact_accuracy = "Exact accuracy", adjacent_accuracy = "Adjacent accuracy",
    mean_ordinal_error = "Mean ordinal error", pooled_spearman = "Pooled Spearman rho"))
write_source(all_metrics_long, "Supplementary_Figure_S2", "A", "LODO performance across major cell types")
s2a <- ggplot2::ggplot(all_metrics_long, ggplot2::aes(cell_type, value, fill = cell_type)) +
  ggplot2::geom_col(width = 0.72) + ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 2) +
  ggplot2::scale_fill_brewer(palette = "Set2", guide = "none") +
  ggplot2::labs(title = "LODO performance across major cell types", x = NULL, y = NULL) + theme_pub() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))

stable_all <- stable_counts |> dplyr::mutate(stage = factor(stage, levels = stage_levels))
write_source(stable_all, "Supplementary_Figure_S2", "B", "Stable signature counts for all modeled cell types")
s2b <- ggplot2::ggplot(stable_all, ggplot2::aes(stage, cell_type, fill = n_stable_genes)) +
  ggplot2::geom_tile(colour = "white") + ggplot2::geom_text(ggplot2::aes(label = n_stable_genes), size = 2.5, family = BASE_FONT) +
  ggplot2::scale_x_discrete(labels = stage_labels) + ggplot2::scale_fill_gradient(low = "white", high = "#2166AC", name = "Stable genes") +
  ggplot2::labs(title = "Stable gene counts across cell types", x = NULL, y = NULL) + theme_heatmap()

conf_all <- confusion |> dplyr::mutate(true_stage = factor(true_stage, levels = stage_levels), predicted_stage = factor(predicted_stage, levels = stage_levels))
write_source(conf_all, "Supplementary_Figure_S2", "C", "LODO confusion matrices across major cell types")
s2c <- ggplot2::ggplot(conf_all, ggplot2::aes(predicted_stage, true_stage, fill = Freq)) +
  ggplot2::geom_tile(colour = "white") + ggplot2::geom_text(ggplot2::aes(label = Freq), size = 2.3, family = BASE_FONT) +
  ggplot2::facet_wrap(~cell_type, ncol = 3) + ggplot2::scale_x_discrete(labels = stage_labels) + ggplot2::scale_y_discrete(labels = stage_labels) +
  ggplot2::scale_fill_gradient(low = "white", high = "#2166AC") + ggplot2::labs(title = "Cell-type-specific confusion matrices", x = "Predicted", y = "True") + theme_pub() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none")

sig_stab_summary <- sig_stability |> dplyr::group_by(cell_type, stage) |>
  dplyr::summarise(n_any = dplyr::n(), n_all3 = sum(stable_all_3_folds), proportion_all3 = n_all3/n_any, .groups = "drop") |>
  dplyr::mutate(stage = factor(stage, levels = stage_levels))
write_source(sig_stab_summary, "Supplementary_Figure_S2", "D", "Cross-fold stability proportion by cell type and stage")
s2d <- ggplot2::ggplot(sig_stab_summary, ggplot2::aes(stage, proportion_all3, colour = cell_type, group = cell_type)) +
  ggplot2::geom_line(linewidth = 0.65) + ggplot2::geom_point(size = 2.6) + ggplot2::scale_x_discrete(labels = stage_labels) +
  ggplot2::scale_y_continuous(labels = scales::percent) + ggplot2::labs(title = "Cross-fold signature stability", x = NULL, y = "Selected in all three folds", colour = NULL) + theme_pub()

supp2 <- s2a + s2b + s2c + s2d + patchwork::plot_layout(ncol = 2, heights = c(1, 1.45)) +
  patchwork::plot_annotation(tag_levels = "A", theme = theme_patch_tags)
if (SAVE_COMPOSITE_TIFFS) save_tiff(supp2, "Supplementary_Figure_S2.tif", 7.2, 7.7, dir = DIR_SUPP)
save_panel_set(
  list(A = s2a, B = s2b, C = s2c, D = s2d),
  "Supplementary_Figure_S2"
)

# ========================= SUPPLEMENTARY FIGURE S3 ===========================
message("Building Supplementary Figure S3...")
if (!is.null(s5_cell_counts)) {
  s3a_data <- s5_cell_counts
  write_source(s3a_data, "Supplementary_Figure_S3", "A", "GSE163973 harmonized cell-type counts")
  s3a <- ggplot2::ggplot(s3a_data, ggplot2::aes(sample_id, main_cell_type, fill = log10(n_cells + 1))) +
    ggplot2::geom_tile(colour = "white") + ggplot2::scale_fill_gradient(low = "white", high = "#2166AC", name = "log10 cells") +
    ggplot2::labs(title = "Cell-type coverage in GSE163973", x = NULL, y = NULL) + theme_heatmap()
} else s3a <- panel_placeholder("GSE163973 cell-type coverage")

s3b_data <- s5_effects |> dplyr::filter(metric %in% c("remodeling_completion", "ordinal_position", "z_Wound30", "z_Skin", "early_state_persistence")) |>
  dplyr::mutate(metric_label = pretty_metric(metric), oriented_g = ifelse(metric == "z_Skin", -hedges_g, hedges_g))
write_source(s3b_data, "Supplementary_Figure_S3", "B", "Discovery effects across available cell types")
s3b <- ggplot2::ggplot(s3b_data, ggplot2::aes(oriented_g, metric_label, colour = cell_type)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "#888888") +
  ggplot2::geom_point(size = 2.5, position = ggplot2::position_dodge(width = 0.45)) +
  ggplot2::labs(title = "Cell-type-specific direction of wound-state effects", x = "Oriented Hedges g", y = NULL, colour = NULL) + theme_pub() + ggplot2::theme(legend.position = "top")

s3c_data <- s5_fib |> dplyr::select(patient_id, group, early_state_persistence, off_trajectory_ratio, mixed_state_entropy) |>
  tidyr::pivot_longer(c(early_state_persistence, off_trajectory_ratio, mixed_state_entropy), names_to = "metric", values_to = "value") |>
  dplyr::mutate(metric_label = pretty_metric(metric), group_label = pretty_group(as.character(group)))
write_source(s3c_data, "Supplementary_Figure_S3", "C", "Discovery negative and mixed candidate endpoints")
s3c <- ggplot2::ggplot(s3c_data, ggplot2::aes(group_label, value, fill = group)) +
  ggplot2::geom_point(shape = 21, size = 2.8, stroke = 0.5, colour = "black", position = ggplot2::position_jitter(width = 0.07, seed = 4)) +
  ggplot2::facet_wrap(~metric_label, scales = "free_y", nrow = 1) + ggplot2::scale_fill_manual(values = group_colors, guide = "none") +
  ggplot2::labs(title = "Competing explanations were not uniformly supported", x = NULL, y = NULL) + theme_pub()

s3d_data <- s5_fib |> dplyr::select(patient_id, group, n_cells, pseudobulk_library_size, pseudobulk_detected_genes) |>
  tidyr::pivot_longer(c(n_cells, pseudobulk_library_size, pseudobulk_detected_genes), names_to = "metric", values_to = "value")
write_source(s3d_data, "Supplementary_Figure_S3", "D", "Discovery pseudobulk sample information")
s3d <- ggplot2::ggplot(s3d_data, ggplot2::aes(patient_id, value, fill = group)) +
  ggplot2::geom_col(width = 0.72) + ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 1) +
  ggplot2::scale_fill_manual(values = group_colors, labels = c("Normal scar", "Keloid")) +
  ggplot2::labs(title = "Pseudobulk support per patient", x = NULL, y = NULL, fill = NULL) + theme_pub() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "top")

supp3 <- s3a + s3b + s3c + s3d + patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A", theme = theme_patch_tags)
if (SAVE_COMPOSITE_TIFFS) save_tiff(supp3, "Supplementary_Figure_S3.tif", 7.2, 7.4, dir = DIR_SUPP)
save_panel_set(
  list(A = s3a, B = s3b, C = s3c, D = s3d),
  "Supplementary_Figure_S3",
  width = c(6.0, 5.8, 5.8, 5.8),
  height = c(4.8, 4.6, 4.6, 4.6)
)

# ========================= SUPPLEMENTARY FIGURE S4 ===========================
message("Building Supplementary Figure S4...")
s4a_data <- s6_effects |> dplyr::filter(metric %in% c("late_remodeling_state", "ordinal_wound_state_position", "z_Wound30", "z_Skin", "early_state_persistence", "off_trajectory_ratio")) |>
  dplyr::mutate(oriented_g = ifelse(metric == "z_Skin", -hedges_g, hedges_g), metric_label = pretty_metric(metric))
write_source(s4a_data, "Supplementary_Figure_S4", "A", "Validation sensitivity across high-specificity and broad fibroblast definitions")
s4a <- ggplot2::ggplot(s4a_data, ggplot2::aes(oriented_g, metric_label, colour = analysis_set)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "#888888") +
  ggplot2::geom_line(ggplot2::aes(group = metric_label), colour = "#B8B8B8", linewidth = 0.5) + ggplot2::geom_point(size = 2.8) +
  ggplot2::labs(title = "Sensitivity to fibroblast-definition stringency", x = "Oriented Hedges g", y = NULL, colour = NULL) + theme_pub()

s4b_data <- s6_scores |> dplyr::select(patient_id, group, analysis_set, late_remodeling_state, ordinal_wound_state_position) |>
  tidyr::pivot_longer(c(late_remodeling_state, ordinal_wound_state_position), names_to = "metric", values_to = "value") |>
  dplyr::mutate(metric_label = pretty_metric(metric), group = dplyr::recode(group, healthy_skin = "normal_scar", normal_skin = "normal_scar", .default = group), group_label = pretty_group(group))
write_source(s4b_data, "Supplementary_Figure_S4", "B", "Validation patient values under both fibroblast definitions")
s4b <- ggplot2::ggplot(s4b_data, ggplot2::aes(group_label, value, fill = group)) +
  ggplot2::geom_point(shape = 21, size = 2.6, stroke = 0.45, colour = "black", position = ggplot2::position_jitter(width = 0.06, seed = 5)) +
  ggplot2::facet_grid(metric_label ~ analysis_set, scales = "free_y") + ggplot2::scale_fill_manual(values = group_colors, guide = "none") +
  ggplot2::labs(title = "Patient-level sensitivity analysis", x = NULL, y = NULL) + theme_pub()

if (!is.null(s6_gene_coverage)) {
  write_source(s6_gene_coverage, "Supplementary_Figure_S4", "C", "Validation signature-gene coverage")
  s4c <- ggplot2::ggplot(dplyr::filter(s6_gene_coverage, analysis_set == "high_specificity"), ggplot2::aes(stage, coverage_fraction, fill = stage)) +
    ggplot2::geom_col(width = 0.68, colour = "black", linewidth = 0.4) + ggplot2::scale_fill_manual(values = stage_colors, guide = "none") +
    ggplot2::scale_x_discrete(labels = stage_labels) + ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0,1)) +
    ggplot2::labs(title = "Frozen signature coverage in GSE181316", x = NULL, y = "Genes measured") + theme_pub()
} else s4c <- panel_placeholder("Frozen signature coverage")

write_source(k3_lr, "Supplementary_Figure_S4", "D", "Full K3 left-right source data")
s4d <- p3h + ggplot2::labs(title = "K3 left-right agreement across all reported metrics")

supp4 <- s4a + s4b + s4c + s4d + patchwork::plot_layout(ncol = 2, widths = c(1,1.35)) +
  patchwork::plot_annotation(tag_levels = "A", theme = theme_patch_tags)
if (SAVE_COMPOSITE_TIFFS) save_tiff(supp4, "Supplementary_Figure_S4.tif", 7.2, 7.4, dir = DIR_SUPP)
save_panel_set(
  list(A = s4a, B = s4b, C = s4c, D = s4d),
  "Supplementary_Figure_S4",
  width = c(6.0, 5.8, 5.8, 5.8),
  height = c(4.8, 4.6, 4.6, 4.6)
)

# ========================= SUPPLEMENTARY FIGURE S5 ===========================
message("Building Supplementary Figure S5...")
s5a_data <- s7_center_periph |> dplyr::filter(analysis_set == "high_specificity", metric %in% c("late_remodeling_state", "ordinal_wound_state_position")) |>
  dplyr::mutate(metric_label = pretty_metric(metric))
write_source(s5a_data, "Supplementary_Figure_S5", "A", "GSE220300 paired center-periphery differences")
s5a <- ggplot2::ggplot(s5a_data, ggplot2::aes(condition_low, low_value, group = pair_id, colour = activity_group)) +
  ggplot2::geom_segment(ggplot2::aes(xend = condition_high, yend = high_value), linewidth = 0.65, arrow = grid::arrow(length = grid::unit(0.07, "inches"))) +
  ggplot2::geom_point(size = 2.7) + ggplot2::facet_wrap(~metric_label, scales = "free_y") +
  ggplot2::labs(title = "Center-periphery gradients were heterogeneous", x = NULL, y = NULL, colour = "Activity") + theme_pub()

s5b_data <- s7_inactive_scar |> dplyr::filter(analysis_set == "high_specificity", metric %in% c("late_remodeling_state", "ordinal_wound_state_position", "z_Wound30", "z_Skin")) |>
  dplyr::mutate(metric_label = pretty_metric(metric))
write_source(s5b_data, "Supplementary_Figure_S5", "B", "Inactive lesion versus paired mature scar")
s5b <- ggplot2::ggplot(s5b_data, ggplot2::aes(condition_low, low_value, group = pair_id)) +
  ggplot2::geom_segment(ggplot2::aes(xend = condition_high, yend = high_value), linewidth = 0.65, colour = "#666666", arrow = grid::arrow(length = grid::unit(0.07, "inches"))) +
  ggplot2::geom_point(size = 2.7, colour = "#B2182B") + ggplot2::facet_wrap(~metric_label, scales = "free_y") +
  ggplot2::labs(title = "Inactive lesions versus paired mature scars", x = NULL, y = NULL) + theme_pub()

s5c_data <- s7_activity |> dplyr::filter(analysis_set == "high_specificity") |> dplyr::mutate(metric_label = pretty_metric(metric))
write_source(s5c_data, "Supplementary_Figure_S5", "C", "Active versus inactive descriptive effects")
s5c <- ggplot2::ggplot(s5c_data, ggplot2::aes(mean_difference, metric_label)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "#888888") + ggplot2::geom_point(size = 2.8, colour = "#B2182B") +
  ggplot2::labs(title = "Activity status did not produce a stable gradient", x = "Active-inactive mean difference", y = NULL) + theme_pub()

s5d_data <- s7_lesion_scar |> dplyr::filter(analysis_set == "high_specificity") |> dplyr::mutate(metric_label = pretty_metric(metric))
write_source(s5d_data, "Supplementary_Figure_S5", "D", "Lesion versus scar-reference descriptive effects")
s5d <- ggplot2::ggplot(s5d_data, ggplot2::aes(mean_difference, metric_label)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "#888888") + ggplot2::geom_point(ggplot2::aes(size = abs(cliffs_delta)), colour = "#B2182B") +
  ggplot2::scale_size_continuous(range = c(2.5,5), name = "|Cliff's delta|") +
  ggplot2::labs(title = "Overall lesion-versus-scar support", x = "Mean difference", y = NULL) + theme_pub()

supp5 <- s5a + s5b + s5c + s5d + patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A", theme = theme_patch_tags)
if (SAVE_COMPOSITE_TIFFS) save_tiff(supp5, "Supplementary_Figure_S5.tif", 7.2, 7.4, dir = DIR_SUPP)
save_panel_set(
  list(A = s5a, B = s5b, C = s5c, D = s5d),
  "Supplementary_Figure_S5",
  width = c(6.0, 5.8, 5.8, 5.8),
  height = c(4.8, 4.6, 4.6, 4.6)
)

# ========================= SUPPLEMENTARY FIGURE S6 ===========================
message("Building Supplementary Figure S6...")
all_spots_high <- s8_spot_preview |>
  dplyr::filter(high_fibroblast_enriched) |>
  dplyr::mutate(condition = factor(condition, levels = stage_levels))

write_source(
  all_spots_high,
  "Supplementary_Figure_S6",
  "A_all_spatial_spots",
  "All high-specificity fibroblast-enriched spatial spots used for donor-split spatial maps"
)

make_spatial_donor_panel <- function(
  donor_id,
  title_prefix = "Normal-wound spatial sections"
) {
  d <- all_spots_high |>
    dplyr::filter(donor == donor_id) |>
    dplyr::mutate(condition = factor(condition, levels = stage_levels))

  if (!nrow(d)) {
    return(
      panel_placeholder(
        paste(title_prefix, donor_id),
        "No high-specificity fibroblast-enriched spot"
      )
    )
  }

  ggplot2::ggplot(
    d,
    ggplot2::aes(
      pxl_col_in_fullres,
      -pxl_row_in_fullres,
      colour = spot_late_remodeling
    )
  ) +
    ggplot2::geom_point(
      size = 0.75,
      alpha = 0.90
    ) +
    ggplot2::facet_wrap(
      ~ condition,
      nrow = 1,
      labeller = ggplot2::as_labeller(stage_labels)
    ) +
    ggplot2::scale_colour_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      name = "Late-remodeling\nspot score",
      guide = ggplot2::guide_colourbar(
        title.position = "top",
        title.hjust = 0.5,
        barheight = grid::unit(36, "pt"),
        barwidth = grid::unit(7, "pt")
      )
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(
      title = paste0(title_prefix, ": ", donor_id),
      x = NULL,
      y = NULL
    ) +
    theme_pub(base_size = 8.4, title_size = 9.0) +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      legend.position = "right",
      strip.text = ggplot2::element_text(size = 8.0, face = "bold"),
      plot.margin = ggplot2::margin(8, 12, 8, 8)
    )
}

spatial_donors <- sort(
  unique(
    as.character(
      all_spots_high$donor
    )
  )
)

s6a_donor_plots <- stats::setNames(
  lapply(
    spatial_donors,
    make_spatial_donor_panel
  ),
  paste0(
    "A_Donor_",
    safe_file_id(spatial_donors)
  )
)

# Save each donor-specific spatial map as a separate TIFF. This avoids the
# unreadable 4-by-4 spatial panel that compressed all 16 sections into one plot.
if (
  SAVE_PANEL_TIFFS &&
    length(
      s6a_donor_plots
    )
) {
  for (
    nm in names(
      s6a_donor_plots
    )
  ) {
    save_tiff(
      s6a_donor_plots[[nm]],
      paste0(
        "Supplementary_Figure_S6_",
        nm,
        ".tif"
      ),
      width = 8.2,
      height = 3.8,
      dir = DIR_PANELS
    )
  }
}

# The default S6A panel is the first donor-specific map. All donors are exported
# separately with explicit donor-specific file names above.
s6a <- if (
  length(
    s6a_donor_plots
  )
) {
  s6a_donor_plots[[1]]
} else {
  panel_placeholder(
    "Normal-wound spatial sections",
    "No high-specificity fibroblast-enriched spot"
  )
}

s6b_data <- s8_sections |> dplyr::mutate(condition = factor(condition, levels = stage_levels))
write_source(s6b_data, "Supplementary_Figure_S6", "B", "Spatial validation under high-specificity and broad spot definitions")
s6b <- ggplot2::ggplot(s6b_data, ggplot2::aes(stage_code, ordinal_wound_state_position, group = donor, colour = donor)) +
  ggplot2::geom_line(linewidth = 0.55) + ggplot2::geom_point(size = 2.3) + ggplot2::facet_wrap(~analysis_set, nrow = 1) +
  ggplot2::scale_x_continuous(breaks = 0:3, labels = stage_labels) +
  ggplot2::labs(title = "Spatial ordinal validation under alternative spot definitions", x = "Actual stage", y = "Predicted position", colour = "Donor") + theme_pub()

s6c_data <- s8_overall
write_source(s6c_data, "Supplementary_Figure_S6", "C", "Overall spatial validation metrics")
s6c <- s8_overall |> dplyr::select(analysis_set, exact_accuracy, adjacent_accuracy, mean_ordinal_error, pooled_spearman) |>
  tidyr::pivot_longer(-analysis_set, names_to = "metric", values_to = "value") |>
  ggplot2::ggplot(ggplot2::aes(analysis_set, value, fill = analysis_set)) + ggplot2::geom_col(width = 0.65) +
  ggplot2::facet_wrap(~metric, scales = "free_y", ncol = 2) + ggplot2::scale_fill_brewer(palette = "Set2", guide = "none") +
  ggplot2::labs(title = "Overall cross-modal performance", x = NULL, y = NULL) + theme_pub() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

supp6 <- s6a / (s6b + s6c) + patchwork::plot_layout(heights = c(2.2,1)) +
  patchwork::plot_annotation(tag_levels = "A", theme = theme_patch_tags)
if (SAVE_COMPOSITE_TIFFS) save_tiff(supp6, "Supplementary_Figure_S6.tif", 7.2, 9.5, dir = DIR_SUPP)
save_panel_set(
  list(A = s6a, B = s6b, C = s6c),
  "Supplementary_Figure_S6"
)

# ========================= SUPPLEMENTARY FIGURE S7 ===========================
message("Building Supplementary Figure S7...")
s7a_data <- s10_scores |> dplyr::filter(stringr::str_detect(analysis_set, "high_specificity"))
write_source(s7a_data, "Supplementary_Figure_S7", "A", "GSE181297 high-specificity cross-modal wound-state scores")
s7a <- heatmap_stage_scores(s7a_data, "sample_key", "group", "GSE181297 high-specificity stage scores")

s7b_data <- s10_direction |> dplyr::select(analysis_set, modality, late_remodeling_mean_difference, ordinal_mean_difference) |>
  tidyr::pivot_longer(c(late_remodeling_mean_difference, ordinal_mean_difference), names_to = "metric", values_to = "effect") |>
  dplyr::mutate(metric = dplyr::recode(metric, late_remodeling_mean_difference = "Late-remodeling state", ordinal_mean_difference = "Ordinal position"))
write_source(s7b_data, "Supplementary_Figure_S7", "B", "GSE181297 modality-specific direction summary")
s7b <- ggplot2::ggplot(s7b_data, ggplot2::aes(effect, analysis_set, colour = modality)) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, colour = "#888888") + ggplot2::geom_point(size = 3) +
  ggplot2::facet_wrap(~metric, scales = "free_x") + ggplot2::labs(title = "Direction depends on modality and cell-definition stringency", x = "Keloid-control mean difference", y = NULL, colour = NULL) + theme_pub()

write_source(s10_spot_summary, "Supplementary_Figure_S7", "C", "GSE181297 Visium spot-distribution summary")
s7c <- ggplot2::ggplot(s10_spot_summary, ggplot2::aes(group, high_median_spot_late_remodeling, fill = group)) +
  ggplot2::geom_point(shape = 21, size = 3.3, stroke = 0.5, colour = "black", position = ggplot2::position_jitter(width = 0.06, seed = 6)) +
  ggrepel::geom_text_repel(ggplot2::aes(label = sample_id), size = 2.6, family = BASE_FONT, segment.size = 0.3) +
  ggplot2::scale_fill_manual(values = group_colors, guide = "none") +
  ggplot2::labs(title = "Visium spot-level late-remodeling summaries", x = NULL, y = "Median spot score") + theme_pub()

write_source(s10_qc, "Supplementary_Figure_S7", "D", "GSE181297 Visium QC and fibroblast-enriched spot yield")
s7d <- s10_qc |> dplyr::select(sample_id, group, n_QC_spots, n_high_fibroblast_enriched, n_broad_fibroblast_enriched) |>
  tidyr::pivot_longer(-c(sample_id, group), names_to = "metric", values_to = "n") |>
  ggplot2::ggplot(ggplot2::aes(sample_id, n, fill = metric)) + ggplot2::geom_col(position = "dodge") +
  ggplot2::labs(title = "Visium spot-yield audit", x = NULL, y = "Spots", fill = NULL) + theme_pub() + ggplot2::theme(legend.position = "top")

supp7 <- s7a + s7b + s7c + s7d + patchwork::plot_layout(ncol = 2) +
  patchwork::plot_annotation(tag_levels = "A", theme = theme_patch_tags)
if (SAVE_COMPOSITE_TIFFS) save_tiff(supp7, "Supplementary_Figure_S7.tif", 7.2, 7.5, dir = DIR_SUPP)
save_panel_set(
  list(A = s7a, B = s7b, C = s7c, D = s7d),
  "Supplementary_Figure_S7"
)

# ========================= SUPPLEMENTARY FIGURE S8 ===========================
message("Building Supplementary Figure S8...")

prepare_ora_resource <- function(data, resource, module_name) {
  data |>
    dplyr::filter(
      !is.na(p.adjust),
      p.adjust < 0.05
    ) |>
    dplyr::mutate(
      resource = resource,
      module = module_name,
      pathway = stringr::str_wrap(
        clean_pathway(Description),
        width = 34
      ),
      score = -log10(
        pmax(p.adjust, .Machine$double.xmin)
      )
    )
}

late_ora_significant <- dplyr::bind_rows(
  prepare_ora_resource(
    late_hallmark_ora,
    "Hallmark",
    "Persistent late remodeling"
  ),
  prepare_ora_resource(
    late_gobp_ora,
    "GO BP",
    "Persistent late remodeling"
  ),
  prepare_ora_resource(
    late_reactome_ora,
    "Reactome",
    "Persistent late remodeling"
  )
) |>
  dplyr::arrange(p.adjust, dplyr::desc(FoldEnrichment)) |>
  dplyr::slice_head(n = 18) |>
  dplyr::mutate(
    pathway = forcats::fct_reorder(pathway, score)
  )

write_source(
  late_ora_significant,
  "Supplementary_Figure_S8",
  "A",
  "FDR-significant ORA terms for the persistent late-remodeling module"
)

if (nrow(late_ora_significant)) {
  s8a <- ggplot2::ggplot(
    late_ora_significant,
    ggplot2::aes(
      score,
      pathway,
      size = Count,
      colour = FoldEnrichment
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(
      resource ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_colour_gradient(
      low = "#FDB863",
      high = "#B2182B"
    ) +
    ggplot2::labs(
      title = "Persistent late-remodeling module: FDR-significant ORA",
      x = "-log10 FDR",
      y = NULL,
      colour = "Fold enrichment",
      size = "Genes"
    ) +
    theme_pub(base_size = 8)
} else {
  s8a <- panel_placeholder(
    "Persistent late-remodeling ORA",
    "No FDR-significant term"
  )
}

skin_ora_significant <- dplyr::bind_rows(
  prepare_ora_resource(
    skin_hallmark_ora,
    "Hallmark",
    "Loss of Skin homeostasis"
  ),
  prepare_ora_resource(
    skin_gobp_ora,
    "GO BP",
    "Loss of Skin homeostasis"
  ),
  prepare_ora_resource(
    skin_reactome_ora,
    "Reactome",
    "Loss of Skin homeostasis"
  )
) |>
  dplyr::arrange(p.adjust, dplyr::desc(FoldEnrichment)) |>
  dplyr::slice_head(n = 18) |>
  dplyr::mutate(
    pathway = forcats::fct_reorder(pathway, score)
  )

write_source(
  skin_ora_significant,
  "Supplementary_Figure_S8",
  "B",
  "FDR-significant ORA terms for the Skin-homeostasis-loss module"
)

if (nrow(skin_ora_significant)) {
  s8b <- ggplot2::ggplot(
    skin_ora_significant,
    ggplot2::aes(
      score,
      pathway,
      size = Count,
      colour = FoldEnrichment
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::facet_grid(
      resource ~ .,
      scales = "free_y",
      space = "free_y"
    ) +
    ggplot2::scale_colour_gradient(
      low = "#92C5DE",
      high = "#2166AC"
    ) +
    ggplot2::labs(
      title = "Skin-homeostasis-loss module: FDR-significant ORA",
      x = "-log10 FDR",
      y = NULL,
      colour = "Fold enrichment",
      size = "Genes"
    ) +
    theme_pub(base_size = 8)
} else {
  s8b <- panel_placeholder(
    "Skin-homeostasis-loss ORA",
    "No FDR-significant term"
  )
}

write_source(
  panel_scores,
  "Supplementary_Figure_S8",
  "C",
  "Patient-level predefined fibroblast-panel scores"
)

s8c <- panel_scores |>
  dplyr::mutate(
    panel_label = stringr::str_replace_all(panel, "_", " ") |>
      stringr::str_to_sentence() |>
      stringr::str_wrap(width = 28),
    group_label = pretty_group(group)
  ) |>
  ggplot2::ggplot(
    ggplot2::aes(group_label, module_score, fill = group)
  ) +
  ggplot2::geom_point(
    shape = 21,
    size = 2.3,
    stroke = 0.4,
    colour = "black",
    position = ggplot2::position_jitter(
      width = 0.06,
      seed = 7
    )
  ) +
  ggplot2::facet_grid(
    panel_label ~ cohort,
    scales = "free_y"
  ) +
  ggplot2::scale_fill_manual(
    values = group_colors,
    guide = "none"
  ) +
  ggplot2::labs(
    title = "Patient-level predefined fibroblast-state scores",
    x = NULL,
    y = NULL
  ) +
  theme_pub(base_size = 8)

all_path_top <- dplyr::bind_rows(
  hallmark_gsea |>
    dplyr::mutate(resource = "Hallmark"),
  gobp_gsea |>
    dplyr::mutate(resource = "GO BP"),
  reactome_gsea |>
    dplyr::mutate(resource = "Reactome")
) |>
  dplyr::filter(
    !is.na(NES),
    !is.na(p.adjust),
    p.adjust < 0.05
  ) |>
  dplyr::mutate(
    pathway = clean_pathway_short(Description, width = 36),
    logFDR = -log10(
      pmax(p.adjust, .Machine$double.xmin)
    )
  )

write_source(
  all_path_top,
  "Supplementary_Figure_S8",
  "D_all_FDR_significant_GSEA",
  "All FDR-significant GSEA results before resource-specific plotting"
)

s8d_hallmark_data <- all_path_top |>
  dplyr::filter(resource == "Hallmark") |>
  dplyr::arrange(p.adjust, dplyr::desc(abs(NES))) |>
  dplyr::mutate(pathway_clean = clean_pathway_short(Description, width = 34)) |>
  dplyr::distinct(pathway_clean, .keep_all = TRUE) |>
  dplyr::slice_head(n = 12)

s8d_gobp_data <- all_path_top |>
  dplyr::filter(resource == "GO BP") |>
  dplyr::arrange(p.adjust, dplyr::desc(abs(NES))) |>
  dplyr::mutate(pathway_clean = clean_pathway_short(Description, width = 34)) |>
  dplyr::distinct(pathway_clean, .keep_all = TRUE) |>
  dplyr::slice_head(n = 12)

s8d_reactome_data <- all_path_top |>
  dplyr::filter(resource == "Reactome") |>
  dplyr::arrange(p.adjust, dplyr::desc(abs(NES))) |>
  dplyr::mutate(pathway_clean = clean_pathway_short(Description, width = 34)) |>
  dplyr::distinct(pathway_clean, .keep_all = TRUE) |>
  dplyr::slice_head(n = 12)

write_source(
  s8d_hallmark_data,
  "Supplementary_Figure_S8",
  "D1_Hallmark_GSEA",
  "FDR-significant Hallmark GSEA terms after de-duplication"
)

write_source(
  s8d_gobp_data,
  "Supplementary_Figure_S8",
  "D2_GO_BP_GSEA",
  "FDR-significant GO BP GSEA terms after de-duplication"
)

write_source(
  s8d_reactome_data,
  "Supplementary_Figure_S8",
  "D3_Reactome_GSEA",
  "FDR-significant Reactome GSEA terms after de-duplication"
)

s8d_hallmark <- pathway_resource_plot(
  all_path_top,
  "Hallmark",
  "Hallmark GSEA",
  n_terms = 12,
  wrap_width = 34
)

s8d_gobp <- pathway_resource_plot(
  all_path_top,
  "GO BP",
  "GO biological process GSEA",
  n_terms = 12,
  wrap_width = 34
)

s8d_reactome <- pathway_resource_plot(
  all_path_top,
  "Reactome",
  "Reactome GSEA",
  n_terms = 12,
  wrap_width = 34
)

# Composite Supplementary Figure S8 is disabled by default in the panels-only
# build. For backward compatibility, keep s8d as the Hallmark panel.
s8d <- s8d_hallmark

supp8 <- (s8a + s8b) / (s8c + s8d) +
  patchwork::plot_layout(heights = c(1, 1.5)) +
  patchwork::plot_annotation(tag_levels = "A", theme = theme_patch_tags)

if (SAVE_COMPOSITE_TIFFS) save_tiff(
  supp8,
  "Supplementary_Figure_S8.tif",
  7.2,
  9.2,
  dir = DIR_SUPP
)

save_panel_set(
  list(
    A = s8a,
    B = s8b,
    C = s8c,
    D1_Hallmark = s8d_hallmark,
    D2_GO_BP = s8d_gobp,
    D3_Reactome = s8d_reactome
  ),
  "Supplementary_Figure_S8",
  width = c(5.8, 5.8, 6.2, 7.0, 7.0, 7.0),
  height = c(4.6, 4.6, 5.0, 5.8, 5.8, 5.8)
)

# ================================ TABLES =====================================
message("Exporting manuscript tables...")

# Main Table 1: dataset architecture
main_table1 <- tibble::tribble(
  ~Dataset, ~Modality, ~Role, ~Biological_units, ~Primary_use,
  "GSE241132", "scRNA-seq", "Physiological reference", "3 donors × 4 stages", "Frozen fibroblast wound-state reference and donor-held-out validation",
  "GSE163973", "scRNA-seq", "Discovery", "3 keloids + 3 normal scars", "Discovery of persistent late remodeling",
  "GSE181316", "scRNA-seq", "Independent replication", "3 keloids + 3 normal scars", "Frozen patient-level replication",
  "GSE220300", "scRNA-seq", "External support", "Patient-level lesion/scar units", "Direction support across lesion contexts",
  "GSE241124", "Spatial transcriptomics", "Cross-modal validation", "4 donors × 4 stages", "New-donor and section-level ordinal validation",
  "GSE265972", "scRNA-seq", "Disease contrast", "4 venous ulcers + 5 normal skin", "Distinct wound-resolution failure modes",
  "GSE181297", "scRNA-seq + Visium", "Exploratory supplementary", "Small cross-modal dataset", "Mixed modality-dependent evidence"
)
readr::write_csv(main_table1, file.path(DIR_TABLES, "Table1_Dataset_architecture.csv"))

# Main Table 2: integrated core effects
main_table2 <- core_integrated |>
  dplyr::mutate(metric_display = pretty_metric(metric)) |>
  dplyr::select(
    metric_display, endpoint_role, equal_weight_mean_difference,
    mean_hedges_g, mean_cliffs_delta,
    exact_stratified_two_sided_p, exact_stratified_directional_p,
    all_cohorts_expected_direction, all_cohorts_complete_separation
  )
readr::write_csv(main_table2, file.path(DIR_TABLES, "Table2_Core_effects_and_exact_integration.csv"))

# Supplementary tables: direct, transparent exports
supp_tables <- list(
  "Supplementary_Table_S1_Master_sample_dictionary.csv" = master_dictionary,
  "Supplementary_Table_S2_QC_summary_reference.csv" = ref_qc,
  "Supplementary_Table_S3_Frozen_stage_signatures.csv" = fixed_signatures,
  "Supplementary_Table_S4_Discovery_patient_scores.csv" = s5_fib,
  "Supplementary_Table_S5_Replication_patient_scores.csv" = s6_primary,
  "Supplementary_Table_S6_External_support_scores.csv" = s7_patient,
  "Supplementary_Table_S7_Spatial_section_scores.csv" = s8_sections,
  "Supplementary_Table_S8_Venous_ulcer_scores.csv" = s9_primary,
  "Supplementary_Table_S9_Consensus_genes.csv" = dplyr::bind_rows(
    late_genes |> dplyr::mutate(consensus_class = "Persistent late remodeling"),
    skin_loss_genes |> dplyr::mutate(consensus_class = "Loss of Skin homeostasis")
  ),
  "Supplementary_Table_S10_Pathway_results.csv" = dplyr::bind_rows(
    hallmark_gsea |> dplyr::mutate(resource = "Hallmark"),
    gobp_gsea |> dplyr::mutate(resource = "GO BP"),
    reactome_gsea |> dplyr::mutate(resource = "Reactome")
  )
)
for (nm in names(supp_tables)) {
  x <- supp_tables[[nm]]
  if (!is.null(x)) readr::write_csv(x, file.path(DIR_TABLES, nm), na = "")
}

# ------------------------------- FINAL LOGS ----------------------------------
readr::write_csv(
  source_index,
  file.path(DIR_SOURCE, "Source_Data_Index.csv")
)

figure_blueprint <- tibble::tribble(
  ~figure, ~panel_count, ~width_in, ~height_in, ~scientific_role,
  "Figure1", 7L, 7.2, 9.0, "Donor-validated physiological wound-state reference without a workflow diagram",
  "Figure2", 6L, 7.2, 8.2, "Discovery of predominant late-remodeling persistence and loss of Skin-like homeostasis",
  "Figure3", 8L, 7.2, 9.2, "Frozen independent replication and stratified exact evidence synthesis",
  "Figure4", 8L, 7.2, 9.1, "Third keloid cohort support and cross-modal normal-wound validation",
  "Figure5", 6L, 7.2, 8.0, "Distinct wound-resolution failure modes in keloid and venous ulcer",
  "Figure6", 8L, 7.2, 9.3, "Cross-cohort consensus genes, pathways and fibroblast-state context",
  "Supplementary_Figure_S1", 6L, 7.2, 9.0, "QC, cell yield and biological-unit coverage",
  "Supplementary_Figure_S2", 4L, 7.2, 7.7, "Reference-model performance across cell types",
  "Supplementary_Figure_S3", 4L, 7.2, 7.4, "Discovery coverage and competing hypotheses",
  "Supplementary_Figure_S4", 4L, 7.2, 7.4, "Replication sensitivity and within-patient consistency",
  "Supplementary_Figure_S5", 4L, 7.2, 7.4, "Regional, activity and paired-scar exploratory analyses",
  "Supplementary_Figure_S6", 3L, 7.2, 9.5, "Complete normal-wound spatial validation",
  "Supplementary_Figure_S7", 4L, 7.2, 7.5, "Mixed GSE181297 cross-modal evidence",
  "Supplementary_Figure_S8", 4L, 7.2, 9.2, "Extended FDR-controlled gene-set interpretation"
)

readr::write_csv(
  figure_blueprint,
  file.path(DIR_LOGS, "FIGURE_BLUEPRINT_AND_DIMENSIONS.csv")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(DIR_LOGS, "SESSION_INFO_FIGURE_PIPELINE.txt")
)

# In this final workflow, this script is intentionally run in panels-only mode.
# Composite figures are generated by the second script after all individual panels
# have been written. Therefore, do not require main_tif/ or supplementary_tif/
# composites at this stage.
if (isTRUE(SAVE_COMPOSITE_TIFFS)) {
  expected_main <- file.path(
    DIR_MAIN,
    paste0("Figure", 1:6, ".tif")
  )

  expected_supp <- file.path(
    DIR_SUPP,
    paste0("Supplementary_Figure_S", 1:8, ".tif")
  )

  missing_composites <- c(
    expected_main[!file.exists(expected_main)],
    expected_supp[!file.exists(expected_supp)]
  )

  if (length(missing_composites)) {
    stop(
      "The pipeline ended with missing composite TIFF files:\n",
      paste(missing_composites, collapse = "\n")
    )
  }
}

all_output_files <- list.files(
  OUTPUT_DIR,
  recursive = TRUE,
  full.names = TRUE,
  include.dirs = FALSE
)

checksum_table <- tibble::tibble(
  file = substring(
    normalizePath(all_output_files, winslash = "/", mustWork = TRUE),
    nchar(normalizePath(OUTPUT_DIR, winslash = "/", mustWork = TRUE)) + 2L
  ),
  size_bytes = file.info(all_output_files)$size,
  md5 = unname(tools::md5sum(all_output_files))
)

readr::write_csv(
  checksum_table,
  file.path(DIR_LOGS, "OUTPUT_FILE_CHECKSUMS.csv")
)

completion_text <- c(
  "PUBLICATION FIGURE PIPELINE COMPLETED",
  "Version: FINAL SINGLE PANELS + FIGURE2B DISEASE PROJECTION",
  paste0("Date: ", Sys.time()),
  paste0("Project root: ", PROJECT_ROOT),
  paste0("Output directory: ", OUTPUT_DIR),
  "Individual panel TIFFs exported; Figure2B uses disease-projection emphasis",
  "Supplementary figures: Supplementary_Figure_S1.tif through S8.tif",
  paste0(
    "Individual panel TIFFs saved: ",
    SAVE_PANEL_TIFFS
  ),
  "All panel-level source data are indexed in source_data/Source_Data_Index.csv",
  "All figure image outputs are TIFF only, white background, 600 dpi.",
  "No workflow/flowchart panel is included.",
  "Main figures retain the positive prespecified evidence chain; competing and mixed findings remain supplementary."
)

writeLines(
  completion_text,
  file.path(DIR_LOGS, "FIGURE_PIPELINE_COMPLETED.txt")
)

# Copy the exact script used, when available.
if (nzchar(script_path) && file.exists(script_path)) {
  file.copy(
    script_path,
    file.path(DIR_LOGS, basename(script_path)),
    overwrite = TRUE
  )
}

# Rebuild final output inventory after all logs have been written.
output_inventory <- list.files(
  OUTPUT_DIR,
  recursive = TRUE,
  full.names = FALSE,
  include.dirs = FALSE
)

readr::write_csv(
  tibble::tibble(file = output_inventory),
  file.path(DIR_LOGS, "OUTPUT_FILE_INVENTORY.csv")
)

# Create a review package. Image files inside the package remain TIFF only.
review_zip <- file.path(
  PROJECT_ROOT,
  "瘢痕疙瘩_投稿级图表_TIFF与源数据_FINAL.zip"
)

if (file.exists(review_zip)) {
  file.remove(review_zip)
}

zip_success <- FALSE

if (.Platform$OS.type == "windows") {
  ps_command <- paste0(
    "$ErrorActionPreference='Stop'; ",
    "$items=Get-ChildItem -LiteralPath '",
    gsub("'", "''", normalizePath(OUTPUT_DIR, winslash = "/", mustWork = TRUE), fixed = TRUE),
    "'; ",
    "Compress-Archive -Path $items.FullName -DestinationPath '",
    gsub("'", "''", normalizePath(review_zip, winslash = "/", mustWork = FALSE), fixed = TRUE),
    "' -Force"
  )

  ps_result <- try(
    system2(
      "powershell.exe",
      c(
        "-NoProfile",
        "-ExecutionPolicy",
        "Bypass",
        "-Command",
        ps_command
      ),
      stdout = TRUE,
      stderr = TRUE
    ),
    silent = TRUE
  )

  zip_success <- file.exists(review_zip) &&
    file.info(review_zip)$size > 0
}

if (!zip_success) {
  old_wd <- getwd()
  setwd(PROJECT_ROOT)

  try(
    utils::zip(
      zipfile = review_zip,
      files = basename(OUTPUT_DIR),
      flags = "-r9X"
    ),
    silent = TRUE
  )

  setwd(old_wd)

  zip_success <- file.exists(review_zip) &&
    file.info(review_zip)$size > 0
}

if (!zip_success) {
  warning(
    "All figures were generated, but automatic ZIP creation failed. ",
    "Compress the output folder manually."
  )
}

message(paste(completion_text, collapse = "\n"))

if (zip_success) {
  message("Review package: ", review_zip)
}




# ============================================================================
# PART 2/2: FINAL COMPOSITE FIGURE STITCHING
# ============================================================================

# ============================================================================
# Keloid study - optimized final composite figure stitching script
# Version: FINAL COMPOSITES + F2B + FIG4 RELAYOUT | FIXED HANDOFF | 2026-06-19
# ============================================================================
# PURPOSE
#   Stitch already-generated individual TIFF panels into final composite TIFFs
#   for all main figures (Figure 1-6) and supplementary figures (S1-S8).
#
# DESIGN PRINCIPLE
#   - Do NOT force all figures to the same size.
#   - Use the best canvas ratio for each figure.
#   - Main figures use symmetric/orderly layouts; Figure 4 is relaid out as a clean 2x4 sequence.
#   - Preserve the approved single-panel TIFFs; do not redraw scientific data.
#   - Increase panel-letter size, reduce blank space, reduce gutters, and allow
#     dense panels to occupy larger areas.
#
# INPUT
#   A folder containing:
#     Figure1_A.tif ... Figure6_H.tif
#     Supplementary_Figure_S1_A.tif ... Supplementary_Figure_S8_D3_Reactome.tif
#
# OUTPUT
#   16_COMPOSITE_FIGURES_FINAL_F2B_FIG4/
#     main_tif/
#     supplementary_tif/
#     preview_png/
#     logs/
#
# TIFF SETTINGS
#   600 dpi, LZW compression, white background.
# ============================================================================

options(stringsAsFactors = FALSE, scipen = 999)

# ----------------------------- USER SETTINGS ---------------------------------
USER_PROJECT_ROOT <- "C:/Users/33652/Desktop/LZZ/1、yssj"
OUTPUT_FOLDER_NAME <- "__TEMP_KELOID_COMPOSITE_ENGINE__"

# The script searches these folders first. Add your own path here if needed.
PANEL_DIR_CANDIDATES <- c(
  file.path(USER_PROJECT_ROOT, "15_PUBLICATION_PANELS_FINAL_F2B_DISEASE_PROJECTION", "individual_panel_tif"),
  file.path(USER_PROJECT_ROOT, "15_PUBLICATION_PANELS_ONLY_V1_2_FIXED", "individual_panel_tif"),
  file.path(USER_PROJECT_ROOT, "15_PUBLICATION_PANELS_ONLY_V1_1_FIXED", "individual_panel_tif"),
  file.path(USER_PROJECT_ROOT, "15_PUBLICATION_PANELS_ONLY_FINAL", "individual_panel_tif"),
  file.path(USER_PROJECT_ROOT, "15_PUBLICATION_PANELS_FINAL", "individual_panel_tif"),
  file.path(USER_PROJECT_ROOT, "individual_panel_tif"),
  USER_PROJECT_ROOT
)

DPI <- 600
SAVE_PREVIEW_PNG <- FALSE
PREVIEW_DPI <- 220
AUTO_INSTALL <- TRUE

# Larger, clearer panel letters. If you still feel A/B/C are small, increase to 92-100.
DEFAULT_PANEL_LETTER_SIZE <- 82L
SUPP_PANEL_LETTER_SIZE <- 78L

# Small internal panel padding. Lower values make panels larger.
DEFAULT_PANEL_PAD_PX <- 10L
OUTER_BORDER_PX <- 44L

FONT_FAMILY <- if (.Platform$OS.type == "windows") "Arial" else "sans"

# ----------------------------- PACKAGES --------------------------------------
required_pkgs <- c("magick")

install_and_load <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0L) {
    if (!AUTO_INSTALL) stop("Missing package(s): ", paste(missing, collapse = ", "))
    install.packages(missing, repos = "https://cloud.r-project.org")
  }
  invisible(lapply(pkgs, library, character.only = TRUE))
}
install_and_load(required_pkgs)

# ----------------------------- PATH UTILITIES --------------------------------
normalize2 <- function(x) normalizePath(x, winslash = "/", mustWork = FALSE)

resolve_project_root <- function() {
  candidates <- unique(c(
    USER_PROJECT_ROOT,
    getwd(),
    dirname(getwd()),
    file.path(Sys.getenv("USERPROFILE"), "Desktop", "LZZ", "1、yssj")
  ))
  candidates <- candidates[nzchar(candidates)]
  candidates <- candidates[dir.exists(candidates)]
  if (!length(candidates)) return(normalize2(USER_PROJECT_ROOT))

  score_root <- function(root) {
    files <- tryCatch(
      list.files(root, recursive = TRUE, full.names = FALSE, include.dirs = FALSE),
      error = function(e) character()
    )
    anchors <- c(
      "Figure1_A.tif",
      "Figure6_H.tif",
      "Supplementary_Figure_S8_D3_Reactome.tif"
    )
    sum(tolower(basename(files)) %in% tolower(anchors))
  }

  scores <- vapply(candidates, score_root, numeric(1))
  best <- candidates[which.max(scores)]
  if (max(scores) == 0) best <- candidates[1]
  normalize2(best)
}

PROJECT_ROOT <- resolve_project_root()

OUTPUT_DIR <- file.path(PROJECT_ROOT, OUTPUT_FOLDER_NAME)
DIR_MAIN <- file.path(OUTPUT_DIR, "main_tif")
DIR_SUPP <- file.path(OUTPUT_DIR, "supplementary_tif")
DIR_PREV <- file.path(OUTPUT_DIR, "preview_png")
DIR_LOG  <- file.path(OUTPUT_DIR, "logs")
for (d in c(OUTPUT_DIR, DIR_MAIN, DIR_SUPP, DIR_PREV, DIR_LOG)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

resolve_panel_dir <- function(root, candidates) {
  for (p in candidates) {
    pp <- normalize2(p)
    if (dir.exists(pp)) {
      files <- list.files(pp, pattern = "\\.tif$", full.names = FALSE, ignore.case = TRUE)
      if (
        any(grepl("^Figure1_A\\.tif$", files, ignore.case = TRUE)) &&
        any(grepl("^Supplementary_Figure_S1_A\\.tif$", files, ignore.case = TRUE))
      ) {
        return(pp)
      }
    }
  }

  # Fallback: recursively search under project root only.
  tif_files <- tryCatch(
    list.files(root, recursive = TRUE, full.names = TRUE, pattern = "\\.tif$", ignore.case = TRUE),
    error = function(e) character()
  )
  all_dirs <- unique(dirname(tif_files))
  for (pp in all_dirs) {
    files <- basename(list.files(pp, pattern = "\\.tif$", full.names = FALSE, ignore.case = TRUE))
    if (
      any(grepl("^Figure1_A\\.tif$", files, ignore.case = TRUE)) &&
      any(grepl("^Supplementary_Figure_S1_A\\.tif$", files, ignore.case = TRUE))
    ) {
      return(normalize2(pp))
    }
  }

  stop(
    "Could not find individual panel TIFFs. Please edit PANEL_DIR_CANDIDATES or put the panel TIFFs under the project root."
  )
}

PANEL_DIR <- resolve_panel_dir(PROJECT_ROOT, PANEL_DIR_CANDIDATES)

cat("PROJECT_ROOT = ", PROJECT_ROOT, "\n", file = file.path(DIR_LOG, "00_run_info.txt"))
cat("PANEL_DIR    = ", PANEL_DIR,    "\n", file = file.path(DIR_LOG, "00_run_info.txt"), append = TRUE)
cat("OUTPUT_DIR   = ", OUTPUT_DIR,   "\n", file = file.path(DIR_LOG, "00_run_info.txt"), append = TRUE)

# ----------------------------- IMAGE UTILITIES -------------------------------
read_panel <- function(filename) {
  path <- file.path(PANEL_DIR, filename)
  if (!file.exists(path)) stop("Missing panel file: ", path)
  img <- magick::image_read(path)
  magick::image_background(img, "white", flatten = TRUE)
}

# Conservative trimming: removes excessive white border from the single-panel TIFF
# but adds a small white border back to avoid cropping labels.
trim_panel <- function(img) {
  out <- tryCatch(
    magick::image_trim(img, fuzz = 1),
    error = function(e) img
  )
  out <- magick::image_border(out, color = "white", geometry = "10x10")
  out
}

scale_to_fit <- function(img, target_w_px, target_h_px) {
  info <- magick::image_info(img)[1, ]
  scale <- min(target_w_px / info$width, target_h_px / info$height)
  new_w <- max(1L, floor(info$width  * scale))
  new_h <- max(1L, floor(info$height * scale))
  magick::image_resize(img, paste0(new_w, "x", new_h, "!"))
}

add_letter <- function(canvas, label, x_px, y_px, pointsize) {
  magick::image_annotate(
    canvas,
    text = label,
    size = pointsize,
    font = FONT_FAMILY,
    weight = 700,
    color = "black",
    location = paste0("+", x_px, "+", y_px),
    gravity = "northwest"
  )
}

place_panel <- function(
  canvas,
  filename,
  x,
  y,
  w,
  h,
  label = NULL,
  panel_pad_px = DEFAULT_PANEL_PAD_PX,
  label_dx = 0L,
  label_dy = 0L,
  label_size = DEFAULT_PANEL_LETTER_SIZE
) {
  stopifnot(all(c(x, y, w, h) >= 0), x + w <= 1.0001, y + h <= 1.0001)

  info <- magick::image_info(canvas)[1, ]
  canvas_w <- info$width
  canvas_h <- info$height

  box_x <- round(x * canvas_w)
  box_y <- round(y * canvas_h)
  box_w <- round(w * canvas_w)
  box_h <- round(h * canvas_h)

  img <- read_panel(filename)
  img <- trim_panel(img)
  img <- scale_to_fit(
    img,
    max(1L, box_w - 2L * panel_pad_px),
    max(1L, box_h - 2L * panel_pad_px)
  )
  ii <- magick::image_info(img)[1, ]

  # Center panel image in its allocated area.
  x_off <- box_x + floor((box_w - ii$width) / 2)
  y_off <- box_y + floor((box_h - ii$height) / 2)
  canvas <- magick::image_composite(canvas, img, operator = "over", offset = paste0("+", x_off, "+", y_off))

  if (!is.null(label) && nzchar(label)) {
    canvas <- add_letter(
      canvas,
      label = label,
      x_px = box_x + label_dx,
      y_px = box_y + label_dy,
      pointsize = label_size
    )
  }

  canvas
}

panel_df <- function(
  label,
  file,
  x,
  y,
  w,
  h,
  panel_pad_px = DEFAULT_PANEL_PAD_PX,
  label_dx = 0L,
  label_dy = 0L,
  label_size = DEFAULT_PANEL_LETTER_SIZE
) {
  data.frame(
    label = label,
    file = file,
    x = x,
    y = y,
    w = w,
    h = h,
    panel_pad_px = panel_pad_px,
    label_dx = label_dx,
    label_dy = label_dy,
    label_size = label_size,
    stringsAsFactors = FALSE
  )
}

compose_figure <- function(figure_id, panels, width_in, height_in, out_tif, out_png = NULL) {
  width_px  <- round(width_in  * DPI)
  height_px <- round(height_in * DPI)

  canvas <- magick::image_blank(width = width_px, height = height_px, color = "white")

  for (i in seq_len(nrow(panels))) {
    canvas <- place_panel(
      canvas = canvas,
      filename = panels$file[i],
      x = panels$x[i],
      y = panels$y[i],
      w = panels$w[i],
      h = panels$h[i],
      label = panels$label[i],
      panel_pad_px = panels$panel_pad_px[i],
      label_dx = panels$label_dx[i],
      label_dy = panels$label_dy[i],
      label_size = panels$label_size[i]
    )
  }

  # Remove global external blank space and add a consistent outside border.
  canvas <- tryCatch(magick::image_trim(canvas, fuzz = 1), error = function(e) canvas)
  canvas <- magick::image_border(canvas, color = "white", geometry = paste0(OUTER_BORDER_PX, "x", OUTER_BORDER_PX))

  magick::image_write(
    canvas,
    path = out_tif,
    format = "tiff",
    compression = "LZW",
    density = paste0(DPI, "x", DPI)
  )

  if (!is.null(out_png)) {
    magick::image_write(
      canvas,
      path = out_png,
      format = "png",
      density = paste0(PREVIEW_DPI, "x", PREVIEW_DPI)
    )
  }

  invisible(canvas)
}

# ----------------------------- LAYOUTS ---------------------------------------
# The coordinates below fill the page vertically. Figure size is optimized per
# figure instead of being forced to the same height.
#
# x, y, w, h are relative coordinates on the page.
# Increase w/h for dense panels; decrease white space by keeping rows filled.

figure_layouts <- list(
  # --------------------------------------------------------------------------
  # Main figures: symmetric/orderly layout version.
  # Rule: avoid one panel alone in a full row; panels are arranged in sequence.
  # --------------------------------------------------------------------------

  Figure1 = list(
    width_in = 8.0, height_in = 8.4, dir = DIR_MAIN,
    panels = rbind(
      # A-B
      panel_df("A", "Figure1_A.tif", 0.00, 0.00, 0.50, 0.32),
      panel_df("B", "Figure1_B.tif", 0.50, 0.00, 0.50, 0.32),
      # C-D-E
      panel_df("C", "Figure1_C.tif", 0.00, 0.32, 0.34, 0.34),
      panel_df("D", "Figure1_D.tif", 0.34, 0.32, 0.33, 0.34),
      panel_df("E", "Figure1_E.tif", 0.67, 0.32, 0.33, 0.34),
      # F-G
      panel_df("F", "Figure1_F.tif", 0.00, 0.66, 0.50, 0.34),
      panel_df("G", "Figure1_G.tif", 0.50, 0.66, 0.50, 0.34)
    )
  ),

  Figure2 = list(
    width_in = 8.0, height_in = 8.8, dir = DIR_MAIN,
    panels = rbind(
      # A-B
      panel_df("A", "Figure2_A.tif", 0.00, 0.00, 0.47, 0.34),
      panel_df("B", "Figure2_B.tif", 0.47, 0.00, 0.53, 0.34),
      # C-D
      panel_df("C", "Figure2_C.tif", 0.00, 0.34, 0.50, 0.33),
      panel_df("D", "Figure2_D.tif", 0.50, 0.34, 0.50, 0.33),
      # E-F
      panel_df("E", "Figure2_E.tif", 0.00, 0.67, 0.50, 0.33),
      panel_df("F", "Figure2_F.tif", 0.50, 0.67, 0.50, 0.33)
    )
  ),

  Figure3 = list(
    width_in = 11.2, height_in = 6.8, dir = DIR_MAIN,
    panels = rbind(
      # A-B-C-D
      panel_df("A", "Figure3_A.tif", 0.00, 0.00, 0.25, 0.50),
      panel_df("B", "Figure3_B.tif", 0.25, 0.00, 0.25, 0.50),
      panel_df("C", "Figure3_C.tif", 0.50, 0.00, 0.25, 0.50),
      panel_df("D", "Figure3_D.tif", 0.75, 0.00, 0.25, 0.50),
      # E-F-G-H
      panel_df("E", "Figure3_E.tif", 0.00, 0.50, 0.25, 0.50),
      panel_df("F", "Figure3_F.tif", 0.25, 0.50, 0.25, 0.50),
      panel_df("G", "Figure3_G.tif", 0.50, 0.50, 0.25, 0.50),
      panel_df("H", "Figure3_H.tif", 0.75, 0.50, 0.25, 0.50)
    )
  ),

  Figure4 = list(
    width_in = 13.2, height_in = 7.2, dir = DIR_MAIN,
    panels = rbind(
      # Clean 2x4 layout in strict panel order
      panel_df("A", "Figure4_A.tif", 0.00, 0.00, 0.25, 0.50),
      panel_df("B", "Figure4_B.tif", 0.25, 0.00, 0.25, 0.50),
      panel_df("C", "Figure4_C.tif", 0.50, 0.00, 0.25, 0.50),
      panel_df("D", "Figure4_D.tif", 0.75, 0.00, 0.25, 0.50),
      panel_df("E", "Figure4_E.tif", 0.00, 0.50, 0.25, 0.50),
      panel_df("F", "Figure4_F.tif", 0.25, 0.50, 0.25, 0.50),
      panel_df("G", "Figure4_G.tif", 0.50, 0.50, 0.25, 0.50),
      panel_df("H", "Figure4_H.tif", 0.75, 0.50, 0.25, 0.50)
    )
  ),

  Figure5 = list(
    width_in = 10.8, height_in = 7.2, dir = DIR_MAIN,
    panels = rbind(
      # A-B-C
      panel_df("A", "Figure5_A.tif", 0.00, 0.00, 0.33, 0.48),
      panel_df("B", "Figure5_B.tif", 0.33, 0.00, 0.34, 0.48),
      panel_df("C", "Figure5_C.tif", 0.67, 0.00, 0.33, 0.48),
      # D-E-F; F is the visual core and gets double width
      panel_df("D", "Figure5_D.tif", 0.00, 0.48, 0.25, 0.52),
      panel_df("E", "Figure5_E.tif", 0.25, 0.48, 0.25, 0.52),
      panel_df("F", "Figure5_F.tif", 0.50, 0.48, 0.50, 0.52)
    )
  ),

  Figure6 = list(
    width_in = 11.4, height_in = 8.6, dir = DIR_MAIN,
    panels = rbind(
      # A-B-C
      panel_df("A", "Figure6_A.tif", 0.00, 0.00, 0.27, 0.33),
      panel_df("B", "Figure6_B.tif", 0.27, 0.00, 0.33, 0.33),
      panel_df("C", "Figure6_C.tif", 0.60, 0.00, 0.40, 0.33),
      # D-E-F
      panel_df("D", "Figure6_D.tif", 0.00, 0.33, 0.34, 0.33),
      panel_df("E", "Figure6_E.tif", 0.34, 0.33, 0.33, 0.33),
      panel_df("F", "Figure6_F.tif", 0.67, 0.33, 0.33, 0.33),
      # G-H; G gets wider space for pathway names
      panel_df("G", "Figure6_G.tif", 0.00, 0.66, 0.58, 0.34),
      panel_df("H", "Figure6_H.tif", 0.58, 0.66, 0.42, 0.34)
    )
  ),

  # --------------------------------------------------------------------------
  # Supplementary figures: keep the already approved supplementary layout.
  # --------------------------------------------------------------------------

  Supplementary_Figure_S1 = list(
    width_in = 7.4, height_in = 8.8, dir = DIR_SUPP,
    panels = rbind(
      panel_df("A", "Supplementary_Figure_S1_A.tif", 0.00, 0.00, 0.50, 0.33, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("B", "Supplementary_Figure_S1_B.tif", 0.50, 0.00, 0.50, 0.33, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("C", "Supplementary_Figure_S1_C.tif", 0.00, 0.33, 0.50, 0.33, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("D", "Supplementary_Figure_S1_D.tif", 0.50, 0.33, 0.50, 0.33, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("E", "Supplementary_Figure_S1_E.tif", 0.00, 0.66, 0.50, 0.34, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("F", "Supplementary_Figure_S1_F.tif", 0.50, 0.66, 0.50, 0.34, label_size = SUPP_PANEL_LETTER_SIZE)
    )
  ),

  Supplementary_Figure_S2 = list(
    width_in = 7.4, height_in = 8.3, dir = DIR_SUPP,
    panels = rbind(
      panel_df("A", "Supplementary_Figure_S2_A.tif", 0.00, 0.00, 0.50, 0.34, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("B", "Supplementary_Figure_S2_B.tif", 0.50, 0.00, 0.50, 0.34, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("C", "Supplementary_Figure_S2_C.tif", 0.00, 0.34, 1.00, 0.33, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("D", "Supplementary_Figure_S2_D.tif", 0.00, 0.67, 1.00, 0.33, label_size = SUPP_PANEL_LETTER_SIZE)
    )
  ),

  Supplementary_Figure_S3 = list(
    width_in = 7.4, height_in = 7.9, dir = DIR_SUPP,
    panels = rbind(
      panel_df("A", "Supplementary_Figure_S3_A.tif", 0.00, 0.00, 0.50, 0.50, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("B", "Supplementary_Figure_S3_B.tif", 0.50, 0.00, 0.50, 0.50, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("C", "Supplementary_Figure_S3_C.tif", 0.00, 0.50, 0.50, 0.50, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("D", "Supplementary_Figure_S3_D.tif", 0.50, 0.50, 0.50, 0.50, label_size = SUPP_PANEL_LETTER_SIZE)
    )
  ),

  Supplementary_Figure_S4 = list(
    width_in = 7.4, height_in = 7.9, dir = DIR_SUPP,
    panels = rbind(
      panel_df("A", "Supplementary_Figure_S4_A.tif", 0.00, 0.00, 0.50, 0.50, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("B", "Supplementary_Figure_S4_B.tif", 0.50, 0.00, 0.50, 0.50, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("C", "Supplementary_Figure_S4_C.tif", 0.00, 0.50, 0.50, 0.50, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("D", "Supplementary_Figure_S4_D.tif", 0.50, 0.50, 0.50, 0.50, label_size = SUPP_PANEL_LETTER_SIZE)
    )
  ),

  Supplementary_Figure_S5 = list(
    width_in = 7.4, height_in = 7.9, dir = DIR_SUPP,
    panels = rbind(
      panel_df("A", "Supplementary_Figure_S5_A.tif", 0.00, 0.00, 0.50, 0.53, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("B", "Supplementary_Figure_S5_B.tif", 0.50, 0.00, 0.50, 0.53, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("C", "Supplementary_Figure_S5_C.tif", 0.00, 0.53, 0.50, 0.47, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("D", "Supplementary_Figure_S5_D.tif", 0.50, 0.53, 0.50, 0.47, label_size = SUPP_PANEL_LETTER_SIZE)
    )
  ),

  Supplementary_Figure_S6 = list(
    width_in = 8.8, height_in = 10.5, dir = DIR_SUPP,
    panels = rbind(
      panel_df("A", "Supplementary_Figure_S6_A_Donor_Donor1.tif", 0.00, 0.00, 1.00, 0.18, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("B", "Supplementary_Figure_S6_A_Donor_Donor2.tif", 0.00, 0.18, 1.00, 0.18, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("C", "Supplementary_Figure_S6_A_Donor_Donor3.tif", 0.00, 0.36, 1.00, 0.18, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("D", "Supplementary_Figure_S6_A_Donor_Donor4.tif", 0.00, 0.54, 1.00, 0.18, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("E", "Supplementary_Figure_S6_B.tif", 0.00, 0.72, 0.55, 0.28, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("F", "Supplementary_Figure_S6_C.tif", 0.55, 0.72, 0.45, 0.28, label_size = SUPP_PANEL_LETTER_SIZE)
    )
  ),

  Supplementary_Figure_S7 = list(
    width_in = 7.4, height_in = 7.9, dir = DIR_SUPP,
    panels = rbind(
      panel_df("A", "Supplementary_Figure_S7_A.tif", 0.00, 0.00, 0.50, 0.50, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("B", "Supplementary_Figure_S7_B.tif", 0.50, 0.00, 0.50, 0.50, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("C", "Supplementary_Figure_S7_C.tif", 0.00, 0.50, 0.50, 0.50, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("D", "Supplementary_Figure_S7_D.tif", 0.50, 0.50, 0.50, 0.50, label_size = SUPP_PANEL_LETTER_SIZE)
    )
  ),

  Supplementary_Figure_S8 = list(
    width_in = 8.8, height_in = 11.8, dir = DIR_SUPP,
    panels = rbind(
      panel_df("A", "Supplementary_Figure_S8_A.tif", 0.00, 0.00, 0.50, 0.20, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("B", "Supplementary_Figure_S8_B.tif", 0.50, 0.00, 0.50, 0.20, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("C", "Supplementary_Figure_S8_C.tif", 0.00, 0.20, 1.00, 0.24, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("D", "Supplementary_Figure_S8_D1_Hallmark.tif", 0.00, 0.44, 1.00, 0.19, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("E", "Supplementary_Figure_S8_D2_GO_BP.tif", 0.00, 0.63, 1.00, 0.19, label_size = SUPP_PANEL_LETTER_SIZE),
      panel_df("F", "Supplementary_Figure_S8_D3_Reactome.tif", 0.00, 0.82, 1.00, 0.18, label_size = SUPP_PANEL_LETTER_SIZE)
    )
  )
)

# ----------------------------- VALIDATION ------------------------------------
all_required <- unique(unlist(lapply(figure_layouts, function(x) x$panels$file)))
missing_files <- all_required[!file.exists(file.path(PANEL_DIR, all_required))]
if (length(missing_files) > 0L) {
  writeLines(missing_files, file.path(DIR_LOG, "missing_panel_files.txt"))
  stop(
    "Missing ", length(missing_files), " panel TIFF file(s). See logs/missing_panel_files.txt"
  )
}

# ----------------------------- RENDER ALL ------------------------------------
manifest <- data.frame(
  figure_id = character(),
  output_tif = character(),
  output_preview = character(),
  width_in = numeric(),
  height_in = numeric(),
  stringsAsFactors = FALSE
)

for (fig_id in names(figure_layouts)) {
  spec <- figure_layouts[[fig_id]]
  out_tif <- file.path(spec$dir, paste0(fig_id, ".tif"))
  out_png <- if (SAVE_PREVIEW_PNG) file.path(DIR_PREV, paste0(fig_id, "_preview.png")) else NULL

  cat("Rendering ", fig_id, " ...\n", sep = "")
  compose_figure(
    figure_id = fig_id,
    panels = spec$panels,
    width_in = spec$width_in,
    height_in = spec$height_in,
    out_tif = out_tif,
    out_png = out_png
  )

  manifest <- rbind(
    manifest,
    data.frame(
      figure_id = fig_id,
      output_tif = normalize2(out_tif),
      output_preview = if (is.null(out_png)) "" else normalize2(out_png),
      width_in = spec$width_in,
      height_in = spec$height_in,
      stringsAsFactors = FALSE
    )
  )
}

write.csv(
  manifest,
  file.path(DIR_LOG, "01_output_manifest.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

write.csv(
  do.call(
    rbind,
    lapply(names(figure_layouts), function(nm) {
      cbind(figure_id = nm, figure_layouts[[nm]]$panels)
    })
  ),
  file.path(DIR_LOG, "02_layout_coordinates.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

writeLines(capture.output(sessionInfo()), file.path(DIR_LOG, "03_sessionInfo.txt"))
writeLines(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), file.path(DIR_LOG, "RUN_COMPLETED.txt"))

cat("\nAll optimized composite figures have been generated successfully.\n")
cat("Output folder:\n", OUTPUT_DIR, "\n", sep = "")




# ============================================================================
# PART 3/3: FINAL REARRANGED PUBLICATION FIGURES
# Version: FINAL STORY-LAYOUT V3 | 2026-09-08
#
# FINAL OUTPUT STRUCTURE
#   18_PUBLICATION_FIGURES_TIF_FINAL_V3/
#     main_panel_tif/           Figure1_A.tif ... Figure6_F/H.tif
#     main_composite_tif/       Figure1.tif ... Figure6.tif
#     supplementary_panel_tif/  Supplementary_Figure_S1_*.tif ... S9_*.tif
#     supplementary_composite_tif/ Supplementary_Figure_S1.tif ... S9.tif
#     source_data/              panel-level source data CSV files
#     source_data/raw_inputs_used/ copies of numeric input tables actually used
#     tables/                   main/supplementary manuscript tables as CSV
#     logs/                     manifest, warnings, sessionInfo
#
# FIGURE FILES ARE TIFF ONLY (600 dpi, LZW, white background).
# Numeric source data remain CSV because numerical data cannot meaningfully be
# represented as TIFF. No PNG/PDF/JPEG figure outputs are generated here.
# ============================================================================

message("\n============================================================")
message("Starting FINAL STORY-LAYOUT V3 rearrangement...")
message("============================================================")

FINAL_OUTPUT_FOLDER <- "18_PUBLICATION_FIGURES_TIF_FINAL_V3"
FINAL_DIR <- file.path(PROJECT_ROOT, FINAL_OUTPUT_FOLDER)
FINAL_MAIN_PANEL <- file.path(FINAL_DIR, "main_panel_tif")
FINAL_MAIN_COMP <- file.path(FINAL_DIR, "main_composite_tif")
FINAL_SUPP_PANEL <- file.path(FINAL_DIR, "supplementary_panel_tif")
FINAL_SUPP_COMP <- file.path(FINAL_DIR, "supplementary_composite_tif")
FINAL_SOURCE <- file.path(FINAL_DIR, "source_data")
FINAL_RAW <- file.path(FINAL_SOURCE, "raw_inputs_used")
FINAL_TABLES <- file.path(FINAL_DIR, "tables")
FINAL_LOGS <- file.path(FINAL_DIR, "logs")

if (dir.exists(FINAL_DIR)) unlink(FINAL_DIR, recursive = TRUE, force = TRUE)
for (d in c(FINAL_MAIN_PANEL, FINAL_MAIN_COMP, FINAL_SUPP_PANEL, FINAL_SUPP_COMP,
            FINAL_SOURCE, FINAL_RAW, FINAL_TABLES, FINAL_LOGS)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

FINAL_DPI <- 600
FINAL_WIDTH <- 7.2
FINAL_FONT <- BASE_FONT
FINAL_WARNINGS <- character()

final_warn <- function(...) {
  txt <- paste(..., collapse = " ")
  FINAL_WARNINGS <<- c(FINAL_WARNINGS, txt)
  warning(txt, call. = FALSE)
}

final_safe_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  x
}

save_final_tif <- function(plot, filename, width = 5.8, height = 4.6, dir = FINAL_MAIN_PANEL) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(dir, filename)
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_tiff(
      filename = path,
      width = width,
      height = height,
      units = "in",
      res = FINAL_DPI,
      compression = "lzw",
      background = "white"
    )
    print(plot)
    grDevices::dev.off()
  } else {
    grDevices::tiff(
      filename = path,
      width = width,
      height = height,
      units = "in",
      res = FINAL_DPI,
      compression = "lzw",
      bg = "white",
      type = if (.Platform$OS.type == "windows") "windows" else "cairo"
    )
    print(plot)
    grDevices::dev.off()
  }
  invisible(path)
}

FINAL_SOURCE_INDEX <- tibble::tibble(
  figure = character(), panel = character(), description = character(), file = character()
)

write_final_source <- function(x, figure, panel, description) {
  fn <- paste0(final_safe_name(figure), "_", final_safe_name(panel), "_source_data.csv")
  path <- file.path(FINAL_SOURCE, fn)
  utils::write.csv(as.data.frame(x), path, row.names = FALSE, fileEncoding = "UTF-8")
  FINAL_SOURCE_INDEX <<- dplyr::bind_rows(
    FINAL_SOURCE_INDEX,
    tibble::tibble(figure = figure, panel = panel, description = description, file = fn)
  )
  invisible(path)
}

# ----------------------------- common presentation ---------------------------
final_theme <- theme_pub(base_size = 9.2, title_size = 9.4) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 9.4, hjust = 0),
    axis.title = ggplot2::element_text(size = 8.8),
    axis.text = ggplot2::element_text(size = 8.0),
    legend.title = ggplot2::element_text(size = 8.0),
    legend.text = ggplot2::element_text(size = 7.6),
    strip.text = ggplot2::element_text(size = 8.0, face = "bold"),
    plot.margin = ggplot2::margin(7, 8, 7, 8)
  )

control_col <- "#687078"
keloid_col <- "#B2182B"
vu_col <- "#7B3294"
normal_skin_col <- "#8A8A8A"
stage_cols_final <- c(Skin = "#4D4D4D", Wound1 = "#D55E00", Wound7 = "#E69F00", Wound30 = "#0072B2")

stage_levels_final <- c("Skin", "Wound1", "Wound7", "Wound30")
stage_labels_final <- c(Skin = "Skin", Wound1 = "Day 1", Wound7 = "Day 7", Wound30 = "Day 30")

scatter_group_plot <- function(data, group_col, value_col, title, ylab,
                               levels, labels, fills, jitter_width = 0.07) {
  d <- data |>
    dplyr::mutate(
      .group = factor(.data[[group_col]], levels = levels, labels = labels),
      .value = .data[[value_col]]
    )
  ggplot2::ggplot(d, ggplot2::aes(.group, .value, fill = .group)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, colour = "#B0B0B0") +
    ggplot2::geom_point(
      shape = 21, colour = "black", stroke = 0.55, size = 3.1,
      position = ggplot2::position_jitter(width = jitter_width, height = 0, seed = 20260908)
    ) +
    ggplot2::stat_summary(fun = mean, geom = "point", shape = 95, size = 8.0, colour = "black") +
    ggplot2::scale_fill_manual(values = stats::setNames(fills, labels), guide = "none") +
    ggplot2::labs(title = title, x = NULL, y = ylab) +
    final_theme
}

heatmap_z4 <- function(data, id_col, group_col = NULL, title = NULL, order_ids = NULL) {
  score_cols <- c("z_Skin", "z_Wound1", "z_Wound7", "z_Wound30")
  d <- data
  if (is.null(order_ids)) order_ids <- as.character(d[[id_col]])
  d$.id <- factor(as.character(d[[id_col]]), levels = order_ids)
  long <- d |>
    dplyr::select(.id, dplyr::all_of(score_cols), dplyr::any_of(group_col)) |>
    tidyr::pivot_longer(dplyr::all_of(score_cols), names_to = "stage", values_to = "score") |>
    dplyr::mutate(
      stage = sub("^z_", "", stage),
      stage = factor(stage, levels = rev(stage_levels_final))
    )
  p <- ggplot2::ggplot(long, ggplot2::aes(.id, stage, fill = score)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.35) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, name = "z-score") +
    ggplot2::scale_y_discrete(labels = stage_labels_final) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    final_theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid = ggplot2::element_blank(),
      legend.position = "right"
    )
  list(plot = p, source = long)
}

# ----------------------------- Stage 13 inputs -------------------------------
find_final_input <- function(basename_target) {
  hits <- list.files(PROJECT_ROOT, recursive = TRUE, full.names = TRUE, include.dirs = FALSE)
  hits <- hits[tolower(basename(hits)) == tolower(basename_target)]
  hits <- hits[!grepl("__TEMP_KELOID|18_PUBLICATION_FIGURES_TIF_FINAL_V3", hits, fixed = FALSE)]
  if (!length(hits)) return(NA_character_)
  # Prefer Stage13 folder and latest file.
  score <- ifelse(grepl("13_STAGE13_SIGNATURE_ROBUSTNESS|第十三阶段", hits), 10, 0)
  ord <- order(score, file.info(hits)$mtime, decreasing = TRUE, na.last = TRUE)
  normalizePath(hits[ord[1]], winslash = "/", mustWork = TRUE)
}

stage13_logo_path <- find_final_input("06_leave_one_Wound30_gene_out_results.csv")
stage13_random_path <- find_final_input("07_random_80pct_Wound30_subsampling_results.csv")
stage13_summary_path <- find_final_input("09_signature_robustness_summary.csv")
stage13_influence_path <- find_final_input("10_Wound30_single_gene_influence_ranking.csv")
stage13_baseline_path <- find_final_input("03_baseline_full_signature_reproduction.csv")

if (any(is.na(c(stage13_logo_path, stage13_random_path, stage13_summary_path,
                stage13_influence_path, stage13_baseline_path)))) {
  stop("Stage13 output CSV files are required for the final layout. Please keep 13_STAGE13_SIGNATURE_ROBUSTNESS under PROJECT_ROOT.")
}

stage13_logo <- readr::read_csv(stage13_logo_path, show_col_types = FALSE)
stage13_random <- readr::read_csv(stage13_random_path, show_col_types = FALSE)
stage13_summary <- readr::read_csv(stage13_summary_path, show_col_types = FALSE)
stage13_influence <- readr::read_csv(stage13_influence_path, show_col_types = FALSE)
stage13_baseline <- readr::read_csv(stage13_baseline_path, show_col_types = FALSE)

# ============================================================================
# FIGURE 1 — DETAILED GRAPHICAL SUMMARY / WORKFLOW
# ============================================================================
message("Building final Figure 1...")

# 1A Study question
f1a_df <- data.frame(
  x = c(1,2,3,4,5.15),
  y = 1,
  label = c("Skin", "Day 1", "Day 7", "Day 30", "Keloid"),
  type = c("stage","stage","stage","stage","disease")
)
write_final_source(f1a_df, "Figure1", "A", "Study-question timeline and disease projection target")
pf1a <- ggplot2::ggplot(f1a_df, ggplot2::aes(x, y)) +
  ggplot2::geom_segment(data = data.frame(x=1,xend=4,y=1,yend=1), ggplot2::aes(x=x,xend=xend,y=y,yend=yend),
                        inherit.aes = FALSE, linewidth = 0.8, colour = "#666666", arrow = grid::arrow(length = grid::unit(0.08,"in"))) +
  ggplot2::geom_point(ggplot2::aes(fill = label), shape = 21, size = 5.2, colour = "black", stroke = 0.6) +
  ggplot2::geom_text(ggplot2::aes(label = label), vjust = -1.35, size = 3.0, family = FINAL_FONT) +
  ggplot2::annotate("segment", x=4.25, xend=4.85, y=1, yend=1, linewidth=0.7, linetype="dashed", colour="#777777") +
  ggplot2::annotate("text", x=3, y=1.55, label="Where do keloid fibroblasts map\nalong normal wound progression?", size=3.2, fontface="bold", family=FINAL_FONT) +
  ggplot2::scale_fill_manual(values = c("Skin"=stage_cols_final["Skin"],"Day 1"=stage_cols_final["Wound1"],"Day 7"=stage_cols_final["Wound7"],"Day 30"=stage_cols_final["Wound30"],"Keloid"=keloid_col), guide="none") +
  ggplot2::coord_cartesian(xlim=c(0.6,5.6),ylim=c(0.5,1.8),clip="off") +
  ggplot2::theme_void(base_family = FINAL_FONT) +
  ggplot2::theme(plot.margin=ggplot2::margin(12,8,10,8))

# 1B Dataset roles
f1b_df <- data.frame(
  dataset = c("GSE241132","GSE163973","GSE181316","GSE220300","GSE241124","GSE265972","GSE181297"),
  role = c("Normal-wound reference","Discovery","Independent replication","External keloid support","Spatial support","Venous-ulcer control","Exploratory / mixed"),
  x = 1:7
)
write_final_source(f1b_df, "Figure1", "B", "Datasets and analytical roles")
pf1b <- ggplot2::ggplot(f1b_df, ggplot2::aes(x, 1)) +
  ggplot2::geom_tile(ggplot2::aes(fill = role), width = 0.90, height = 0.72, colour="white", linewidth=0.7) +
  ggplot2::geom_text(ggplot2::aes(label = paste0(dataset,"\n",role)), size=2.45, lineheight=0.95, family=FINAL_FONT, fontface="bold") +
  ggplot2::scale_fill_manual(values = c(
    "Normal-wound reference"="#D9E8F5","Discovery"="#F6D6D4","Independent replication"="#F2B9B5",
    "External keloid support"="#F7E3CF","Spatial support"="#D8ECE8","Venous-ulcer control"="#E8DCF1","Exploratory / mixed"="#E2E2E2"
  ), guide="none") +
  ggplot2::labs(title="Public datasets assigned distinct analytical roles") +
  ggplot2::coord_cartesian(ylim=c(0.55,1.45),clip="off") +
  ggplot2::theme_void(base_family=FINAL_FONT) +
  ggplot2::theme(plot.title=ggplot2::element_text(size=9.2,face="bold"), plot.margin=ggplot2::margin(8,5,8,5))

# 1C analysis path
f1c_df <- data.frame(
  x = 1:5,
  label = c("Fibroblast\nprofiles","Donor / patient\npseudobulk","Frozen wound-state\nsignatures","Disease\nprojection","Exact patient-level\ninference")
)
write_final_source(f1c_df, "Figure1", "C", "Primary analysis path")
pf1c <- ggplot2::ggplot(f1c_df, ggplot2::aes(x,1)) +
  ggplot2::geom_segment(data=data.frame(x=1.2,xend=4.8,y=1,yend=1),ggplot2::aes(x=x,xend=xend,y=y,yend=yend),inherit.aes=FALSE,
                        linewidth=0.7,colour="#777777",arrow=grid::arrow(length=grid::unit(0.07,"in"))) +
  ggplot2::geom_label(ggplot2::aes(label=label), size=2.7, family=FINAL_FONT, label.size=0.25, label.r=grid::unit(2,"pt"), fill="white") +
  ggplot2::labs(title="Biological units are preserved throughout the analysis") +
  ggplot2::coord_cartesian(xlim=c(0.6,5.4),ylim=c(0.65,1.35),clip="off") +
  ggplot2::theme_void(base_family=FINAL_FONT) +
  ggplot2::theme(plot.title=ggplot2::element_text(size=9.2,face="bold"),plot.margin=ggplot2::margin(8,5,8,5))

# 1D endpoint definition
f1d_df <- data.frame(term=c("zWound30","mean(zWound1,zWound7)"), sign=c(1,-1))
write_final_source(f1d_df, "Figure1", "D", "Definition of the late-remodeling state")
pf1d <- ggplot2::ggplot() +
  ggplot2::annotate("text",x=0.5,y=0.70,label="Late-remodeling state",size=4.0,fontface="bold",family=FINAL_FONT) +
  ggplot2::annotate("text",x=0.5,y=0.48,label="=  zWound30  −  mean(zWound1, zWound7)",size=3.5,family=FINAL_FONT) +
  ggplot2::annotate("text",x=0.5,y=0.25,label="Prioritized in discovery → frozen before independent replication",size=2.8,colour="#555555",family=FINAL_FONT) +
  ggplot2::xlim(0,1)+ggplot2::ylim(0,1)+ggplot2::theme_void(base_family=FINAL_FONT)

# 1E core numeric results
disc_primary_f1 <- s5_effects |>
  dplyr::filter(cell_type == "Fibroblast", metric == "remodeling_completion") |>
  dplyr::slice(1)
rep_primary_f1 <- s6_effects |>
  dplyr::filter(analysis_set == "high_specificity", metric == "late_remodeling_state") |>
  dplyr::slice(1)
f1e_df <- data.frame(
  cohort=c("GSE163973\nDiscovery","GSE181316\nReplication","Integrated"),
  difference=c(disc_primary_f1$mean_difference[1], rep_primary_f1$mean_difference[1], stage13_baseline$integrated_mean_difference[1]),
  exact_p=c(disc_primary_f1$exact_permutation_p[1], rep_primary_f1$exact_two_sided_p[1], stage13_baseline$integrated_exact_two_sided_p[1])
)
write_final_source(f1e_df, "Figure1", "E", "Core late-remodeling effects")
pf1e <- ggplot2::ggplot(f1e_df,ggplot2::aes(cohort,difference,fill=cohort))+
  ggplot2::geom_col(width=0.62,colour="black",linewidth=0.45)+
  ggplot2::geom_text(ggplot2::aes(label=sprintf("%.3f",difference)),vjust=-0.5,size=3.0,family=FINAL_FONT)+
  ggplot2::scale_fill_manual(values=c("GSE163973\nDiscovery"="#C95A52","GSE181316\nReplication"="#D8817B","Integrated"="#8C2D24"),guide="none")+
  ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0.02,0.16)))+
  ggplot2::labs(title="Late-remodeling effect reproduces across core cohorts",x=NULL,y="Mean difference")+
  final_theme+ggplot2::theme(axis.text.x=ggplot2::element_text(size=7.5))

# 1F disease-control failure mode plane (actual aggregate effects)
f1f_df <- failure_plane |>
  dplyr::mutate(label=dplyr::case_when(
    grepl("keloid",tolower(disease)) ~ "Keloid",
    grepl("venous",tolower(disease)) ~ "Venous ulcer",
    TRUE ~ disease
  ))
write_final_source(f1f_df, "Figure1", "F", "Observed cross-disease early-versus-late effect coordinates")
pf1f <- ggplot2::ggplot(f1f_df,ggplot2::aes(early_state_effect,late_remodeling_effect,fill=label))+
  ggplot2::geom_hline(yintercept=0,colour="#AAAAAA",linewidth=0.35)+ggplot2::geom_vline(xintercept=0,colour="#AAAAAA",linewidth=0.35)+
  ggplot2::geom_abline(slope=1,intercept=0,linetype="dotted",colour="#777777")+
  ggplot2::geom_point(shape=21,size=4.5,colour="black",stroke=0.65)+
  ggrepel::geom_text_repel(ggplot2::aes(label=label),size=3.0,family=FINAL_FONT,show.legend=FALSE)+
  ggplot2::scale_fill_manual(values=c("Keloid"=keloid_col,"Venous ulcer"=vu_col),guide="none")+
  ggplot2::labs(title="Different chronic lesions show different failure balances",x="Early-state effect",y="Late-remodeling effect")+final_theme

# 1G gene-program concordance
f1g_df <- data.frame(program=c("Wound30\nup in both","Skin\ndown in both"),concordant=c(51,29),measurable=c(52,34)) |>
  dplyr::mutate(prop=concordant/measurable,label=paste0(concordant,"/",measurable))
write_final_source(f1g_df, "Figure1", "G", "Cross-cohort frozen-signature directional concordance")
pf1g <- ggplot2::ggplot(f1g_df,ggplot2::aes(program,prop,fill=program))+
  ggplot2::geom_col(width=0.62,colour="black",linewidth=0.45)+
  ggplot2::geom_text(ggplot2::aes(label=label),vjust=-0.45,size=3.1,family=FINAL_FONT)+
  ggplot2::scale_fill_manual(values=c("Wound30\nup in both"=stage_cols_final["Wound30"],"Skin\ndown in both"=stage_cols_final["Skin"]),guide="none")+
  ggplot2::scale_y_continuous(labels=scales::percent,limits=c(0,1.10),expand=c(0,0))+
  ggplot2::labs(title="Frozen wound-state genes show cross-cohort concordance",x=NULL,y="Concordant fraction")+final_theme

# 1H robustness summary
f1h_df <- stage13_summary |>
  dplyr::transmute(
    analysis=dplyr::case_when(grepl("leave",analysis,ignore.case=TRUE) ~ "Leave-one-gene-out",TRUE ~ "Random 80% subsampling"),
    direction=100*fraction_both_cohorts_expected_direction,
    separation=100*fraction_both_cohorts_complete_separation,
    effect80=100*fraction_integrated_effect_ge_80pct_baseline
  ) |>
  tidyr::pivot_longer(c(direction,separation,effect80),names_to="criterion",values_to="percent") |>
  dplyr::mutate(criterion=dplyr::recode(criterion,direction="Expected direction",separation="Complete separation",effect80=">=80% baseline effect"))
write_final_source(f1h_df, "Figure1", "H", "Stage13 signature robustness summary")
pf1h <- ggplot2::ggplot(f1h_df,ggplot2::aes(criterion,percent,fill=analysis))+
  ggplot2::geom_col(position=ggplot2::position_dodge(width=0.72),width=0.62,colour="black",linewidth=0.3)+
  ggplot2::scale_y_continuous(limits=c(0,105),breaks=c(0,50,100),labels=function(x)paste0(x,"%"),expand=c(0,0))+
  ggplot2::scale_fill_manual(values=c("Leave-one-gene-out"="#4C78A8","Random 80% subsampling"="#72B7B2"))+
  ggplot2::labs(title="Late-remodeling signal is not driven by one marker gene",x=NULL,y="Iterations meeting criterion",fill=NULL)+
  final_theme+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=25,hjust=1),legend.position="top")

fig1_final <- (pf1a | pf1b) / (pf1c | pf1d) / (pf1e | pf1f) / (pf1g | pf1h) +
  patchwork::plot_annotation(tag_levels="A", theme=theme_patch_tags) +
  patchwork::plot_layout(heights=c(1.0,0.9,1.0,1.0))

for (x in list(A=pf1a,B=pf1b,C=pf1c,D=pf1d,E=pf1e,F=pf1f,G=pf1g,H=pf1h)) {
  # loop handled below with explicit names
}
fig1_panels <- list(A=pf1a,B=pf1b,C=pf1c,D=pf1d,E=pf1e,F=pf1f,G=pf1g,H=pf1h)
fig1_dims <- list(A=c(6.0,3.5),B=c(7.0,3.7),C=c(6.4,3.4),D=c(5.7,3.4),E=c(5.6,4.2),F=c(5.6,4.4),G=c(5.3,4.1),H=c(6.0,4.3))
for (nm in names(fig1_panels)) save_final_tif(fig1_panels[[nm]],paste0("Figure1_",nm,".tif"),fig1_dims[[nm]][1],fig1_dims[[nm]][2],FINAL_MAIN_PANEL)
save_final_tif(fig1_final,"Figure1.tif",FINAL_WIDTH,11.2,FINAL_MAIN_COMP)

# ============================================================================
# FIGURE 2 — NORMAL-WOUND REFERENCE
# ============================================================================
message("Building final Figure 2...")
pf2a <- p1a
pf2b <- p1b
pf2c <- p1c
pf2d <- p1d
pf2e <- p1f

lodo_fib_final <- lodo_pred |>
  dplyr::filter(cell_type == "Fibroblast") |>
  dplyr::mutate(true_stage=factor(true_stage,levels=stage_levels_final))
write_final_source(lodo_fib_final, "Figure2", "F", "Donor-wise held-out ordinal trajectories")
pf2f <- ggplot2::ggplot(lodo_fib_final,ggplot2::aes(true_stage,predicted_position,group=held_out_donor,colour=held_out_donor))+
  ggplot2::geom_line(linewidth=0.7)+ggplot2::geom_point(size=2.8)+
  ggplot2::scale_x_discrete(labels=stage_labels_final)+
  ggplot2::labs(title="Held-out donor trajectories follow wound-stage order",x=NULL,y="Predicted ordinal position",colour="Donor")+
  final_theme+ggplot2::theme(legend.position="top",axis.text.x=ggplot2::element_text(angle=25,hjust=1))

# new source mappings for reused panels
write_final_source(ref_cell_counts |> dplyr::filter(main_cell_type=="Fibroblast"),"Figure2","A","Fibroblast yield across reference donor-stage samples")
write_final_source(ref_calibration |> dplyr::filter(cell_type=="Fibroblast"),"Figure2","B","Reference state-space coordinates derived from frozen stage scores")
write_final_source(ref_calibration |> dplyr::filter(cell_type=="Fibroblast") |> dplyr::select(sample_id,donor,true_stage,z_Skin,z_Wound1,z_Wound7,z_Wound30),"Figure2","C","Frozen reference stage scores")
write_final_source(confusion |> dplyr::filter(cell_type=="Fibroblast"),"Figure2","D","Leave-one-donor-out confusion matrix")
write_final_source(lodo_fib_final,"Figure2","E","True versus predicted ordinal wound-state positions")

fig2_final <- (pf2a | pf2b) / (pf2c | pf2d) / (pf2e | pf2f) +
  patchwork::plot_annotation(tag_levels="A",theme=theme_patch_tags)
fig2_panels <- list(A=pf2a,B=pf2b,C=pf2c,D=pf2d,E=pf2e,F=pf2f)
for (nm in names(fig2_panels)) save_final_tif(fig2_panels[[nm]],paste0("Figure2_",nm,".tif"),6.0,4.8,FINAL_MAIN_PANEL)
save_final_tif(fig2_final,"Figure2.tif",FINAL_WIDTH,9.6,FINAL_MAIN_COMP)

# ============================================================================
# FIGURE 3 — DISCOVERY + INDEPENDENT REPLICATION + EXACT INTEGRATION
# ============================================================================
message("Building final Figure 3...")
pf3a <- p2a
pf3b <- p2b
pf3c <- p2d

write_final_source(s5_fib |> dplyr::select(patient_id,group,z_Skin),"Figure3","D","Discovery-cohort patient-level Skin scores")
pf3d <- scatter_group_plot(
  s5_fib,"group","z_Skin","Discovery: loss of Skin-like homeostasis","Skin score",
  levels=c("normal_scar","keloid"),labels=c("Normal scar","Keloid"),fills=c(control_col,keloid_col)
)

pf3e <- p3a
pf3f <- p3c
pf3g <- p3d
pf3h <- p3g

write_final_source(s5_fib |> dplyr::select(patient_id,group,z_Skin,z_Wound1,z_Wound7,z_Wound30),"Figure3","A","Discovery frozen stage scores")
write_final_source(s5_fib |> dplyr::select(patient_id,group,late_remodeling_state),"Figure3","C","Discovery late-remodeling state")
write_final_source(s6_primary |> dplyr::select(patient_id,group,z_Skin,z_Wound1,z_Wound7,z_Wound30),"Figure3","E","Independent replication frozen stage scores")
write_final_source(s6_primary |> dplyr::select(patient_id,group,late_remodeling_state),"Figure3","F","Independent replication late-remodeling state")
write_final_source(core_patient,"Figure3","G","All 12 core patient-level observations")
write_final_source(core_integrated,"Figure3","H","Stratified exact integrated wound-state effects")

fig3_final <- (pf3a | pf3b | pf3c | pf3d) / (pf3e | pf3f | pf3g | pf3h) +
  patchwork::plot_annotation(tag_levels="A",theme=theme_patch_tags)
fig3_panels <- list(A=pf3a,B=pf3b,C=pf3c,D=pf3d,E=pf3e,F=pf3f,G=pf3g,H=pf3h)
for (nm in names(fig3_panels)) save_final_tif(fig3_panels[[nm]],paste0("Figure3_",nm,".tif"),5.8,4.7,FINAL_MAIN_PANEL)
save_final_tif(fig3_final,"Figure3.tif",FINAL_WIDTH,8.3,FINAL_MAIN_COMP)

# ============================================================================
# FIGURE 4 — EXTERNAL KELIOD SUPPORT + SPATIAL TRANSFERABILITY
# ============================================================================
message("Building final Figure 4...")
pf4a <- p4b
pf4b <- p4a
pf4c <- p4c
pf4d <- p4f
pf4e <- p4d
pf4f <- p4e

write_final_source(s7_support_group,"Figure4","A","External GSE220300 lesion-versus-scar late-remodeling scores")
write_final_source(s7_high |> dplyr::select(sample_key,patient_id,aggregate_state,z_Skin,z_Wound1,z_Wound7,z_Wound30),"Figure4","B","External GSE220300 frozen stage scores")
write_final_source(cross_support,"Figure4","C","Cross-cohort direction summary")
write_final_source(spot_rep,"Figure4","D",paste0("Representative GSE241124 spatial spots from ",rep_donor))
write_final_source(sp_conf,"Figure4","E","Spatial section true-versus-predicted stages")
write_final_source(s8_primary,"Figure4","F","Spatial donor-stage ordinal trajectories")

fig4_final <- (pf4a | pf4b | pf4c) / (pf4d | pf4e | pf4f) +
  patchwork::plot_annotation(tag_levels="A",theme=theme_patch_tags)
fig4_panels <- list(A=pf4a,B=pf4b,C=pf4c,D=pf4d,E=pf4e,F=pf4f)
fig4_dims <- list(A=c(5.8,4.7),B=c(6.2,4.8),C=c(5.8,4.7),D=c(7.2,5.4),E=c(5.8,4.8),F=c(6.0,4.9))
for (nm in names(fig4_panels)) save_final_tif(fig4_panels[[nm]],paste0("Figure4_",nm,".tif"),fig4_dims[[nm]][1],fig4_dims[[nm]][2],FINAL_MAIN_PANEL)
save_final_tif(fig4_final,"Figure4.tif",FINAL_WIDTH,8.9,FINAL_MAIN_COMP)

# ============================================================================
# FIGURE 5 — VENOUS-ULCER DISEASE CONTROL
# ============================================================================
message("Building final Figure 5...")
pf5a <- p5a
write_final_source(s9_primary |> dplyr::select(sample_key,group,z_Skin,z_Wound1,z_Wound7,z_Wound30),"Figure5","A","Venous-ulcer frozen stage scores")
write_final_source(s9_primary |> dplyr::select(sample_key,group,early_state_persistence),"Figure5","B","Venous-ulcer early-state persistence")
pf5b <- scatter_group_plot(
  s9_primary,"group","early_state_persistence","Venous ulcer retains early-state activity","Early-state persistence",
  levels=c("normal_skin","venous_ulcer"),labels=c("Healthy skin","Venous ulcer"),fills=c(normal_skin_col,vu_col)
)
write_final_source(s9_primary |> dplyr::select(sample_key,group,late_remodeling_state),"Figure5","C","Venous-ulcer late-remodeling state")
pf5c <- scatter_group_plot(
  s9_primary,"group","late_remodeling_state","Late remodeling is less dominant in venous ulcer","Late-remodeling state",
  levels=c("normal_skin","venous_ulcer"),labels=c("Healthy skin","Venous ulcer"),fills=c(normal_skin_col,vu_col)
)
pf5d <- p5f
write_final_source(patient_failure_plane,"Figure5","D","Patient-level wound-state failure plane")

balance_core <- core_patient |>
  dplyr::transmute(sample_id=patient_id,disease="Keloid cohort",group=ifelse(group=="keloid","Keloid","Normal scar"),
                   balance=late_remodeling_state-early_state_persistence)
balance_vu <- s9_primary |>
  dplyr::transmute(sample_id=sample_key,disease="Venous-ulcer cohort",group=ifelse(group=="venous_ulcer","Venous ulcer","Healthy skin"),
                   balance=late_remodeling_state-early_state_persistence)
balance_all <- dplyr::bind_rows(balance_core,balance_vu) |>
  dplyr::mutate(group=factor(group,levels=c("Normal scar","Keloid","Healthy skin","Venous ulcer")))
write_final_source(balance_all,"Figure5","E","Patient-level late-minus-early balance across keloid and venous-ulcer cohorts")
pf5e <- ggplot2::ggplot(balance_all,ggplot2::aes(group,balance,fill=group))+
  ggplot2::geom_hline(yintercept=0,linetype="dashed",linewidth=0.4,colour="#888888")+
  ggplot2::geom_point(shape=21,colour="black",stroke=0.5,size=3.0,position=ggplot2::position_jitter(width=0.07,seed=20260908))+
  ggplot2::stat_summary(fun=mean,geom="point",shape=95,size=7.5,colour="black")+
  ggplot2::facet_wrap(~disease,scales="free_x",nrow=1)+
  ggplot2::scale_fill_manual(values=c("Normal scar"=control_col,"Keloid"=keloid_col,"Healthy skin"=normal_skin_col,"Venous ulcer"=vu_col),guide="none")+
  ggplot2::labs(title="Late-minus-early balance distinguishes the two failure patterns",x=NULL,y="Late-remodeling − early-state persistence")+
  final_theme+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=25,hjust=1))

fig5_final <- (pf5a | pf5b | pf5c) / (pf5d | pf5e) +
  patchwork::plot_annotation(tag_levels="A",theme=theme_patch_tags) +
  patchwork::plot_layout(heights=c(1,1.05))
fig5_panels <- list(A=pf5a,B=pf5b,C=pf5c,D=pf5d,E=pf5e)
fig5_dims <- list(A=c(6.0,4.6),B=c(5.4,4.5),C=c(5.4,4.5),D=c(6.4,5.2),E=c(6.0,4.8))
for (nm in names(fig5_panels)) save_final_tif(fig5_panels[[nm]],paste0("Figure5_",nm,".tif"),fig5_dims[[nm]][1],fig5_dims[[nm]][2],FINAL_MAIN_PANEL)
save_final_tif(fig5_final,"Figure5.tif",FINAL_WIDTH,8.1,FINAL_MAIN_COMP)

# ============================================================================
# FIGURE 6 — CROSS-COHORT GENES, GSEA, PATHWAYS, ROBUSTNESS
# ============================================================================
message("Building final Figure 6...")
pf6a <- p6c
write_final_source(concord,"Figure6","A","All-gene cross-cohort effect concordance with frozen-signature classes")

# ---- 6B patient-level expression heatmap when Stage5/6B RDS objects exist ----
canonical_counts_final <- function(counts) {
  if (is.null(rownames(counts))) stop("Count matrix lacks row names")
  genes <- toupper(trimws(as.character(rownames(counts))))
  keep <- !is.na(genes) & nzchar(genes)
  counts <- counts[keep,,drop=FALSE]
  genes <- genes[keep]
  if (inherits(counts,"sparseMatrix")) counts <- as.matrix(counts)
  storage.mode(counts) <- "double"
  rowsum(counts,group=genes,reorder=FALSE)
}

build_patient_gene_heatmap_final <- function() {
  desired <- c("POSTN","COL8A1","COL10A1","COL11A1","ADAM12","ADAMTS14","ASPN","MMP11","MDK","PIEZO2")
  selected <- desired[desired %in% toupper(late_genes$gene)]
  if (length(selected) < 8) selected <- toupper(head(late_genes$gene,12))

  all_rds <- list.files(PROJECT_ROOT,recursive=TRUE,full.names=TRUE,pattern="\\.rds$",ignore.case=TRUE)
  r5 <- all_rds[tolower(basename(all_rds))==tolower("GSE163973_fixed_wound_state_mapping.rds")]
  r6 <- all_rds[tolower(basename(all_rds))==tolower("GSE181316_fixed_wound_state_replication.rds")]
  if (!length(r5) || !length(r6) || !requireNamespace("edgeR",quietly=TRUE)) return(NULL)

  o5 <- readRDS(r5[which.max(file.info(r5)$mtime)])
  o6 <- readRDS(r6[which.max(file.info(r6)$mtime)])
  if (is.null(o5$disease_pseudobulk$counts) || is.null(o6$high_specificity$patient_counts)) return(NULL)

  i5 <- o5$disease_pseudobulk$sample_info
  i5 <- i5[i5$main_cell_type=="Fibroblast",,drop=FALSE]
  c5 <- o5$disease_pseudobulk$counts[,match(i5$sample_key,colnames(o5$disease_pseudobulk$counts)),drop=FALSE]
  c5 <- canonical_counts_final(c5)
  colnames(c5) <- i5$sample_key
  g5 <- tolower(i5$group); g5[g5 %in% c("kl") ] <- "keloid"; g5[g5 %in% c("ns") ] <- "normal_scar"

  i6 <- o6$high_specificity$patient_info
  c6 <- canonical_counts_final(o6$high_specificity$patient_counts)
  id6 <- if ("patient_id" %in% names(i6)) as.character(i6$patient_id) else colnames(c6)
  idx6 <- match(id6,colnames(c6)); if (anyNA(idx6)) return(NULL)
  c6 <- c6[,idx6,drop=FALSE]; colnames(c6) <- id6
  g6 <- tolower(as.character(i6$group)); g6[g6 %in% c("healthy_skin","normal_skin","ns") ] <- "normal_scar"

  norm_one <- function(counts) {
    d <- edgeR::DGEList(counts=counts)
    cpm0 <- edgeR::cpm(d,log=FALSE)
    keep <- rowSums(cpm0>=1)>=2
    d <- d[keep,,keep.lib.sizes=FALSE]
    d <- edgeR::calcNormFactors(d,method="TMM")
    edgeR::cpm(d,log=TRUE,prior.count=2)
  }
  e5 <- norm_one(c5); e6 <- norm_one(c6)
  genes <- selected[selected %in% rownames(e5) & selected %in% rownames(e6)]
  if (length(genes)<6) return(NULL)

  zlong <- function(expr,ids,groups,cohort) {
    m <- expr[genes,,drop=FALSE]
    z <- t(scale(t(m)))
    z[!is.finite(z)] <- 0
    df <- as.data.frame(z); df$gene <- rownames(df)
    long <- tidyr::pivot_longer(df,-gene,names_to="sample_id",values_to="z_expr")
    map <- data.frame(sample_id=ids,group=groups,stringsAsFactors=FALSE)
    dplyr::left_join(long,map,by="sample_id") |> dplyr::mutate(cohort=cohort)
  }
  long <- dplyr::bind_rows(
    zlong(e5,colnames(e5),g5,"GSE163973"),
    zlong(e6,colnames(e6),g6,"GSE181316")
  ) |>
    dplyr::mutate(
      group=factor(group,levels=c("normal_scar","keloid")),
      gene=factor(gene,levels=rev(genes)),
      sample_label=paste(cohort,sample_id,sep=" | ")
    )
  ord <- long |> dplyr::distinct(cohort,sample_id,group,sample_label) |>
    dplyr::arrange(cohort,group,sample_id) |> dplyr::pull(sample_label)
  long$sample_label <- factor(long$sample_label,levels=ord)
  p <- ggplot2::ggplot(long,ggplot2::aes(sample_label,gene,fill=z_expr))+
    ggplot2::geom_tile(colour="white",linewidth=0.25)+
    ggplot2::facet_grid(~cohort,scales="free_x",space="free_x")+
    ggplot2::scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,limits=c(-2.5,2.5),oob=scales::squish,name="Within-cohort\ngene z-score")+
    ggplot2::labs(title="Representative consensus late-remodeling genes",x=NULL,y=NULL)+
    final_theme+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=60,hjust=1,size=6.7),panel.grid=ggplot2::element_blank())
  list(plot=p,source=long)
}

heat6b <- tryCatch(build_patient_gene_heatmap_final(),error=function(e){final_warn("Figure6B patient heatmap fallback:",conditionMessage(e));NULL})
if (!is.null(heat6b)) {
  pf6b <- heat6b$plot
  write_final_source(heat6b$source,"Figure6","B","Within-cohort standardized patient-level expression of representative consensus genes")
} else {
  fallback_genes <- late_genes |> dplyr::arrange(dplyr::desc(minimum_oriented_core_effect)) |> dplyr::slice_head(n=12)
  f6b_src <- fallback_genes |> dplyr::select(gene,mean_difference_GSE163973,mean_difference_GSE181316) |>
    tidyr::pivot_longer(c(mean_difference_GSE163973,mean_difference_GSE181316),names_to="cohort",values_to="mean_difference") |>
    dplyr::mutate(cohort=dplyr::recode(cohort,mean_difference_GSE163973="GSE163973",mean_difference_GSE181316="GSE181316"),gene=factor(gene,levels=rev(unique(fallback_genes$gene))))
  pf6b <- ggplot2::ggplot(f6b_src,ggplot2::aes(cohort,gene,fill=mean_difference))+
    ggplot2::geom_tile(colour="white",linewidth=0.3)+ggplot2::scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,name="Mean\ndifference")+
    ggplot2::labs(title="Consensus late-remodeling gene effects",x=NULL,y=NULL)+final_theme
  write_final_source(f6b_src,"Figure6","B","Fallback cohort-level consensus gene effect heatmap")
}

# ---- weighted running enrichment curves from actual ranked core gene effects ----
make_running_es_final <- function(all_effects, gene_set, set_name) {
  d <- all_effects |>
    dplyr::filter(is.finite(equal_weight_mean_difference)) |>
    dplyr::arrange(dplyr::desc(equal_weight_mean_difference)) |>
    dplyr::mutate(rank=dplyr::row_number(),hit=toupper(gene) %in% toupper(gene_set),weight=abs(equal_weight_mean_difference))
  hit_sum <- sum(d$weight[d$hit]); miss_n <- sum(!d$hit)
  if (hit_sum<=0 || miss_n<=0) stop("Invalid running ES inputs for ",set_name)
  d <- d |> dplyr::mutate(step=ifelse(hit,weight/hit_sum,-1/miss_n),running_ES=cumsum(step),set=set_name)
  d
}

w30_set <- gene_contrib |> dplyr::filter(stage=="Wound30") |> dplyr::pull(gene)
skin_set <- gene_contrib |> dplyr::filter(stage=="Skin") |> dplyr::pull(gene)
es_w30 <- make_running_es_final(all_gene_effects,w30_set,"Wound30 fixed signature")
es_skin <- make_running_es_final(all_gene_effects,skin_set,"Skin fixed signature")
write_final_source(es_w30,"Figure6","C","Weighted running enrichment curve for the frozen Wound30 signature")
write_final_source(es_skin,"Figure6","D","Weighted running enrichment curve for the frozen Skin signature")

pf6c <- ggplot2::ggplot(es_w30,ggplot2::aes(rank,running_ES))+
  ggplot2::geom_hline(yintercept=0,colour="#999999",linewidth=0.35)+ggplot2::geom_line(colour=stage_cols_final["Wound30"],linewidth=0.9)+
  ggplot2::geom_rug(data=dplyr::filter(es_w30,hit),sides="b",alpha=0.45,colour=stage_cols_final["Wound30"])+
  ggplot2::labs(title=paste0("Wound30 signature enrichment (NES ",sprintf("%.2f",stage_gsea$NES[grepl("WOUND30",stage_gsea$ID)]),")"),x="Genes ranked by integrated keloid effect",y="Running enrichment score")+final_theme
pf6d <- ggplot2::ggplot(es_skin,ggplot2::aes(rank,running_ES))+
  ggplot2::geom_hline(yintercept=0,colour="#999999",linewidth=0.35)+ggplot2::geom_line(colour="#2166AC",linewidth=0.9)+
  ggplot2::geom_rug(data=dplyr::filter(es_skin,hit),sides="b",alpha=0.45,colour="#2166AC")+
  ggplot2::labs(title=paste0("Skin signature depletion (NES ",sprintf("%.2f",stage_gsea$NES[grepl("SKIN",stage_gsea$ID)]),")"),x="Genes ranked by integrated keloid effect",y="Running enrichment score")+final_theme

# ---- selected ECM/collagen pathways ----
ora_comb <- dplyr::bind_rows(
  late_gobp_ora |> dplyr::mutate(resource="GO BP"),
  late_reactome_ora |> dplyr::mutate(resource="Reactome")
) |>
  dplyr::filter(is.finite(p.adjust),p.adjust>0) |>
  dplyr::mutate(key=tolower(Description),neglogFDR=-log10(p.adjust))
ora_focus <- ora_comb |> dplyr::filter(grepl("extracellular|matrix|collagen|fibril|matrix organization|ecm",key)) |>
  dplyr::arrange(p.adjust)
if (nrow(ora_focus)<6) ora_focus <- ora_comb |> dplyr::arrange(p.adjust)
ora_focus <- ora_focus |> dplyr::distinct(Description,.keep_all=TRUE) |> dplyr::slice_head(n=10) |>
  dplyr::mutate(Description=stringr::str_wrap(Description,38),Description=factor(Description,levels=rev(unique(Description))))
write_final_source(ora_focus,"Figure6","E","Selected ECM/collagen over-representation results")
pf6e <- ggplot2::ggplot(ora_focus,ggplot2::aes(RichFactor,Description,size=Count,colour=neglogFDR))+
  ggplot2::geom_point(alpha=0.9)+ggplot2::scale_colour_viridis_c(option="C",name="-log10 FDR")+
  ggplot2::scale_size_continuous(range=c(2.5,6.5),name="Genes")+
  ggplot2::labs(title="Persistent late remodeling is enriched for ECM and collagen programs",x="Rich factor",y=NULL)+
  final_theme+ggplot2::theme(legend.position="right",axis.text.y=ggplot2::element_text(size=7.0))

# ---- Stage13 robustness compact panel ----
robust_long <- dplyr::bind_rows(
  stage13_logo |> dplyr::transmute(analysis="Leave-one-gene-out",relative_effect=integrated_effect_fraction_of_baseline),
  stage13_random |> dplyr::transmute(analysis="Random 80% subsampling",relative_effect=integrated_effect_fraction_of_baseline)
)
write_final_source(robust_long,"Figure6","F","Relative integrated effect across Stage13 robustness iterations")
pf6f <- ggplot2::ggplot(robust_long,ggplot2::aes(analysis,relative_effect,fill=analysis))+
  ggplot2::geom_violin(width=0.75,alpha=0.65,colour="black",linewidth=0.35,trim=TRUE)+
  ggplot2::geom_boxplot(width=0.16,outlier.shape=NA,fill="white",colour="black",linewidth=0.4)+
  ggplot2::geom_hline(yintercept=1,linetype="dashed",linewidth=0.45,colour="#555555")+
  ggplot2::geom_hline(yintercept=0.8,linetype="dotted",linewidth=0.45,colour="#999999")+
  ggplot2::scale_fill_manual(values=c("Leave-one-gene-out"="#4C78A8","Random 80% subsampling"="#72B7B2"),guide="none")+
  ggplot2::labs(title="Late-remodeling effect is robust to signature perturbation",x=NULL,y="Integrated effect / full-signature effect")+
  final_theme+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=18,hjust=1))

fig6_final <- (pf6a | pf6b) / (pf6c | pf6d) / (pf6e | pf6f) +
  patchwork::plot_annotation(tag_levels="A",theme=theme_patch_tags)
fig6_panels <- list(A=pf6a,B=pf6b,C=pf6c,D=pf6d,E=pf6e,F=pf6f)
fig6_dims <- list(A=c(6.0,5.2),B=c(6.5,5.0),C=c(6.0,4.5),D=c(6.0,4.5),E=c(6.2,5.1),F=c(5.8,4.6))
for (nm in names(fig6_panels)) save_final_tif(fig6_panels[[nm]],paste0("Figure6_",nm,".tif"),fig6_dims[[nm]][1],fig6_dims[[nm]][2],FINAL_MAIN_PANEL)
save_final_tif(fig6_final,"Figure6.tif",FINAL_WIDTH,10.4,FINAL_MAIN_COMP)

# ============================================================================
# SUPPLEMENTARY FIGURES S1-S8: COPY VALIDATED LEGACY PANELS/COMPOSITES
# ============================================================================
message("Copying validated supplementary figures S1-S8...")
legacy_panel_dir <- file.path(PROJECT_ROOT,"__TEMP_KELOID_PANEL_ENGINE__","individual_panel_tif")
legacy_comp_dir <- file.path(PROJECT_ROOT,"__TEMP_KELOID_COMPOSITE_ENGINE__","supplementary_tif")
if (!dir.exists(legacy_panel_dir)) legacy_panel_dir <- PANEL_DIR

for (s in 1:8) {
  patt <- paste0("^Supplementary_Figure_S",s,"_.*\\.tif$")
  pfiles <- list.files(legacy_panel_dir,pattern=patt,full.names=TRUE,ignore.case=TRUE)
  if (length(pfiles)) file.copy(pfiles,FINAL_SUPP_PANEL,overwrite=TRUE)
  cfile <- file.path(legacy_comp_dir,paste0("Supplementary_Figure_S",s,".tif"))
  if (file.exists(cfile)) {
    file.copy(cfile,file.path(FINAL_SUPP_COMP,basename(cfile)),overwrite=TRUE)
  } else {
    final_warn("Legacy supplementary composite missing:",cfile)
  }
}

# ============================================================================
# SUPPLEMENTARY FIGURE S9 — FULL STAGE13 ROBUSTNESS
# ============================================================================
message("Building Supplementary Figure S9...")
s9a_src <- stage13_logo |> dplyr::arrange(integrated_effect_fraction_of_baseline) |>
  dplyr::mutate(omitted_gene=factor(omitted_gene,levels=omitted_gene))
write_final_source(s9a_src,"Supplementary_Figure_S9","A","All leave-one-Wound30-gene-out integrated effects")
ps9a <- ggplot2::ggplot(s9a_src,ggplot2::aes(omitted_gene,integrated_mean_difference))+
  ggplot2::geom_hline(yintercept=stage13_baseline$integrated_mean_difference[1],linetype="dashed",colour="#555555")+
  ggplot2::geom_point(size=1.8,colour="#4C78A8")+
  ggplot2::labs(title="Leave-one-Wound30-gene-out stability",x="Omitted gene",y="Integrated mean difference")+
  final_theme+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90,hjust=1,size=5.8))

s9b_src <- stage13_influence |> dplyr::arrange(dplyr::desc(absolute_loss_fraction)) |> dplyr::slice_head(n=15) |>
  dplyr::mutate(omitted_gene=factor(omitted_gene,levels=rev(omitted_gene)))
write_final_source(s9b_src,"Supplementary_Figure_S9","B","Top single-gene influence ranking")
ps9b <- ggplot2::ggplot(s9b_src,ggplot2::aes(absolute_loss_fraction,omitted_gene))+
  ggplot2::geom_col(fill="#4C78A8",colour="black",linewidth=0.25)+
  ggplot2::scale_x_continuous(labels=scales::percent_format(accuracy=0.1))+
  ggplot2::labs(title="Largest effect loss after removing one Wound30 gene",x="Fractional loss of integrated effect",y=NULL)+final_theme

write_final_source(stage13_random,"Supplementary_Figure_S9","C","All 1000 random 80-percent Wound30 subsampling results")
ps9c <- ggplot2::ggplot(stage13_random,ggplot2::aes(integrated_mean_difference))+
  ggplot2::geom_histogram(bins=35,fill="#72B7B2",colour="white",linewidth=0.25)+
  ggplot2::geom_vline(xintercept=stage13_baseline$integrated_mean_difference[1],linetype="dashed",linewidth=0.6,colour="black")+
  ggplot2::labs(title="Random 80% Wound30 subsampling",x="Integrated mean difference",y="Iterations")+final_theme

write_final_source(stage13_random |> dplyr::select(iteration_id,integrated_effect_fraction_of_baseline,both_cohorts_expected_direction,both_cohorts_complete_separation,integrated_exact_two_sided_p),
                   "Supplementary_Figure_S9","D","Relative effects and directional stability in random subsampling")
ps9d <- ggplot2::ggplot(stage13_random,ggplot2::aes(y=integrated_effect_fraction_of_baseline,x="Random 80% subsampling"))+
  ggplot2::geom_violin(fill="#72B7B2",alpha=0.7,colour="black")+
  ggplot2::geom_boxplot(width=0.14,outlier.shape=NA,fill="white")+
  ggplot2::geom_hline(yintercept=c(0.8,1),linetype=c("dotted","dashed"),colour=c("#888888","black"))+
  ggplot2::labs(title="Relative effect retention",x=NULL,y="Integrated effect / baseline")+final_theme

supp9 <- (ps9a | ps9b) / (ps9c | ps9d) + patchwork::plot_annotation(tag_levels="A",theme=theme_patch_tags)
s9_panels <- list(A=ps9a,B=ps9b,C=ps9c,D=ps9d)
for (nm in names(s9_panels)) save_final_tif(s9_panels[[nm]],paste0("Supplementary_Figure_S9_",nm,".tif"),6.0,4.8,FINAL_SUPP_PANEL)
save_final_tif(supp9,"Supplementary_Figure_S9.tif",FINAL_WIDTH,8.2,FINAL_SUPP_COMP)

# ============================================================================
# TABLES + RAW/SOURCE DATA EXPORT
# ============================================================================
message("Exporting tables and raw numeric inputs...")

# Copy all legacy manuscript tables if available.
legacy_tables <- file.path(PROJECT_ROOT,"__TEMP_KELOID_PANEL_ENGINE__","tables")
if (dir.exists(legacy_tables)) {
  tabfiles <- list.files(legacy_tables,full.names=TRUE,include.dirs=FALSE)
  if (length(tabfiles)) file.copy(tabfiles,FINAL_TABLES,overwrite=TRUE)
}

# Explicit final Table 1 and Table 2 to match the new narrative.
table1_final <- data.frame(
  Dataset=c("GSE241132","GSE163973","GSE181316","GSE220300","GSE241124","GSE265972","GSE181297"),
  Data_type=c("scRNA-seq","scRNA-seq","scRNA-seq","Lesion-context transcriptomics","Spatial transcriptomics","Transcriptomics","scRNA-seq + Visium"),
  Biological_unit=c("Donor-stage pseudobulk","Patient-level fibroblast pseudobulk","Patient-level fibroblast pseudobulk","Lesion/scar unit","Section/donor","Patient-level fibroblast pseudobulk","Small scRNA + section-level"),
  Main_structure=c("3 donors × 4 stages","3 keloids vs 3 normal scars","3 keloids vs 3 normal scars","External lesion/scar units","16 sections from 4 donors","4 ulcers vs 5 healthy skin","Exploratory mixed data"),
  Analytical_role=c("Reference construction + donor-held-out validation","Discovery","Independent replication","External support","Moderate cross-modal spatial support","Disease-control comparison","Exploratory / mixed evidence"),
  stringsAsFactors=FALSE
)
utils::write.csv(table1_final,file.path(FINAL_TABLES,"Table1_Public_datasets_and_analytical_roles.csv"),row.names=FALSE,fileEncoding="UTF-8")

table2_final <- core_integrated |>
  dplyr::left_join(
    core_effects |> dplyr::select(cohort,metric,mean_difference) |> tidyr::pivot_wider(names_from=cohort,values_from=mean_difference,names_prefix="effect_"),
    by="metric"
  )
utils::write.csv(table2_final,file.path(FINAL_TABLES,"Table2_Core_patient_level_wound_state_effects.csv"),row.names=FALSE,fileEncoding="UTF-8")
utils::write.csv(stage13_summary,file.path(FINAL_TABLES,"Supplementary_Table_S11_Wound30_signature_robustness_summary.csv"),row.names=FALSE,fileEncoding="UTF-8")
utils::write.csv(stage13_influence,file.path(FINAL_TABLES,"Supplementary_Table_S11B_single_gene_influence.csv"),row.names=FALSE,fileEncoding="UTF-8")

# Copy every numeric CSV resolved by the validated legacy input manifest.
if (exists("input_manifest") && nrow(input_manifest)) {
  raw_manifest <- input_manifest |> dplyr::filter(status=="found",nzchar(resolved_path),file.exists(resolved_path))
  for (i in seq_len(nrow(raw_manifest))) {
    src <- raw_manifest$resolved_path[i]
    dest <- file.path(FINAL_RAW,paste0(sprintf("%03d",i),"_",final_safe_name(raw_manifest$logical_name[i]),"__",basename(src)))
    file.copy(src,dest,overwrite=TRUE)
  }
  utils::write.csv(raw_manifest,file.path(FINAL_LOGS,"resolved_numeric_input_manifest.csv"),row.names=FALSE,fileEncoding="UTF-8")
}
# Stage13 numeric inputs
for (p in c(stage13_logo_path,stage13_random_path,stage13_summary_path,stage13_influence_path,stage13_baseline_path)) {
  file.copy(p,file.path(FINAL_RAW,paste0("Stage13__",basename(p))),overwrite=TRUE)
}

utils::write.csv(FINAL_SOURCE_INDEX,file.path(FINAL_SOURCE,"Source_Data_Index.csv"),row.names=FALSE,fileEncoding="UTF-8")

# Final file manifest
all_final_files <- list.files(FINAL_DIR,recursive=TRUE,full.names=TRUE,include.dirs=FALSE)
final_root_norm <- normalizePath(FINAL_DIR, winslash="/", mustWork=TRUE)
final_files_norm <- normalizePath(all_final_files, winslash="/", mustWork=TRUE)
relative_final <- ifelse(
  startsWith(final_files_norm, paste0(final_root_norm, "/")),
  substring(final_files_norm, nchar(final_root_norm) + 2L),
  basename(final_files_norm)
)
final_manifest <- data.frame(
  relative_path=relative_final,
  extension=tolower(tools::file_ext(all_final_files)),
  size_bytes=file.info(all_final_files)$size,
  stringsAsFactors=FALSE
)
utils::write.csv(final_manifest,file.path(FINAL_LOGS,"FINAL_OUTPUT_MANIFEST.csv"),row.names=FALSE,fileEncoding="UTF-8")

capture.output(sessionInfo(),file=file.path(FINAL_LOGS,"SESSION_INFO.txt"))
if (length(FINAL_WARNINGS)) writeLines(unique(FINAL_WARNINGS),file.path(FINAL_LOGS,"WARNINGS.txt"),useBytes=TRUE)

# Verify figure outputs are TIFF only.
figure_dirs <- c(FINAL_MAIN_PANEL,FINAL_MAIN_COMP,FINAL_SUPP_PANEL,FINAL_SUPP_COMP)
non_tiff <- unlist(lapply(figure_dirs,function(d) list.files(d,full.names=TRUE,include.dirs=FALSE)))
non_tiff <- non_tiff[tolower(tools::file_ext(non_tiff)) != "tif"]
if (length(non_tiff)) stop("Non-TIFF figure outputs detected: ",paste(non_tiff,collapse=" | "))

required_main <- c(paste0("Figure",1:6,".tif"))
missing_main <- required_main[!file.exists(file.path(FINAL_MAIN_COMP,required_main))]
required_supp <- paste0("Supplementary_Figure_S",1:9,".tif")
missing_supp <- required_supp[!file.exists(file.path(FINAL_SUPP_COMP,required_supp))]
if (length(missing_main) || length(missing_supp)) {
  stop("Final composite verification failed. Missing: ",paste(c(missing_main,missing_supp),collapse=", "))
}

writeLines(c(
  "FINAL_STORY_LAYOUT_V3_COMPLETED",
  paste0("time=",format(Sys.time(),"%Y-%m-%d %H:%M:%S")),
  paste0("project_root=",PROJECT_ROOT),
  paste0("final_output=",FINAL_DIR),
  "figure_format=TIFF only",
  paste0("dpi=",FINAL_DPI),
  "compression=LZW",
  "main_figures=6",
  "supplementary_figures=9",
  "source_data=panel-level CSV + copies of raw numeric inputs used",
  "Figure1=graphical summary/workflow; Figure2=normal wound reference; Figure3=discovery+replication; Figure4=external+spatial; Figure5=venous-ulcer control; Figure6=genes+pathways+robustness"
),file.path(FINAL_LOGS,"FINAL_COMPLETED.txt"),useBytes=TRUE)

# Optional cleanup of intermediate legacy figure-engine folders.
CLEAN_TEMP_ENGINES <- TRUE
if (isTRUE(CLEAN_TEMP_ENGINES)) {
  for (d in c(file.path(PROJECT_ROOT,"__TEMP_KELOID_PANEL_ENGINE__"),file.path(PROJECT_ROOT,"__TEMP_KELOID_COMPOSITE_ENGINE__"))) {
    if (dir.exists(d)) unlink(d,recursive=TRUE,force=TRUE)
  }
}

message("============================================================")
message("FINAL STORY-LAYOUT V3 COMPLETED")
message("Output: ",FINAL_DIR)
message("All figure files are TIFF only, 600 dpi, LZW.")
message("Panel source data and raw numeric inputs are exported as CSV.")
message("============================================================")
