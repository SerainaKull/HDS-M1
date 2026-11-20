# Final-dataset checks for the supplied LACDHS project.
# Source this file before calling validate_ard(). No validate package required.
# Thresholds below encode YOUR existing specification; they are not clinical
# recommendations. patient_seq_id identifies a source record, not a person.

.ard_results <- function() {
  data.frame(Check = character(), Requirement = character(),
             Evaluated = integer(), Failed = integer(),
             Status = character(), Detail = character(),
             stringsAsFactors = FALSE)
}

.ard_result <- function(id, description, result) {
  error_text <- ""
  result <- tryCatch(force(result), error = function(e) {
    error_text <<- conditionMessage(e)
    FALSE
  })
  valid <- is.logical(result) && length(result) > 0L
  failed <- if (valid) sum(is.na(result) | !result) else 1L
  data.frame(Check = id, Requirement = description,
             Evaluated = if (valid) length(result) else 0L,
             Failed = failed,
             Status = if (failed == 0L) "PASS" else "FAIL",
             Detail = error_text, stringsAsFactors = FALSE)
}

.ard_gate <- function(results) {
  failed <- results[results$Status != "PASS", , drop = FALSE]
  if (nrow(failed)) {
    print(failed, row.names = FALSE)
    stop("ARD validation failed; see the checks printed above.", call. = FALSE)
  }
  invisible(results)
}

# Elementwise equality, with missing values equal only to missing values.
.ard_equal <- function(x, y) {
  if (length(x) != length(y)) return(FALSE)
  x <- if (is.factor(x)) as.character(x) else x
  y <- if (is.factor(y)) as.character(y) else y
  out <- is.na(x) & is.na(y)
  observed <- !is.na(x) & !is.na(y)
  out[observed] <- x[observed] == y[observed]
  out
}

validate_ard <- function(raw_orig, ard_df, expected_files) {
  results <- .ard_results()
  add <- function(id, description, result) {
    results <<- rbind(results, .ard_result(id, description, result))
  }

  # The 29-column output contract, including column order and storage types.
  schema <- c(
    patient_seq_id = "integer", cohort_split = "factor",
    calendar_period = "factor", dr_status = "integer", dr_factor = "factor",
    age_group = "ordered", sex_code = "factor", ethnicity_code = "factor",
    duration_yrs_clean = "double", hba1c_clean = "double",
    glycemic_control_cat = "factor", bun_clean = "double",
    renal_azotemia_flag = "logical", sbp_clean = "double",
    dbp_clean = "double", stage2_htn_flag = "logical",
    hemoglobin_clean = "double", triglycerides_clean = "double",
    insulin_dep_code = "factor", nephropathy_code = "factor",
    neuropathy_code = "factor", stroke_code = "factor",
    miss_hba1c = "logical", miss_bun = "logical", miss_bp = "logical",
    miss_hemoglobin = "logical", miss_triglycerides = "logical",
    miss_duration = "logical", lab_missing_count = "integer"
  )
  clean_spec <- data.frame(
    source = c("DURATION_OF_DIABETES_YRS", "HEMOGLOBIN_A1C",
               "BLOOD_UREA_NITROGEN", "SYSTOLIC_BLOOD_PRESSURE",
               "DIASTOLIC_BLOOD_PRESSURE", "HEMOGLOBIN", "TRIGLYCERIDES"),
    output = c("duration_yrs_clean", "hba1c_clean", "bun_clean",
               "sbp_clean", "dbp_clean", "hemoglobin_clean",
               "triglycerides_clean"),
    lower = c(0, 3, 0, 50, 30, 3, 0),
    upper = c(85, 25, 250, 300, 200, 25, 5000),
    include_lower = c(TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE),
    stringsAsFactors = FALSE
  )
  code_sources <- c(sex_code = "SEX", ethnicity_code = "ETHNICITY",
                    insulin_dep_code = "INSULIN_DEPENDENCE",
                    nephropathy_code = "NEPHROPATHY",
                    neuropathy_code = "NEUROPATHY", stroke_code = "STROKE")
  required_source <- c("source_file", "cohort_split", "calendar_period",
                       "RETINOPATHY", "AGE", clean_spec$source,
                       unname(code_sources))
  add("STRUCT01", "Expected source and output columns exist",
      is.data.frame(raw_orig) && is.data.frame(ard_df) &&
        all(required_source %in% names(raw_orig)) &&
        identical(names(ard_df), names(schema)))
  .ard_gate(results)  # Later checks require the expected columns.

  type_ok <- vapply(names(schema), function(nm) {
    x <- ard_df[[nm]]
    switch(schema[[nm]],
           factor = identical(class(x), "factor"),
           ordered = identical(class(x), c("ordered", "factor")),
           identical(typeof(x), schema[[nm]]) && !is.object(x))
  }, logical(1))
  add("STRUCT02", "Each output column has its specified type", type_ok)
  add("ROWS01", "All 49,931 source records are retained",
      nrow(raw_orig) == 49931L && nrow(ard_df) == nrow(raw_orig))
  add("ROWS02", "Record IDs cover source rows exactly once",
      !anyNA(ard_df$patient_seq_id) &&
        !anyDuplicated(ard_df$patient_seq_id) &&
        identical(sort(ard_df$patient_seq_id), seq_len(nrow(raw_orig))))
  add("SOURCE01", "Source measurement columns are numeric",
      vapply(raw_orig[clean_spec$source], is.numeric, logical(1)))
  .ard_gate(results)

  # Current project contract: patient_seq_id was assigned by row_number()
  # directly on raw_orig. Align by that ID so ARD reordering is permitted.
  src <- raw_orig[match(ard_df$patient_seq_id, seq_len(nrow(raw_orig))), ]
  split_levels <- c("train_2015_17", "test_2015_17", "ext_val_2018")
  expected_counts <- c(27223L, 13408L, 9300L)
  period_levels <- c("2015-2017", "2018")
  file_keys <- c("train_orig", "test_orig", "ext_orig")
  add("SOURCE02", "Three distinct original source filenames are specified",
      all(file_keys %in% names(expected_files)) &&
        !anyNA(expected_files[file_keys]) &&
        length(unique(unname(expected_files[file_keys]))) == 3L)
  .ard_gate(results)
  source_index <- match(src$source_file, unname(expected_files[file_keys]))
  add("COHORT01", "Each record has the correct original source and split",
      !is.na(source_index) &
        .ard_equal(src$cohort_split, split_levels[source_index]) &
        .ard_equal(ard_df$cohort_split, split_levels[source_index]))
  add("COHORT02", "Each record has the correct calendar period",
      !is.na(source_index) &
        .ard_equal(src$calendar_period,
                   c("2015-2017", "2015-2017", "2018")[source_index]) &
        .ard_equal(ard_df$calendar_period,
                   c("2015-2017", "2015-2017", "2018")[source_index]))
  add("COHORT03", "All three split counts match the specification",
      vapply(split_levels, function(s) sum(ard_df$cohort_split == s),
             integer(1)) == expected_counts)

  outcome <- toupper(trimws(as.character(src$RETINOPATHY)))
  expected_outcome <- match(outcome, c("N", "Y")) - 1L
  add("OUTCOME01", "Source outcomes are complete and recognised",
      !is.na(outcome) & outcome %in% c("N", "Y"))
  add("OUTCOME02", "Outcome integers and labels preserve source meaning",
      .ard_equal(ard_df$dr_status, expected_outcome) &
        .ard_equal(ard_df$dr_factor,
                   c("No_DR", "DR_Present")[expected_outcome + 1L]))

  age_levels <- c("18-29", "30-39", "40-49", "50-59", "60-69", "70-79", "80+")
  age <- trimws(as.character(src$AGE))
  add("AGE01", "Age groups are complete, recognised and preserved",
      !is.na(age) & age %in% age_levels & .ard_equal(ard_df$age_group, age))
  factor_levels <- list(cohort_split = split_levels,
                        calendar_period = period_levels,
                        dr_factor = c("No_DR", "DR_Present"),
                        age_group = age_levels,
                        glycemic_control_cat = c("Controlled_lt7",
                                                 "Suboptimal_7_to_8.9",
                                                 "Uncontrolled_ge9"))
  for (nm in names(code_sources)) {
    domain <- if (nm %in% c("sex_code", "ethnicity_code")) c("0", "1") else c("1", "2")
    original <- as.character(src[[code_sources[[nm]]]])
    target <- ifelse(is.na(original), NA_character_, paste0("Code_", original))
    add(paste0("CODE_", nm), paste("Recognised source codes preserved in", nm),
        (is.na(original) | original %in% domain) & .ard_equal(ard_df[[nm]], target))
    factor_levels[[nm]] <- paste0("Code_", domain)
  }
  add("FACTORS01", "All factor levels and their order match the specification",
      vapply(names(factor_levels), function(nm) {
        identical(levels(ard_df[[nm]]), factor_levels[[nm]])
      }, logical(1)))

  # Independently calculate expected values from the source and specification.
  # DBP also depends on the expected cleaned SBP, not the ARD's SBP.
  expected_clean <- list()
  audit_rows <- list()
  for (i in seq_len(nrow(clean_spec))) {
    spec <- clean_spec[i, ]
    x <- as.numeric(src[[spec$source]])
    keep <- !is.na(x) & is.finite(x) & x <= spec$upper &
      (if (spec$include_lower) x >= spec$lower else x > spec$lower)
    if (spec$output == "dbp_clean") {
      sbp <- expected_clean$sbp_clean
      keep <- keep & (is.na(sbp) | x < sbp)
    }
    expected <- rep(NA_real_, length(x))
    expected[which(keep)] <- x[which(keep)]
    expected_clean[[spec$output]] <- expected
    add(paste0("CLEAN_", spec$output),
        paste("Valid source values preserved and specified exclusions applied:", spec$output),
        .ard_equal(ard_df[[spec$output]], expected))
    audit_rows[[i]] <- data.frame(
      Variable = spec$output, Source_missing = sum(is.na(x)),
      Newly_set_to_NA = sum(!is.na(x) & is.na(expected)),
      Final_missing = sum(is.na(ard_df[[spec$output]])))
  }
  e <- expected_clean
  gly <- rep(NA_character_, nrow(ard_df))
  gly[which(e$hba1c_clean < 7)] <- "Controlled_lt7"
  gly[which(e$hba1c_clean >= 7 & e$hba1c_clean < 9)] <- "Suboptimal_7_to_8.9"
  gly[which(e$hba1c_clean >= 9)] <- "Uncontrolled_ge9"
  htn <- ifelse(is.na(e$sbp_clean) | is.na(e$dbp_clean), NA,
                e$sbp_clean >= 140 | e$dbp_clean >= 90)
  add("PHENO01", "Glycaemic category matches the source-based specification",
      .ard_equal(ard_df$glycemic_control_cat, gly))
  add("PHENO02", "Hypertension flag and its missingness match the specification",
      .ard_equal(ard_df$stage2_htn_flag, htn))
  add("PHENO03", "Renal flag and its missingness match the specification",
      .ard_equal(ard_df$renal_azotemia_flag, e$bun_clean > 20))

  expected_flags <- data.frame(
    miss_hba1c = is.na(e$hba1c_clean), miss_bun = is.na(e$bun_clean),
    miss_bp = is.na(e$sbp_clean) | is.na(e$dbp_clean),
    miss_hemoglobin = is.na(e$hemoglobin_clean),
    miss_triglycerides = is.na(e$triglycerides_clean),
    miss_duration = is.na(e$duration_yrs_clean))
  for (nm in names(expected_flags)) {
    add(paste0("MISSING_", nm), paste("Correct missingness flag:", nm),
        .ard_equal(ard_df[[nm]], expected_flags[[nm]]))
  }
  add("MISSING_TOTAL", "Missingness count equals the six correct flags",
      .ard_equal(ard_df$lab_missing_count, as.integer(rowSums(expected_flags))))
  .ard_gate(results)

  # Sidecar provenance, retaining the 29-column ARD contract.
  # source_row counts data records after CSV import (not physical text lines).
  source_rows <- ave(seq_len(nrow(raw_orig)), raw_orig$source_file, FUN = seq_along)
  lineage <- data.frame(patient_seq_id = ard_df$patient_seq_id,
                        source_file = src$source_file,
                        source_row = as.integer(source_rows[ard_df$patient_seq_id]))
  list(summary = results, cleaning_audit = do.call(rbind, audit_rows),
       lineage = lineage)
}

verify_ard_exports <- function(ard_df, csv_path, rds_path) {
  results <- .ard_results()
  add <- function(id, description, result) {
    results <<- rbind(results, .ard_result(id, description, result))
  }
  add("EXPORT_RDS", "RDS restores the exact validated R object",
      identical(readRDS(rds_path), ard_df))

  # CSV has no factor-level metadata. Read with explicit types and compare
  # values, allowing only a small floating-point round-trip tolerance.
  collectors <- lapply(ard_df, function(x) {
    if (is.factor(x)) return(readr::col_character())
    switch(typeof(x), integer = readr::col_integer(),
           double = readr::col_double(), logical = readr::col_logical(),
           stop("Unexpected output column type."))
  })
  csv <- readr::read_csv(csv_path, col_types = do.call(readr::cols, collectors),
                         na = "NA", show_col_types = FALSE, progress = FALSE)
  add("EXPORT_CSV_STRUCTURE", "CSV retains row count and column names/order",
      nrow(csv) == nrow(ard_df) && identical(names(csv), names(ard_df)))
  add("EXPORT_CSV_PARSE", "CSV has no parsing problems",
      nrow(readr::problems(csv)) == 0L)
  .ard_gate(results)
  add("EXPORT_CSV_VALUES", "CSV retains every column's values and missingness",
      vapply(names(ard_df), function(nm) {
        expected <- ard_df[[nm]]
        if (is.factor(expected)) expected <- as.character(expected)
        observed <- csv[[nm]]
        identical(is.na(observed), is.na(expected)) &&
          isTRUE(all.equal(observed, expected, tolerance = 1e-12,
                           check.attributes = FALSE))
      }, logical(1)))
  .ard_gate(results)
  results
}
