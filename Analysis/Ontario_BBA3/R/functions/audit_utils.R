# === audit_utils.R ===
# Run from your project root directory
setwd("./R/")

rm(list=ls())

# --- Config ---
utils_files <- c(
  "functions/covariate_processing_utils.R",
  "functions/figure_utils.R",
  "functions/inla_model_utils2.R",
  "functions/cv_utils.R",
  "functions/spatial_utils.R",
  "functions/survey_processing_utils.R"
)


workflow_files <- paste0(list.files(pattern = "^[0-9]+_.*\\.R$"))

# --- Extract function definitions from utils ---
extract_defined_functions <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  # Matches: my_func <- function(...) or my_func = function(...)
  matches <- regmatches(lines, regexpr("^([a-zA-Z0-9_.]+)\\s*(<-|=)\\s*function", lines))
  fns <- sub("\\s*(<-|=)\\s*function", "", matches)
  return(fns)
}

# --- Search for calls in workflow scripts ---
find_calls <- function(fn_names, search_files) {
  all_code <- unlist(lapply(search_files, readLines, warn = FALSE))
  
  sapply(fn_names, function(fn) {
    # Look for fn( in the workflow code
    any(grepl(paste0("\\b", fn, "\\s*\\("), all_code))
  })
}

# --- Run audit ---
results <- lapply(utils_files, function(uf) {
  fns <- extract_defined_functions(uf)
  if (length(fns) == 0) return(NULL)
  called <- find_calls(fns, workflow_files)
  data.frame(
    utils_file = uf,
    function_name = fns,
    used_in_workflow = called,
    stringsAsFactors = FALSE
  )
})

audit <- do.call(rbind, results)

# --- Report ---
cat("\n=== UNUSED FUNCTIONS (safe to review for deletion) ===\n")
unused <- audit[!audit$used_in_workflow, ]
print(unused[, c("utils_file", "function_name")], row.names = FALSE)

cat("\n=== USED FUNCTIONS (keep) ===\n")
used <- audit[audit$used_in_workflow, ]
print(used[, c("utils_file", "function_name")], row.names = FALSE)

# Optional: save full report
write.csv(audit, "utils_audit.csv", row.names = FALSE)
cat("\nFull report saved to utils_audit.csv\n")