# Download count matrices (SingleCellExperiment / .rds) for all 6
# pediatric high-grade glioma samples used in this analysis
#
# SCPCS000002, SCPCS000014, SCPCS000019, SCPCS000020, SCPCS000016, SCPCS000022.

library(ScPCAr)
source("R/00_config.R")

my_email <- Sys.getenv("SCPCA_EMAIL")   # set before running: export SCPCA_EMAIL="you@example.org"
auth_token <- get_auth(email = my_email, agree = TRUE)

# ---- patients ID -------------------------------------------------
hgg6_samples <- c(  # = sample_ids in 00_config.R
  "SCPCS000002",  
  "SCPCS000014", 
  "SCPCS000019",
  "SCPCS000020",
  "SCPCS000016",
  "SCPCS000022" 
)

ds_id <- create_dataset(
  samples    = hgg6_samples,
  format     = "sce",
  email      = my_email,
  auth_token = auth_token
)
cat("Dataset ID:", ds_id, "\n")

# --- Dataset processing -----------------------------
start_dataset_processing(ds_id, auth_token = auth_token)

# --- Download ----------------
hgg6_file_paths <- download_dataset(
  ds_id,
  destination      = "scpca_data_hgg6",
  await_processing = TRUE,
  auth_token       = auth_token
)
