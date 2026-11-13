# ============================================================
# Replicate-Split Validation - Sparse-GY Training Design
# Partially Observed Candidates: GS model fits
# ============================================================
#
# Purpose:
# Refits the multi-trait GS models for partially observed
# candidates so that predictive ability can be estimated
# without shared plot error between the secondary trait and
# the validation GY (Runcie & Cheng, 2019).
#
# Each candidate plot is a "predictor" or "validation" plot,
# assigned from the phenotype data (see make_roles() below;
# no separate role file is needed). As in the original
# scenario, observed GY is masked on every candidate plot and
# in the sparse-GY zones. In addition, the secondary trait is
# masked on the candidates' validation plots, so the model
# only sees secondary-trait data from plots that are not used
# for validation.
#
# For a candidate with two plots in an environment (always in
# different blocks):
# Split A: the plot in the higher block is the validation plot.
# Split B: the plot in the lower block is the validation plot.
# Unreplicated candidate plots are predictor plots in both.
#
# Everything else matches the original
# *_SparseGY_PartiallyObserved.R scripts.
#
# Requires an ASReml-R installation that fits vm() models.
#
# ============================================================

# 1. Load packages -------------------------------------------

library(tidyverse)
library(asreml)
library(here)

# 2. Load supporting functions -------------------------------

source(here("R", "Functions_earlyGS.R"))

# 3. Settings ------------------------------------------------

# One row per secondary trait: data file and FA rank, as in
# the original scripts.
configs <- tribble(
  ~sec,     ~file,                          ~k,
  "IGYMLR", "Pheno_SparseGY_Observed.rds",  4,
  "IGYMF",  "Pheno_SparseGY_Observed.rds",  5,
  "NDRE",   "Pheno_ndre.rds",               5,
  "PC1",    "Pheno_PC1.rds",                5
)

# Runs to fit. Remove rows to fit a subset.
runs <- expand_grid(
  sec = c("IGYMLR", "NDRE", "PC1", "IGYMF"),
  split = c("A", "B")
)

# 4. Load shared inputs --------------------------------------

Ginv <- readRDS(here("Data", "Ginv.rds"))

val <- read.csv(here("Data", "selection_candidates.csv"))

candidates <- unique(as.character(val$Gen))

# Assign each candidate plot a role, from the phenotype data
# itself. Within an environment, a candidate with two plots
# has them in different blocks: ordered by block, the first
# is plot 1 and the second is plot 2.
#   Split A: plot 1 = predictor,  plot 2 = validation
#   Split B: plot 1 = validation, plot 2 = predictor
# A candidate with one plot in an environment has a predictor
# plot only. The rule depends only on Env, Gen, Block, Col
# and Row, so every data file gives the same assignment.
make_roles <- function(d, candidate_names) {
  out <- d |>
    filter(Trait == "GY", as.character(Gen) %in% candidate_names) |>
    transmute(
      plot_key,
      Env = as.character(Env),
      Gen = as.character(Gen),
      b = as.integer(as.character(Block)),
      c = as.integer(as.character(Col)),
      r = as.integer(as.character(Row))
    ) |>
    arrange(Env, Gen, b, c, r) |>
    group_by(Env, Gen) |>
    mutate(n_plots_in_env = n(), idx = row_number()) |>
    ungroup() |>
    mutate(
      role_splitA = if_else(
        n_plots_in_env == 2 & idx == 2, "validation", "predictor"
      ),
      role_splitB = if_else(
        n_plots_in_env == 2 & idx == 1, "validation", "predictor"
      )
    )

  stopifnot(
    !anyDuplicated(out$plot_key),
    max(out$n_plots_in_env) <= 2
  )
  out
}

dir.create(here("Results"), showWarnings = FALSE, recursive = TRUE)

# 5. Fit each run --------------------------------------------

for (i in seq_len(nrow(runs))) {

  sec <- runs$sec[i]
  split <- runs$split[i]
  cfg <- configs[configs$sec == sec, ]
  k <- cfg$k  # Factor-analytic rank

  run_name <- paste0(
    sec, "_SparseGY_RepSplit", split, "_FA", k
  )
  message("\n===== ", run_name, " =====")

  # 5a. Load and prepare phenotypic data ---------------------

  Pheno <- readRDS(here("Data", cfg$file)) |>
    ungroup()

  stopifnot(all(c("Env", "Row", "Col") %in% names(Pheno)))

  Pheno <- Pheno |>
    mutate(
      plot_key = paste(Env, Row, Col, sep = "_"),
      # Block is stored as an integer in Pheno_ndre.rds; it
      # must be a factor to be fitted as a block effect.
      Block = factor(Block)
    )

  roles <- make_roles(Pheno, candidates)
  message(
    "Candidate plots: ", nrow(roles),
    " | validation plots per split: ",
    sum(roles$role_splitA == "validation")
  )

  # Pheno_SparseGY_Observed.rds holds both IGY traits; keep
  # GY plus the one used in this run.
  # IGYMLR: drop the "mf" trait-environments (as in the
  # original script).
  # IGYMF: keep only GY and the "mf" trait-environments.
  te <- as.character(Pheno$TraitEnv)
  if (sec == "IGYMLR") {
    Pheno <- Pheno[!str_starts(te, "mf"), ]
  } else if (sec == "IGYMF") {
    Pheno <- Pheno[str_starts(te, "GY") | str_starts(te, "mf"), ]
  }

  Pheno_f <- Pheno |>
    mutate(
      # Genotypes with genomic relationship information
      Gkeep = factor(
        if_else(
          as.character(Gen) %in% colnames(Ginv),
          as.character(Gen),
          NA_character_
        )
      ),

      # Genotypes without genomic relationship information
      Gdrop = factor(
        if_else(
          as.character(Gen) %in% colnames(Ginv),
          NA_character_,
          as.character(Gen)
        )
      )
    ) |>
    droplevels()

  message("TraitEnv levels in this run:")
  print(levels(factor(Pheno_f$TraitEnv)))

  # 5b. Mask phenotypic observations -------------------------

  role_col <- paste0("role_split", split)
  validation_plots <- roles$plot_key[roles[[role_col]] == "validation"]

  df_seed <- Pheno_f |>
    mutate(
      is_val_plot = plot_key %in% validation_plots,

      # Mask observed GY for selection candidates
      Pheno = if_else(
        Trait == "GY" & Gen %in% val$Gen,
        NA_real_,
        Pheno
      ),

      # Mask observed GY in sparse-GY training zones
      Pheno = if_else(
        Trait == "GY" & mask_zone == TRUE,
        NA_real_,
        Pheno
      ),

      # NEW: mask the secondary trait on the candidates'
      # validation plots
      sec_masked = Trait != "GY" & is_val_plot,
      Pheno = if_else(sec_masked, NA_real_, Pheno)
    ) |>
    group_by(TraitEnv) |>
    mutate(
      # Standardize phenotypes within each trait-environment
      Pheno_z = as.vector(
        scale(Pheno, center = TRUE, scale = TRUE)
      ),
      Pheno_mean = mean(Pheno, na.rm = TRUE),
      Pheno_sd = sd(Pheno, na.rm = TRUE)
    ) |>
    ungroup() |>
    droplevels()

  # Every validation plot must have exactly one secondary-
  # trait record masked. If this fails, environment + Row +
  # Col does not identify the same plot across traits in this
  # data file.
  n_masked <- sum(df_seed$sec_masked)
  message(
    "Secondary-trait records masked: ", n_masked,
    " (validation plots: ", length(validation_plots), ")"
  )
  if (n_masked != length(validation_plots)) {
    stop(
      "Masked ", n_masked, " secondary-trait records but ",
      "expected ", length(validation_plots), "."
    )
  }

  observation_summary <- df_seed |>
    mutate(
      candidate = !is.na(Gen) &
        as.character(Gen) %in% as.character(val$Gen)
    ) |>
    group_by(Trait, candidate) |>
    summarise(
      n_total = n(),
      n_available = sum(!is.na(Pheno)),
      .groups = "drop"
    )

  print(observation_summary)

  # 5c. Fit MTME model ---------------------------------------

  fit <- asreml(
    Pheno_z ~ TraitEnv,
    random = ~
      rr(TraitEnv, k):vm(Gkeep, Ginv) +
      diag(TraitEnv):vm(Gkeep, Ginv) +
      diag(TraitEnv):ide(Gkeep) +
      diag(TraitEnv):Block,
    residual = ~ dsum(
      ~ ar1(Col):ar1(Row) | TraitEnv
    ),
    sparse = ~ TraitEnv:Gdrop,
    data = df_seed,
    na.action = na.method(x = "include"),
    maxit = 13,
    workspace = "8gb"
  )

  # 5d. Update model -----------------------------------------

  fit <- update_asreml(fit, max_updates = 20)

  if (!isTRUE(fit$converge)) {
    warning(run_name, " did not converge.")
  }

  # 5e. Save fitted model ------------------------------------

  saveRDS(fit, here("Results", paste0(run_name, ".rds")))

  rm(fit, df_seed, Pheno_f, Pheno)
  gc()
}
