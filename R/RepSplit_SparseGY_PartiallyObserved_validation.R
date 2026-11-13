# ============================================================
# Replicate-Split Validation - Sparse-GY Training Design
# Partially Observed Candidates: validation BLUPs and
# predictive ability
# ============================================================
#
# Purpose:
# For each split (A, B):
#   1. Fits the validation BLUP model (no relationship
#      matrix) with the candidates' GY kept only on their
#      validation plots.
#   2. Extracts candidate GEBVs for GY from the replicate-
#      split GS models and from the GY-only baseline.
#   3. Correlates GEBVs with the validation BLUPs, within and
#      across environments.
# Results are then averaged over the two splits.
#
# Run RepSplit_SparseGY_PartiallyObserved_models.R first.
# The BLUP models here have no vm() term.
#
# ============================================================

# 1. Load packages -------------------------------------------

library(tidyverse)
library(asreml)
library(here)

# 2. Load supporting functions -------------------------------

source(here("R", "Functions_earlyGS.R"))

# 3. Settings ------------------------------------------------

k_blup <- 2  # FA rank of the validation BLUP model

# GY-only baseline for sparse-GY, candidates' GY masked.
# EDIT if your file has a different name.
baseline_file <- "GYonly_SparseGY_Unobserved_FA2.rds"

# Replicate-split GS models (from the models script)
gs_models <- tribble(
  ~Model,   ~k,
  "IGYMLR", 4,
  "IGYMF",  5,
  "NDRE",   5,
  "PC1",    5
)

# Original partially observed models, which also used the
# secondary trait on the validation plots. Scored against the
# same BLUPs for reference: the difference from the replicate-
# split model shows what the validation-plot records added.
original_files <- c(
  IGYMLR = "IGYMLR_SparseGY_PartiallyObserved_FA4.rds",
  IGYMF  = "IGYMF_SparseGY_PartiallyObserved_FA5.rds",
  NDRE   = "NDRE_SparseGY_PartiallyObserved_FA5.rds",
  PC1    = "PC1_SparseGY_PartiallyObserved_FA5.rds"
)

# 4. Load inputs ---------------------------------------------

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

# All observed GY (the sparse-GY zones are NOT masked here:
# validation uses every available GY record)
GY <- readRDS(here("Data", "Pheno_SparseGY_Observed.rds")) |>
  ungroup() |>
  filter(Trait == "GY") |>
  mutate(plot_key = paste(Env, Row, Col, sep = "_")) |>
  droplevels()

roles <- make_roles(GY, candidates)

message(
  "Candidate plots: ", nrow(roles),
  " | validation plots per split: ",
  sum(roles$role_splitA == "validation")
)

dir.create(here("Results"), showWarnings = FALSE, recursive = TRUE)
write.csv(
  roles,
  here("Results", "RepSplit_plot_roles.csv"),
  row.names = FALSE
)

# 5. Functions -----------------------------------------------

# Validation BLUPs for one split
fit_validation_blups <- function(split) {

  role_col <- paste0("role_split", split)
  predictor_plots <- roles$plot_key[roles[[role_col]] == "predictor"]
  validation_plots <- roles$plot_key[roles[[role_col]] == "validation"]

  df <- GY |>
    mutate(
      # Mask candidate GY on predictor plots; keep it on
      # validation plots
      Pheno = if_else(
        plot_key %in% predictor_plots,
        NA_real_,
        Pheno
      ),
      Gen = factor(Gen),
      TraitEnv = factor(TraitEnv),
      Block = factor(Block),
      Row = factor(Row),
      Col = factor(Col)
    ) |>
    group_by(TraitEnv) |>
    mutate(
      Pheno_z = as.vector(scale(Pheno, center = TRUE, scale = TRUE))
    ) |>
    ungroup() |>
    arrange(TraitEnv, Col, Row) |>
    droplevels() |>
    as.data.frame()

  n_cand_gy <- sum(
    !is.na(df$Pheno) & as.character(df$Gen) %in% candidates
  )
  message(
    "Split ", split, ": candidate GY records kept = ", n_cand_gy,
    " (validation plots: ", length(validation_plots), ")"
  )

  m <- asreml(
    Pheno_z ~ TraitEnv,
    random = ~
      rr(TraitEnv, k_blup):Gen +
      diag(TraitEnv):Gen +
      diag(TraitEnv):Block,
    residual = ~ dsum(
      ~ ar1(Col):ar1(Row) | TraitEnv
    ),
    data = df,
    na.action = na.method(x = "include"),
    maxit = 13,
    workspace = "4gb"
  )
  m <- update_asreml(m, max_updates = 20)
  if (!isTRUE(m$converge)) {
    warning("Validation BLUP model, split ", split,
            ", did not converge.")
  }
  saveRDS(
    m,
    here("Results",
         paste0("ValidationBLUP_RepSplit", split, "_FA", k_blup, ".rds"))
  )

  # Line effect in an environment = rr() effect + diag() effect
  rand <- m$coefficients$random
  rand_df <- tibble(
    term = rownames(rand),
    effect = as.numeric(rand[, 1])
  )

  blups <- map_dfr(levels(df$TraitEnv), \(e) {
    key <- paste0("_", e, ":Gen_")
    rand_df |>
      filter(str_detect(term, fixed(key))) |>
      mutate(
        TraitEnv = e,
        Gen = str_remove(term, paste0("^.*", str_escape(key))),
        part = if_else(str_starts(term, fixed("rr(")), "rr", "diag")
      )
  }) |>
    group_by(TraitEnv, Gen) |>
    summarise(
      n_terms = n_distinct(part),
      BLUP_z = sum(effect),
      .groups = "drop"
    )

  if (nrow(blups) == 0 || any(blups$n_terms != 2)) {
    print(head(rownames(rand), 10))
    stop("Validation BLUP effects did not parse as rr + diag.")
  }

  # Which candidates have a validation plot in each environment
  has_val <- df |>
    filter(plot_key %in% validation_plots) |>
    distinct(TraitEnv, Gen) |>
    mutate(
      TraitEnv = as.character(TraitEnv),
      Gen = as.character(Gen),
      has_validation_plot = TRUE
    )

  blups |>
    filter(Gen %in% candidates) |>
    left_join(has_val, by = c("TraitEnv", "Gen")) |>
    mutate(
      has_validation_plot = replace_na(has_validation_plot, FALSE),
      split = split
    )
}

# GY GEBVs of the candidates from one fitted GS model.
# detect_fa_rank() and extract_gy_gebvs() are copied unchanged
# from Evaluate_PredictiveAbility.R, so GEBVs are rebuilt
# exactly as in the manuscript: factor scores x loadings plus
# the environment-specific genomic effect.

detect_fa_rank <- function(mod) {

  names_vp <- names(mod$vparameters)

  loading_names <- grep(
    "^rr\\(.*\\):vm\\(Gkeep, *Ginv\\)!.*!fa[0-9]+$",
    names_vp,
    value = TRUE
  )

  if (!length(loading_names)) {
    stop("No genomic FA loadings found in model.")
  }

  factor_numbers <- as.integer(
    sub(".*!fa([0-9]+)$", "\\1", loading_names)
  )

  K <- max(factor_numbers)

  if (!identical(
    sort(unique(factor_numbers)),
    seq_len(K)
  )) {
    stop("Unexpected factor numbering in model.")
  }

  K
}

extract_gy_gebvs <- function(mod) {

  K <- detect_fa_rank(mod)

  vparams <- mod$vparameters
  coefs <- mod$coefficients$random
  cn <- rownames(coefs)

  if (is.null(cn)) {
    stop("ASReml random coefficients lack row names.")
  }

  # FA loadings
  load_i <- grep(
    "^rr\\(.*\\):vm\\(Gkeep, *Ginv\\)!.*!fa[0-9]+$",
    names(vparams)
  )

  load_names <- names(vparams)[load_i]

  te_load <- sub(
    "^.*!([^!]+)!fa[0-9]+$",
    "\\1",
    load_names
  )

  comp_load <- as.integer(
    sub(
      "^.*!fa([0-9]+)$",
      "\\1",
      load_names
    )
  )

  te <- unique(te_load)

  if (
    length(load_i) != length(te) * K ||
    !identical(
      comp_load,
      rep(seq_len(K), each = length(te))
    ) ||
    !identical(
      te_load,
      rep(te, times = K)
    )
  ) {
    stop("FA loadings have unexpected structure.")
  }

  Lam <- matrix(
    as.numeric(vparams[load_i]),
    nrow = length(te),
    ncol = K
  )

  # Genotype scores
  score_i <- grep(
    "^rr\\(.*\\)_Comp[0-9]+:vm\\(Gkeep, *Ginv\\)_",
    cn
  )

  score_names <- cn[score_i]

  comp_score <- as.integer(
    sub(
      "^.*_Comp([0-9]+):vm.*$",
      "\\1",
      score_names
    )
  )

  score_gen <- sub(
    "^.*:vm\\(Gkeep, *Ginv\\)_",
    "",
    score_names
  )

  gen <- score_gen[comp_score == 1L]

  if (
    !length(gen) ||
    length(score_i) != length(gen) * K ||
    !identical(
      comp_score,
      rep(seq_len(K), each = length(gen))
    ) ||
    !identical(
      score_gen,
      rep(gen, times = K)
    )
  ) {
    stop("Unexpected FA genotype score ordering.")
  }

  F <- matrix(
    as.numeric(coefs[score_i]),
    nrow = length(gen),
    ncol = K
  )

  # Environment-specific genomic effects
  delta_i <- grep(
    "^TraitEnv_.*:vm\\(Gkeep, *Ginv\\)_",
    cn
  )

  delta_names <- cn[delta_i]

  delta_gen <- sub(
    "^.*:vm\\(Gkeep, *Ginv\\)_",
    "",
    delta_names
  )

  delta_te <- sub(
    "^TraitEnv_(.*):vm\\(Gkeep, *Ginv\\)_.*$",
    "\\1",
    delta_names
  )

  if (
    length(delta_i) != length(te) * length(gen) ||
    !identical(
      delta_gen,
      rep(gen, times = length(te))
    ) ||
    !identical(
      delta_te,
      rep(te, each = length(gen))
    )
  ) {
    stop("Unexpected specific-effect ordering.")
  }

  # Reconstruct genomic breeding values
  gebv <- as.vector(
    F %*% t(Lam)
  ) + as.numeric(coefs[delta_i])

  predictions <- tibble(
    TraitEnv = delta_te,
    TraitEnv_key = str_to_lower(delta_te),
    Gen = delta_gen,
    gebv = gebv
  ) |>
    filter(
      str_detect(TraitEnv_key, "^gy(?:-|$)")
    )

  list(
    K = K,
    predictions = predictions
  )
}

get_gy_gebv <- function(file) {
  mod <- readRDS(here("Results", file))
  if (!isTRUE(mod$converge)) {
    warning("Model is not converged: ", file)
  }
  extract_gy_gebvs(mod)$predictions |>
    filter(Gen %in% candidates) |>
    transmute(TraitEnv, Gen, GEBV_z = gebv)
}

# Predictive ability of one set of GEBVs against one split's
# validation BLUPs
score <- function(gebv, blups) {
  d <- inner_join(gebv, blups, by = c("TraitEnv", "Gen"))

  overall <- d |>
    group_by(Gen) |>
    summarise(
      GEBV = mean(GEBV_z),
      BLUP = mean(BLUP_z),
      .groups = "drop"
    ) |>
    summarise(TraitEnv = "Across", n = n(), r = cor(GEBV, BLUP))

  # Within environment: only candidates with a validation
  # plot in that environment
  within <- d |>
    filter(has_validation_plot) |>
    group_by(TraitEnv) |>
    summarise(n = n(), r = cor(GEBV_z, BLUP_z), .groups = "drop")

  # Mean of the within-environment correlations, as in the
  # manuscript's "Within Env" metric
  within_mean <- within |>
    summarise(TraitEnv = "Within (mean)", n = sum(n), r = mean(r))

  bind_rows(overall, within_mean, within)
}

# 6. Validation BLUPs for both splits ------------------------

blups_A <- fit_validation_blups("A")
blups_B <- fit_validation_blups("B")
blups_split <- list(A = blups_A, B = blups_B)

write.csv(
  bind_rows(blups_A, blups_B),
  here("Results", "ValidationBLUP_RepSplit_candidates.csv"),
  row.names = FALSE
)

# 7. Score the models ----------------------------------------

results <- map_dfr(c("A", "B"), \(s) {

  blups <- blups_split[[s]]

  # Baseline: no secondary trait, so one fit serves both splits
  out <- score(get_gy_gebv(baseline_file), blups) |>
    mutate(Model = "GYonly", version = "baseline")

  for (j in seq_len(nrow(gs_models))) {
    mod <- gs_models$Model[j]

    f_split <- paste0(
      mod, "_SparseGY_RepSplit", s, "_FA", gs_models$k[j], ".rds"
    )
    if (file.exists(here("Results", f_split))) {
      out <- bind_rows(
        out,
        score(get_gy_gebv(f_split), blups) |>
          mutate(Model = mod, version = "replicate_split")
      )
    } else {
      message("Missing: ", f_split)
    }

    f_orig <- original_files[[mod]]
    if (file.exists(here("Results", f_orig))) {
      out <- bind_rows(
        out,
        score(get_gy_gebv(f_orig), blups) |>
          mutate(Model = mod, version = "original_all_plots")
      )
    }
  }

  out |> mutate(split = s)
})

# 8. Gain over baseline, averaged over splits ----------------

baseline <- results |>
  filter(Model == "GYonly") |>
  select(split, TraitEnv, r_baseline = r)

results_gain <- results |>
  filter(Model != "GYonly") |>
  left_join(baseline, by = c("split", "TraitEnv")) |>
  mutate(gain = r - r_baseline)

summary_tab <- results_gain |>
  group_by(TraitEnv, Model, version) |>
  summarise(
    n = mean(n),
    r = mean(r),
    r_baseline = mean(r_baseline),
    gain = mean(gain),
    .groups = "drop"
  ) |>
  mutate(gain_pct = 100 * gain / r_baseline) |>
  arrange(
    match(TraitEnv, c("Across", "Within (mean)"), nomatch = 3L),
    TraitEnv, version, Model
  )

print(summary_tab, n = Inf, width = Inf)

# 9. Save ----------------------------------------------------

write.csv(
  results_gain,
  here("Results", "RepSplit_predictive_ability_bySplit.csv"),
  row.names = FALSE
)

write.csv(
  summary_tab,
  here("Results", "RepSplit_predictive_ability_summary.csv"),
  row.names = FALSE
)
