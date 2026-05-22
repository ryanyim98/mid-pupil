MIDAFFEMO_ROOT <- Sys.getenv("MIDAFFEMO_ROOT", unset = "~/Desktop/midaffemo")

THRESH_STR <- c("0", "05", "1", "15", "2", "25", "3", "35")
THRESH <- c(0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5)
SPLITS <- c("run12", "run34")
TR_DELAY <- 17L
MAX_TR <- 20L
TAIL_NA_FROM_TR <- 18L
PLOT_TR_ANT <- 14L
PLOT_TR_OUT <- 13L
TRIAL_TR_SEC <- 2L
CUE_TR <- 4L

trial_tr_scale <- function(max_tr, by = 3L, cue_tr = CUE_TR, tr_sec = TRIAL_TR_SEC) {
  breaks <- seq(1L, max_tr, by = by)
  list(
    limits = c(1, max_tr),
    breaks = breaks,
    labels = (breaks - cue_tr) * tr_sec
  )
}

# Column names used across template-similarity pipelines.
TS_COL <- list(
  subject = "subject",
  subject_id = "X1",
  trial = "trial",
  probe = "probe",
  hit = "hit",
  hit_num = "hit_num",
  cue_value = "cue_value",
  cue_value_abs = "cue_value_abs",
  arousal_scaled = "arousal_scaled",
  rt = "rt",
  rt_log = "rt_log",
  tr = "tr",
  tr_in_trial = "tr_in_trial",
  trial_tr = "trial_tr",
  threshold = "threshold",
  split = "split",
  half = "half",
  time = "time",
  censor = "censor",
  type = "type",
  val = "val",
  glm_arousal = "glmArousal",
  glm_pupil = "glmPupil",
  pupil = "pupil",
  pupil_lag2 = "pupil_lag2",
  mean_pupil_nobc = "mean_pupil_nobc",
  trial_value = "trial_value",
  trial_outcome = "trial_outcome",
  run = "run",
  block = "block",
  iti = "iti"
)

# Grid column names (in expand_grid tibbles) mapped to data column names.
TS_GRID <- list(
  tr = "tr",
  threshold = "threshold",
  split = "split",
  probe = "probe"
)

TS_GRID_MAP <- list(
  tr_in_trial = TS_GRID$tr,
  threshold = TS_GRID$threshold,
  split = TS_GRID$split,
  probe = TS_GRID$probe
)

TS_GRID_MAP_POOLED <- TS_GRID_MAP[c("tr_in_trial", "threshold")]

TS_VAL <- list(
  hit = "Hit",
  probe_ant = "anticipation probe",
  probe_out = "outcome probe",
  probe_ant_short = "anticipation",
  probe_out_short = "outcome",
  half_first = "first",
  half_second = "second",
  glm_type_arousal = "glmArousal",
  glm_type_pupil = "glmPupil",
  glm_type_pupil_ts = "pupil",
  timing_subject = "ab240324",
  timing_voi = "nacc",
  timing_coreg = "ants"
)

TS_COEF <- list(
  estimate = "Estimate",
  se = "Std. Error",
  stat_glm = "z value",
  stat_lm = "t value",
  p_glm = "Pr(>|z|)",
  p_lm = "Pr(>|t|)"
)

col_name <- function(col) {
  if (is.list(col)) unlist(col, use.names = FALSE) else col
}

col_sym <- function(col) rlang::sym(col_name(col))

filter_by_grid_row <- function(data, row, col_map = TS_GRID_MAP) {
  out <- data
  for (data_col in names(col_map)) {
    grid_col <- col_map[[data_col]]
    if (!grid_col %in% names(row)) next
    out <- out %>% dplyr::filter(.data[[data_col]] == row[[grid_col]])
  }
  out
}

mid_pupil_resolve_data <- function(...) {
  parts <- file.path(...)
  tails <- c(
    file.path("data", parts),
    file.path("..", "data", parts),
    file.path("..", "..", "data", parts),
    file.path("..", "..", "..", "data", parts)
  )
  wd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  for (rel in tails) {
    cand <- normalizePath(file.path(wd, rel), winslash = "/", mustWork = FALSE)
    if (!is.na(cand) && file.exists(cand)) return(cand)
  }
  stop("Cannot find data file ", parts, " (getwd = ", getwd(), ")", call. = FALSE)
}

make_model_grid <- function(
    splits = SPLITS,
    thresholds = THRESH,
    tr_seq = seq_len(TR_DELAY),
    include_split = TRUE
) {
  grid_args <- list(
    threshold = thresholds,
    tr = tr_seq
  )
  if (include_split) {
    grid_args <- c(list(split = splits), grid_args)
  }
  do.call(tidyr::expand_grid, grid_args)
}

load_fmri3_subjects <- function(exclude = "ag240330") {
  id_col <- TS_COL$subject_id
  read_delim(
    file.path(MIDAFFEMO_ROOT, "analyses/subjects-pupil-fmri.txt"),
    delim = "\t",
    col_names = FALSE
  ) %>%
    filter(.data[[id_col]] != exclude)
}

load_motion_censor <- function(subjects) {
  id_col <- TS_COL$subject_id
  sub_col <- TS_COL$subject
  tr_col <- TS_COL$tr
  tr_in_trial_col <- TS_COL$tr_in_trial
  censor_col <- TS_COL$censor

  purrr::map_dfr(subjects[[id_col]], function(sub) {
    read_csv(file.path(MIDAFFEMO_ROOT, "svr/data/beh_index", paste0(sub, ".csv"))) %>%
      mutate(
        !!tr_in_trial_col := .data[[tr_col]],
        !!tr_col := row_number(),
        !!sub_col := sub
      )
  }) %>%
    select(all_of(c(tr_col, sub_col, censor_col)))
}

load_fmri3_glm_arrays <- function(
    template_half,
    pupil_tr_slice,
    subjects,
    thresh_str = THRESH_STR,
    thresh = THRESH,
    root = MIDAFFEMO_ROOT
) {
  subject_ids <- subjects[[TS_COL$subject_id]]
  type_col <- TS_COL$type
  tr_col <- TS_COL$tr
  threshold_col <- TS_COL$threshold
  subject_col <- TS_COL$subject
  val_col <- TS_COL$val

  purrr::map_dfr(seq_along(thresh_str), function(t) {
    arousal_path <- file.path(
      root, "glm_maps",
      paste0("glm27c_", template_half, "_array_threshold", thresh_str[t], ".csv")
    )
    pupil_path <- file.path(
      root, "glm_maps",
      paste0("glm33b_", template_half, "_array_threshold", thresh_str[t], ".csv")
    )

    glm_arousal <- as.data.frame(t(read_csv(arousal_path, col_names = FALSE)[c(1, 3:24), ])) %>%
      mutate(
        !!type_col := TS_VAL$glm_type_arousal,
        !!tr_col := row_number(),
        !!threshold_col := thresh[t]
      )

    glm_pupil <- as.data.frame(t(read_csv(pupil_path, col_names = FALSE))) %>%
      mutate(
        !!type_col := TS_VAL$glm_type_pupil,
        !!tr_col := row_number(),
        !!threshold_col := thresh[t]
      )

    pupil_ts <- as.data.frame(t(read_csv(
      file.path(root, "glm_maps/pupil_array.csv"),
      col_names = FALSE
    ))) %>%
      slice(pupil_tr_slice) %>%
      mutate(
        !!type_col := TS_VAL$glm_type_pupil_ts,
        !!tr_col := row_number(),
        !!threshold_col := thresh[t]
      )

    bind_rows(glm_arousal, glm_pupil, pupil_ts)
  }) %>%
    setNames(c(subject_ids, type_col, tr_col, threshold_col)) %>%
    pivot_longer(all_of(subject_ids), names_to = subject_col, values_to = val_col) %>%
    pivot_wider(names_from = all_of(type_col), values_from = all_of(val_col))
}

build_trial_timing_template <- function(
    runs,
    tail_runs = NULL,
    root = MIDAFFEMO_ROOT
) {
  sub_col <- TS_COL$subject
  probe_col <- TS_COL$probe
  iti_col <- TS_COL$iti
  tr_col <- TS_COL$tr
  trial_tr_col <- TS_COL$trial_tr
  run_col <- TS_COL$run
  trial_col <- TS_COL$trial
  block_col <- TS_COL$block

  if (is.null(tail_runs)) {
    tail_runs <- list(
      list(run = runs[1], tr_offset = 10L, after_row = 240L),
      list(
        run = runs[2],
        tr_offset = if (runs[2] == "run-02") 9L else 11L,
        after_row = 485L
      )
    )
  }

  df <- read_csv(file.path(
    root,
    "analyses/timecourses/data/realsubj/MIDaffemo/timecourses_b4_woOutliers_long.csv"
  )) %>%
    filter(
      .data[[sub_col]] == TS_VAL$timing_subject,
      voi == TS_VAL$timing_voi,
      coreg == TS_VAL$timing_coreg
    ) %>%
    mutate(
      trial_tr_noiti = ifelse(.data[[probe_col]] == 1, 9, 8),
      !!trial_tr_col := trial_tr_noiti + .data[[iti_col]] / 2
    ) %>%
    filter(.data[[tr_col]] <= .data[[trial_tr_col]], .data[[run_col]] %in% runs)

  for (tail in tail_runs) {
    end_run <- df %>%
      filter(.data[[run_col]] == tail$run, .data[[trial_col]] == 24) %>%
      dplyr::slice(1:5) %>%
      mutate(!!tr_col := .data[[tr_col]] + tail$tr_offset)
    df <- df %>% add_row(end_run, .after = tail$after_row)
  }

  df %>%
    mutate(
      !!block_col := as.numeric(substr(.data[[run_col]], 6, 6)),
      !!trial_col := .data[[trial_col]] + (.data[[block_col]] - 1) * 24,
      !!tr_col := row_number()
    )
}

join_glm_arrays_with_timing <- function(glm_arrays, timing, split_label, tr_offset = 0L) {
  sub_col <- TS_COL$subject
  trial_col <- TS_COL$trial
  threshold_col <- TS_COL$threshold
  tr_col <- TS_COL$tr
  tr_in_trial_col <- TS_COL$tr_in_trial
  split_col <- TS_COL$split
  trial_tr_col <- TS_COL$trial_tr
  time_col <- TS_COL$time

  out <- glm_arrays %>%
    left_join(
      timing %>% select(trial:erating_t, all_of(c(time_col, tr_col, trial_tr_col))),
      by = tr_col
    ) %>%
    group_by(.data[[sub_col]], .data[[trial_col]], .data[[threshold_col]]) %>%
    mutate(!!tr_in_trial_col := row_number()) %>%
    relocate(all_of(tr_in_trial_col), .after = all_of(tr_col)) %>%
    mutate(!!split_col := split_label) %>%
    ungroup()

  if (tr_offset != 0L) {
    out <- out %>% mutate(!!tr_col := .data[[tr_col]] + tr_offset)
  }
  out
}

template_beh_cols <- function() {
  c(
    TS_COL$hit,
    TS_COL$cue_value,
    TS_COL$arousal_scaled,
    TS_COL$rt,
    TS_COL$trial_value,
    TS_COL$trial_outcome,
    TS_COL$mean_pupil_nobc
  )
}

template_join_keys <- function() {
  c(TS_COL$subject, TS_COL$trial, TS_COL$probe)
}

complete_trial_tr_grid <- function(
    long_df,
    max_tr = MAX_TR,
    group_keys = c(TS_COL$subject, TS_COL$threshold, TS_COL$split, TS_COL$trial, TS_COL$trial_tr)
) {
  value_cols <- c(TS_COL$glm_arousal, TS_COL$glm_pupil, TS_COL$pupil)
  meta_cols <- c(group_keys, TS_COL$tr_in_trial, value_cols, TS_COL$tr, TS_COL$time)
  fill_cols <- setdiff(names(long_df), meta_cols)

  long_df %>%
    group_by(across(all_of(group_keys))) %>%
    tidyr::complete(!!rlang::sym(TS_COL$tr_in_trial) := seq_len(max_tr)) %>%
    tidyr::fill(all_of(fill_cols), .direction = "downup") %>%
    ungroup()
}

backfill_short_trial_trs <- function(
    long_df,
    value_col = TS_COL$glm_arousal,
    tr_delay = TR_DELAY,
    short_trial_trs = 9:11,
    group_keys = c(TS_COL$subject, TS_COL$threshold, TS_COL$split),
    tail_na_from_tr = TAIL_NA_FROM_TR
) {
  trial_col <- TS_COL$trial
  tr_in_trial_col <- TS_COL$tr_in_trial
  trial_tr_col <- TS_COL$trial_tr
  value_col <- col_name(value_col)

  long_df %>%
    arrange(across(all_of(c(group_keys, trial_col, tr_in_trial_col)))) %>%
    group_by(across(all_of(group_keys))) %>%
    group_modify(function(.x, .y) {
      trial_ids <- sort(unique(.x[[trial_col]]))
      for (j in seq_along(trial_ids)) {
        if (j == 1L) next
        curr_trial <- trial_ids[j]
        prev_trial <- trial_ids[j - 1L]
        curr <- .x[.x[[trial_col]] == curr_trial, , drop = FALSE]
        prev <- .x[.x[[trial_col]] == prev_trial, , drop = FALSE]
        tr_end <- unique(curr[[trial_tr_col]])[1L]
        if (!tr_end %in% short_trial_trs) next
        marker <- curr[[value_col]][curr[[tr_in_trial_col]] == tr_end + 1L]
        if (length(marker) == 0L || !is.na(marker[1L])) next
        n_carry <- tr_delay - tr_end
        for (k in seq_len(n_carry)) {
          dst_tr <- tr_end + k
          src_tr <- k
          dst_idx <- which(
            .x[[trial_col]] == curr_trial & .x[[tr_in_trial_col]] == dst_tr
          )
          src_val <- prev[[value_col]][prev[[tr_in_trial_col]] == src_tr]
          if (length(dst_idx) && length(src_val)) {
            .x[[value_col]][dst_idx] <- src_val[1L]
          }
        }
      }
      .x
    }) %>%
    ungroup() %>%
    mutate(
      !!value_col := if_else(
        .data[[tr_in_trial_col]] >= tail_na_from_tr,
        NA_real_,
        .data[[value_col]]
      )
    )
}

prepare_template_long <- function(
    my_array,
    value_col,
    fmri3_data_2probes,
    tr_max = 975L,
    half_trial_cutoff = 48L,
    max_tr = MAX_TR
) {
  join_keys <- template_join_keys()
  beh_cols <- template_beh_cols()
  beh_data <- fmri3_data_2probes %>% select(all_of(c(join_keys, beh_cols)))

  trial_col <- TS_COL$trial
  tr_col <- TS_COL$tr
  probe_col <- TS_COL$probe
  half_col <- TS_COL$half
  tr_in_trial_col <- TS_COL$tr_in_trial
  time_col <- TS_COL$time
  value_col <- col_name(value_col)

  my_array %>%
    ungroup() %>%
    select(-any_of(beh_cols)) %>%
    complete_trial_tr_grid(max_tr = max_tr) %>%
    backfill_short_trial_trs(value_col = value_col) %>%
    left_join(beh_data, by = join_keys) %>%
    filter(is.na(.data[[tr_col]]) | .data[[tr_col]] <= tr_max) %>%
    mutate(
      !!half_col := ifelse(
        .data[[trial_col]] >= half_trial_cutoff,
        TS_VAL$half_second,
        TS_VAL$half_first
      ),
      !!probe_col := ifelse(
        .data[[probe_col]] == 1,
        TS_VAL$probe_ant,
        TS_VAL$probe_out
      ),
      !!tr_in_trial_col := as.numeric(.data[[tr_in_trial_col]]),
      !!time_col := (.data[[tr_in_trial_col]] - 1) * 2,
      !!TS_COL$rt_log := ifelse(
        .data[[TS_COL$rt]] > 0,
        log(.data[[TS_COL$rt]]),
        NA_real_
      )
    )
}

add_pupil_lag2 <- function(template_df, my_array) {
  join_keys <- c(
    TS_COL$subject, TS_COL$threshold, TS_COL$split,
    TS_COL$trial, TS_COL$tr_in_trial
  )
  pupil_col <- TS_COL$pupil
  pupil_lag_col <- TS_COL$pupil_lag2
  tr_in_trial_col <- TS_COL$tr_in_trial

  pupil_lag <- my_array %>%
    group_by(
      .data[[TS_COL$subject]],
      .data[[TS_COL$threshold]],
      .data[[TS_COL$split]]
    ) %>%
    mutate(!!pupil_lag_col := lag(.data[[pupil_col]], 2)) %>%
    ungroup() %>%
    mutate(!!tr_in_trial_col := as.numeric(.data[[tr_in_trial_col]])) %>%
    select(all_of(c(pupil_lag_col, join_keys)))

  template_df %>% left_join(pupil_lag, by = join_keys)
}

cue_value_to_abs <- function(x) {
  dplyr::case_when(
    x %in% c("+$5", "-$5") ~ 5,
    x %in% c("+$1", "-$1") ~ 1,
    x %in% c("+$0", "-$0") ~ 0,
    TRUE ~ NA_real_
  )
}

predict_ci <- function(model, newdata) {
  X <- model.matrix(delete.response(terms(model)), newdata)
  beta <- lme4::fixef(model)
  V <- stats::vcov(model)

  fit <- as.vector(X %*% beta)
  se <- sqrt(diag(X %*% V %*% t(X)))

  newdata %>%
    mutate(
      fit = fit,
      lwr = fit - 1.96 * se,
      upr = fit + 1.96 * se
    )
}

extract_coef_row <- function(fit_summary, term = NULL) {
  coef_mat <- fit_summary$coefficients
  if (!is.null(term)) {
    term_row <- grep(term, rownames(coef_mat), fixed = TRUE)
    if (length(term_row) == 0L) {
      return(na_coef_row())
    }
    coef_row <- term_row[1L]
  } else {
    coef_row <- 2L
  }

  stat_col <- intersect(
    c(TS_COEF$stat_glm, TS_COEF$stat_lm),
    colnames(coef_mat)
  )[1]
  p_col <- intersect(
    c(TS_COEF$p_glm, TS_COEF$p_lm),
    colnames(coef_mat)
  )[1]

  tibble(
    Estimate = coef_mat[coef_row, TS_COEF$estimate],
    se = coef_mat[coef_row, TS_COEF$se],
    stat = if (!is.na(stat_col)) coef_mat[coef_row, stat_col] else NA_real_,
    p = if (!is.na(p_col)) coef_mat[coef_row, p_col] else NA_real_
  )
}

na_coef_row <- function() {
  tibble(Estimate = NA_real_, se = NA_real_, stat = NA_real_, p = NA_real_)
}

safe_coef_row <- function(fit_fn, term = NULL) {
  tryCatch(
    extract_coef_row(summary(fit_fn()), term = term),
    error = function(e) na_coef_row()
  )
}

scaled_term <- function(predictor) {
  paste0("scale(", col_name(predictor), ")")
}

fit_hit_grid <- function(
    data,
    grid,
    predictor,
    use_glmer = FALSE,
    col_map = TS_GRID_MAP
) {
  hit_col <- TS_COL$hit
  hit_num_col <- TS_COL$hit_num
  probe_col <- TS_COL$probe
  predictor <- col_name(predictor)
  term <- scaled_term(predictor)

  purrr::pmap_dfr(grid, function(...) {
    row <- list(...)
    dat <- filter_by_grid_row(data, row, col_map) %>%
      mutate(!!hit_num_col := as.integer(.data[[hit_col]] == TS_VAL$hit))

    purrr::imap_dfr(
      stats::setNames(
        c(TS_VAL$probe_ant, TS_VAL$probe_out),
        c(TS_VAL$probe_ant_short, TS_VAL$probe_out_short)
      ),
      function(probe_val, probe_label) {
        probe_dat <- dat %>% filter(.data[[probe_col]] == probe_val)
        form <- stats::as.formula(
          paste(TS_COL$hit_num, "~", term)
        )
        safe_coef_row(function() {
          if (use_glmer) {
            lme4::glmer(form, family = binomial, data = probe_dat)
          } else {
            stats::glm(form, family = binomial, data = probe_dat)
          }
        }, term = term) %>%
          mutate(
            probe = probe_label,
            split = row[[TS_GRID$split]],
            threshold = row[[TS_GRID$threshold]],
            tr = row[[TS_GRID$tr]]
          )
      }
    )
  }) %>%
    rename(z = stat)
}

fit_rt_grid <- function(
    data,
    grid,
    predictor,
    use_glmer = FALSE,
    col_map = TS_GRID_MAP
) {
  hit_col <- TS_COL$hit
  probe_col <- TS_COL$probe
  rt_col <- TS_COL$rt_log
  predictor <- col_name(predictor)
  term <- scaled_term(predictor)

  purrr::pmap_dfr(grid, function(...) {
    row <- list(...)
    dat <- filter_by_grid_row(data, row, col_map) %>%
      filter(.data[[hit_col]] == TS_VAL$hit)

    purrr::imap_dfr(
      stats::setNames(
        c(TS_VAL$probe_ant, TS_VAL$probe_out),
        c(TS_VAL$probe_ant_short, TS_VAL$probe_out_short)
      ),
      function(probe_val, probe_label) {
        probe_dat <- dat %>% filter(.data[[probe_col]] == probe_val)
        form <- stats::as.formula(
          paste("scale(", rt_col, ") ~", term)
        )
        safe_coef_row(function() {
          if (use_glmer) {
            lmer(form, data = probe_dat, REML = FALSE)
          } else {
            stats::lm(form, data = probe_dat)
          }
        }, term = term) %>%
          mutate(
            probe = probe_label,
            split = row[[TS_GRID$split]],
            threshold = row[[TS_GRID$threshold]],
            tr = row[[TS_GRID$tr]]
          )
      }
    )
  }) %>%
    rename(z = stat)
}

fit_hit_grid_pooled <- function(
    data,
    grid,
    predictor,
    col_map = TS_GRID_MAP_POOLED
) {
  hit_col <- TS_COL$hit
  hit_num_col <- TS_COL$hit_num
  probe_col <- TS_COL$probe
  split_col <- TS_COL$split
  predictor <- col_name(predictor)
  term <- scaled_term(predictor)

  purrr::pmap_dfr(grid, function(...) {
    row <- list(...)
    dat <- filter_by_grid_row(data, row, col_map) %>%
      mutate(!!hit_num_col := as.integer(.data[[hit_col]] == TS_VAL$hit))

    purrr::imap_dfr(
      stats::setNames(
        c(TS_VAL$probe_ant, TS_VAL$probe_out),
        c(TS_VAL$probe_ant_short, TS_VAL$probe_out_short)
      ),
      function(probe_val, probe_label) {
        probe_dat <- dat %>% filter(.data[[probe_col]] == probe_val)
        form <- stats::as.formula(
          paste(TS_COL$hit_num, "~", term, "+ (1|", split_col, ")", sep = "")
        )
        safe_coef_row(function() {
          lme4::glmer(form, family = binomial, data = probe_dat)
        }, term = term) %>%
          mutate(
            probe = probe_label,
            threshold = row[[TS_GRID$threshold]],
            tr = row[[TS_GRID$tr]]
          )
      }
    )
  }) %>%
    rename(z = stat)
}

fit_rt_grid_pooled <- function(
    data,
    grid,
    predictor,
    col_map = TS_GRID_MAP_POOLED
) {
  hit_col <- TS_COL$hit
  probe_col <- TS_COL$probe
  rt_col <- TS_COL$rt_log
  split_col <- TS_COL$split
  predictor <- col_name(predictor)
  term <- scaled_term(predictor)

  purrr::pmap_dfr(grid, function(...) {
    row <- list(...)
    dat <- filter_by_grid_row(data, row, col_map) %>%
      filter(.data[[hit_col]] == TS_VAL$hit)

    purrr::imap_dfr(
      stats::setNames(
        c(TS_VAL$probe_ant, TS_VAL$probe_out),
        c(TS_VAL$probe_ant_short, TS_VAL$probe_out_short)
      ),
      function(probe_val, probe_label) {
        probe_dat <- dat %>% filter(.data[[probe_col]] == probe_val)
        form <- stats::as.formula(
          paste("scale(", rt_col, ") ~", term, "+ (1|", split_col, ")", sep = "")
        )
        safe_coef_row(function() {
          lmer(form, data = probe_dat, REML = FALSE)
        }, term = term) %>%
          mutate(
            probe = probe_label,
            threshold = row[[TS_GRID$threshold]],
            tr = row[[TS_GRID$tr]]
          )
      }
    )
  }) %>%
    rename(z = stat)
}

fit_cue_grid <- function(
    data,
    grid,
    predictor,
    use_glmer = FALSE,
    col_map = if ("split" %in% names(grid)) TS_GRID_MAP else TS_GRID_MAP_POOLED
) {
  probe_col <- TS_COL$probe
  cue_col <- TS_COL$cue_value
  cue_abs_col <- TS_COL$cue_value_abs
  split_col <- TS_COL$split
  predictor <- col_name(predictor)
  term <- scaled_term(predictor)

  grid_with_probe <- tidyr::crossing(
    grid,
    !!rlang::sym(probe_col) := c(TS_VAL$probe_ant, TS_VAL$probe_out)
  )

  purrr::pmap_dfr(grid_with_probe, function(...) {
    row <- list(...)
    row_col_map <- c(col_map, probe = probe_col)
    dat <- filter_by_grid_row(data, row, row_col_map) %>%
      mutate(!!cue_abs_col := cue_value_to_abs(.data[[cue_col]]))

    form <- if (use_glmer) {
      stats::as.formula(
        paste(cue_abs_col, "~", term, "+ (1|", split_col, ")", sep = "")
      )
    } else {
      stats::as.formula(paste("scale(", cue_abs_col, ") ~", term))
    }

    out <- safe_coef_row(function() {
      if (use_glmer) {
        lmer(form, data = dat, REML = FALSE)
      } else {
        stats::lm(form, data = dat)
      }
    }, term = term)

    out[[probe_col]] <- sub(" probe$", "", row[[probe_col]])
    for (grid_col in setdiff(names(row), probe_col)) {
      out[[grid_col]] <- row[[grid_col]]
    }
    out
  }) %>%
    rename(t = stat)
}

fit_rating_grid <- function(
    data,
    grid,
    predictor,
    use_glmer = FALSE,
    col_map = if ("split" %in% names(grid)) TS_GRID_MAP else TS_GRID_MAP_POOLED
) {
  probe_col <- TS_COL$probe
  arousal_col <- TS_COL$arousal_scaled
  split_col <- TS_COL$split
  predictor <- col_name(predictor)
  term <- scaled_term(predictor)

  grid_with_probe <- tidyr::crossing(
    grid,
    !!rlang::sym(probe_col) := c(TS_VAL$probe_ant, TS_VAL$probe_out)
  )

  purrr::pmap_dfr(grid_with_probe, function(...) {
    row <- list(...)
    row_col_map <- c(col_map, probe = probe_col)
    dat <- filter_by_grid_row(data, row, row_col_map)

    form <- if (use_glmer) {
      stats::as.formula(
        paste("scale(", arousal_col, ") ~", term, "+ (1|", split_col, ")", sep = "")
      )
    } else {
      stats::as.formula(paste("scale(", arousal_col, ") ~", term))
    }

    out <- safe_coef_row(function() {
      if (use_glmer) {
        lmer(form, data = dat, REML = FALSE)
      } else {
        stats::lm(form, data = dat)
      }
    }, term = term)

    out[[probe_col]] <- sub(" probe$", "", row[[probe_col]])
    for (grid_col in setdiff(names(row), probe_col)) {
      out[[grid_col]] <- row[[grid_col]]
    }
    out
  }) %>%
    rename(t = stat)
}
