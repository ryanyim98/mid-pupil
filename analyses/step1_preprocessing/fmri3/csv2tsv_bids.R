library(tidyverse)

bids_root <- "~/Desktop/VRMID-analysis/mid-pupil/data/fmri3/bids-events"

make_trial_events <- function(trialonset, iti, probe, rating_type, arating) {
  
  if (probe == 1 & rating_type == 1) {
    tibble(
      onset = c(
        trialonset,
        trialonset + 2,
        trialonset + 4,
        trialonset + 8,
        trialonset + 12,
        trialonset + 14,
        trialonset + 16,
        trialonset + 18
      ),
      duration = c(
        2,
        2,
        4,
        4,
        2,
        2,
        2,
        iti
      ),
      stage = c(
        "cue",
        "fixation",
        "vrating",
        "arating",
        "fixation",
        "target",
        "outcome",
        "iti"
      )
    )
    
  } else if (probe == 2 & rating_type == 2) {
    tibble(
      onset = c(
        trialonset,
        trialonset + 2,
        trialonset + 4,
        trialonset + 12,
        trialonset + 14,
        trialonset + 16,
        trialonset + 18
      ),
      duration = c(
        2,
        2,
        8,
        2,
        2,
        2,
        iti
      ),
      stage = c(
        "cue",
        "fixation",
        "erating",
        "fixation",
        "target",
        "outcome",
        "iti"
      )
    )
    
  } else if (probe == 2 & rating_type == 1) {
    tibble(
      onset = c(
        trialonset,
        trialonset + 2,
        trialonset + 4,
        trialonset + 6,
        trialonset + 8,
        trialonset + 12,
        trialonset + 16
      ),
      duration = c(
        2,
        2,
        2,
        2,
        4,
        4,
        iti
      ),
      stage = c(
        "cue",
        "fixation",
        "target",
        "outcome",
        "vrating",
        "arating",
        "iti"
      )
    )
    
  } else if (probe == 1 & rating_type == 2) {
    tibble(
      onset = c(
        trialonset,
        trialonset + 2,
        trialonset + 4,
        trialonset + 6,
        trialonset + 8,
        trialonset + 16
      ),
      duration = c(
        2,
        2,
        2,
        2,
        8,
        iti
      ),
      stage = c(
        "cue",
        "fixation",
        "target",
        "outcome",
        "erating",
        "iti"
      )
    )
  }
}


for (sub in fmri3_subjects) {
  print(sub)
  
  # convert to BIDS-style subject label if needed
  # e.g., "ab240324" -> "sub-ab240324"
  sub_bids <- paste0("sub-", sub)
  
  for (run in 1:4) {
    
    filestr <- paste0(sub, "_b", run)
    
    df.beh.temp <- read_csv(
      paste0("~/Desktop/VRMID-analysis/mid-pupil/data/fmri3/raw/behavior/",
             filestr, ".csv")
    ) %>%
      filter(TR == 1)
    
    # # scale ratings conditionally
    # if (sum(is.na(df.beh.temp$arating)) > 22) {
    #   df.beh.temp <- df.beh.temp %>%
    #     mutate(
    #       arousal_scaled = as.numeric(scale(arating)),
    #       valence_scaled = as.numeric(scale(vrating))
    #     )
    # } else {
    #   df.beh.temp <- df.beh.temp %>%
    #     mutate(
    #       arousal_scaled = NA_real_,
    #       valence_scaled = NA_real_
    #     )
    # }
    
    df.beh.temp <- df.beh.temp %>%
      select(-TR) %>%
      mutate(
        last_iti = lag(iti),
        run = run
      )
    
    if ("trial_tr" %in% names(df.beh.temp)) {
      df.beh.temp <- df.beh.temp %>% select(-trial_tr)
    }
    
    # -------------------------
    # BIDS-COMPLIANT EVENTS DF
    # -------------------------
    
    df.events <- df.beh.temp %>%
      rowwise() %>%
      mutate(
        events = list(
          make_trial_events(
            trialonset  = trialonset,
            iti         = iti,
            probe       = probe,
            rating_type = rating_type,
            arating     = arating
          )
        )
      ) %>%
      unnest(events) %>%
      ungroup() %>% 
      rename(trial_type = trialtype) %>% 
      relocate(onset,duration,stage) %>% 
      mutate(onset = round(as.numeric(onset),2))
    
    df.events <- df.events %>%
      mutate(across(
        everything(),
        ~ ifelse(is.na(.) | is.nan(.), "n/a", as.character(.))
      ))
    
    # output path
    out_dir <- file.path(
      bids_root
    )
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    out_file <- file.path(
      out_dir,
      paste0(
        sub_bids,
        "_task-midaffemo_run-",
        run,
        "_events.tsv"
      )
    )
    
    write_tsv(df.events, out_file, na = "n/a")
  }
}



participants_tsv <- fmri3_demographics %>%
  mutate(
    participant_id = paste0("sub-", subject)  # change `sub` if needed
  ) %>%
  select(
    participant_id,
    randomization,
    gender,
    age,
    race,
    `MID compensation sum`,
    -subject  # drop raw subject column
  ) %>% 
  rename(sex=gender,mid_bonus = `MID compensation sum`)

write_tsv(
  participants_tsv,
  "~/Desktop/VRMID-analysis/mid-pupil/data/fmri3/participants.tsv"
)
