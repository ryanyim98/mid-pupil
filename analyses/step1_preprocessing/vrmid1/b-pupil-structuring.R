source("~/Desktop/VRMID-analysis/mid-pupil/analyses/load_libraries.R")

bad_participants <- read.csv("~/Desktop/VRMID-analysis/mid-pupil/data/vrmid1/subjects_list/subject_to_drop_pupil.txt",
                             header = F)


# data w baseline correction
data_pupil <- read_csv("~/Desktop/VRMID-analysis/mid-pupil/data/vrmid1/derivatives/pupillometry.csv")%>%
  filter(!Subject %in% bad_participants)%>%
  group_by(Subject, Time)%>%
  mutate(sample_in_sec = row_number())%>%
  ungroup()%>%
  group_by(Subject,trial)%>%
  mutate(sample_in_trial_n = row_number(),
         sample_in_trial_t = sample_in_trial_n/120)

data_pupil <- data_pupil %>%
  group_by(Subject) %>%
  mutate(AvgPPos_x_d = c(NA, diff(AvgPPos_x)),
         AvgPPos_y_d = c(NA, diff(AvgPPos_y))) %>%
  ungroup()

data_pupil$pupil_L_xyvar <- NA_real_
data_pupil$pupil_R_xyvar <- NA_real_
data_pupil$pupil_Avg_xyvar <- NA_real_
xyvar_fit_summary <- list()

for (sub in unique(data_pupil$Subject)){
  ind <- which(data_pupil$Subject == sub)
  d <- data_pupil[ind, ]
  for (signal in c("L","R","Avg")){
    y_col <- paste0("pupil_", signal)
    res_col <- paste0("pupil_", signal, "_xyvar")
    fit_formula <- reformulate(
      c("AvgPPos_x", "AvgPPos_y", "AvgPPos_x_d", "AvgPPos_y_d"),
      response = y_col
    )
    fit <- tryCatch(
      lm(fit_formula, data = d, na.action = na.exclude),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      warning(paste0("xyvar regression failed for subject ", sub,
                     " signal ", signal, "; leaving NAs"))
      next
    }
    gs <- summary(fit)
    fstat <- gs$fstatistic
    if (is.null(fstat)) {
      f_val <- NA_real_
      f_p   <- NA_real_
    } else {
      f_val <- unname(fstat["value"])
      f_p   <- unname(pf(fstat["value"], fstat["numdf"], fstat["dendf"], lower.tail = FALSE))
    }
    xyvar_fit_summary[[paste(sub, signal, sep = "_")]] <- tibble(
      Subject = sub,
      pupil_signal = signal,
      n_obs = nobs(fit),
      r2 = gs$r.squared,
      adj_r2 = gs$adj.r.squared,
      F = f_val,
      F_pvalue = f_p
    )
    data_pupil[[res_col]][ind] <- as.numeric(residuals(fit))
  }
}

xyvar_fit_summary <- bind_rows(xyvar_fit_summary)
write_csv(xyvar_fit_summary,
          "~/Desktop/VRMID-analysis/mid-pupil/data/vrmid1/derivatives/pupil_xyvar_fit_summary.csv")

data_pupil <- data_pupil %>%
  group_by(Subject) %>%
  mutate(pupil_L_scaled = as.numeric(scale(pupil_L_xyvar)),
         pupil_R_scaled = as.numeric(scale(pupil_R_xyvar)),
         pupil_Avg_scaled = as.numeric(scale(pupil_Avg_xyvar)))

data_pupil%>%
  head(100)

data_pupil$trial_type <- factor(data_pupil$trial_type,
                                levels = c("-$5","-$1","-$0",
                                           "+$0","+$1","+$5"))

#baseline stuff
data_pupil_out <- {}
# get the last two seconds of previous ITI and paste it to this trial
for (i in 1:length(unique(data_pupil$Subject))){
  sub = unique(data_pupil$Subject)[i]
  print(i)
  temp_data <- data_pupil%>%
    filter(Subject == sub)%>%
    ungroup()
  for (j in 1:length(unique(temp_data$trial))){
    temp_trial_data <- temp_data%>%
      filter(trial == j)
    
    if (j == 1){
      d <- temp_trial_data
    } else{
      
      last_trial_iti_data <- temp_data%>%
        filter(trial == (j-1),
               current_stimulus == "ITI")
      
      # last 2 second
      last_trial_iti_data <-last_trial_iti_data%>%
        filter(Time_sec >= max(last_trial_iti_data$Time_sec)-1)
      
      # change timing info
      last_trial_iti_data <- last_trial_iti_data%>%
        mutate(Time_sec = Time_sec-max(last_trial_iti_data$Time_sec),
               Time_in_trial_sec = Time_in_trial_sec - max(last_trial_iti_data$Time_in_trial_sec),
               sample_in_trial_n = sample_in_trial_n - max(last_trial_iti_data$sample_in_trial_n),
               sample_in_trial_t = sample_in_trial_t - max(last_trial_iti_data$sample_in_trial_t),
               current_stimulus = "prestim_baseline")%>%
        select(Subject,Time,Time_sec,current_stimulus,Time_str,Time_str_start,Time_in_trial_sec,sample_in_sec,sample_in_trial_n:sample_in_trial_t,
               pupil_L,pupil_R,pupil_Avg,
               pupil_L_xyvar,pupil_R_xyvar,pupil_Avg_xyvar,
               pupil_L_scaled,pupil_R_scaled,pupil_Avg_scaled)
      
      d <- merge(temp_trial_data,last_trial_iti_data,
                 all.x = T, all.y = T)%>%
        arrange(sample_in_trial_n)%>%
        fill(trial,trial_type:NegA_scaled, .direction = "downup")%>%
        relocate(trial,.after = Subject)
    }
    
    
    data_pupil_out <- rbind(data_pupil_out,d)
  }
}


#get the mean value of the last second in a trial (ITI)
baseline <- data_pupil_out%>%
  group_by(Subject, trial)%>%
  filter(Time_sec == 0)%>%
  summarise_at(vars(pupil_L,pupil_R,pupil_Avg,pupil_L_scaled,pupil_R_scaled,pupil_Avg_scaled), 
               ~mean(.x, na.rm = TRUE))%>%
  ungroup()

data_pupil_out$pupil_L_baseline <- NA
data_pupil_out$pupil_R_baseline <- NA
data_pupil_out$pupil_Avg_baseline <- NA

data_pupil_out$pupil_L_baseline_scaled <- NA
data_pupil_out$pupil_R_baseline_scaled <- NA
data_pupil_out$pupil_Avg_baseline_scaled <- NA

for (i in 1:nrow(baseline)){
  ind = which(data_pupil_out$Subject == baseline$Subject[i] &
                data_pupil_out$trial == baseline$trial[i])
  #baseline info
  if (baseline$trial[i] == 1){
    data_pupil_out$pupil_L_baseline[ind] = NA
    data_pupil_out$pupil_R_baseline[ind] = NA
    data_pupil_out$pupil_Avg_baseline[ind] =  NA
  } else {
    #baseline is from last trial's ITI
    data_pupil_out$pupil_L_baseline[ind] = baseline$pupil_L[i]
    data_pupil_out$pupil_R_baseline[ind] = baseline$pupil_R[i]
    data_pupil_out$pupil_Avg_baseline[ind] = baseline$pupil_Avg[i]
    data_pupil_out$pupil_L_baseline_scaled[ind] = baseline$pupil_L_scaled[i]
    data_pupil_out$pupil_R_baseline_scaled[ind] = baseline$pupil_R_scaled[i]
    data_pupil_out$pupil_Avg_baseline_scaled[ind] = baseline$pupil_Avg_scaled[i]
  }
}

#baseline correction; first trial of each participant was excluded
data_pupil_out <- data_pupil_out%>%
  mutate(pupil_L_bc = pupil_L - pupil_L_baseline,
         pupil_R_bc = pupil_R - pupil_R_baseline,
         pupil_Avg_bc = pupil_Avg - pupil_Avg_baseline,
         pupil_L_bc_scaled = pupil_L_scaled - pupil_L_baseline_scaled,
         pupil_R_bc_scaled = pupil_R_scaled - pupil_R_baseline_scaled,
         pupil_Avg_bc_scaled = pupil_Avg_scaled - pupil_Avg_baseline_scaled)

write_csv(data_pupil_out, "~/Desktop/VRMID-analysis/mid-pupil/data/vrmid1/derivatives/pupillometry_lowpass_baselineCorrected.csv")

remove(data_pupil)
