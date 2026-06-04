# get luminance time course for regression (downsample to TR = 2 sec)

fmri3_pupil_data.for.regs <- fmri3_pupil_data%>%
  mutate(tr_in_trial = floor((Time_sec+2)/2))  %>% 
  # filter(trial_tr > 0)%>% 
  group_by(subject,block,trial,probe,rating_type,
           tr_in_trial)%>%
  summarise(luminance= mean(luminance, na.rm = T),
            pupil_x = mean(pupil_x_preproc,na.rm = T),
            pupil_y = mean(pupil_y_preproc,na.rm = T),
            blink = mean(blink,na.rm = T),
            sacc = mean(sacc,na.rm = T),
            mean_pupil_nobc_scaled= mean(pupilDiameter_scaled, na.rm = T))

df.1sub <- read_csv("~/Desktop/midaffemo/analyses/timecourses/data/realsubj/MIDaffemo/timecourses_b0_woOutliers_long.csv") %>% 
  filter(subject == "nq240330", voi == "nacc",coreg == "ants") %>%  #participants have the same trial structure despite different randomization
  mutate(trial_tr_noiti = ifelse(probe == 1, 9,8),
         trial_tr = trial_tr_noiti + iti/2) %>% 
  select(trial,trial_tr,tr,run) %>% 
  filter(tr <= trial_tr) %>% 
  rename(tr_in_trial = tr)

end_run1 <- df.1sub %>%
  subset(trial == 24) %>%
  .[1:5,] %>%
  mutate(tr_in_trial = tr_in_trial + 10)

# create the 5 TRs at the tail of run 2
end_run2 <- df.1sub %>%
  filter(run == "run-02")%>%
  subset(trial == 24) %>%
  .[1:5,] %>%
  mutate(tr_in_trial = tr_in_trial + 9)

# create the 5 TRs at the tail of run 3
end_run3 <- df.1sub %>%
  filter(run == "run-03")%>%
  subset(trial == 24) %>%
  .[1:5,] %>%
  mutate(tr_in_trial = tr_in_trial + 10)

# create the 5 TRs at the tail of run 4
end_run4 <- df.1sub %>%
  filter(run == "run-04")%>%
  subset(trial == 24) %>%
  .[1:5,] %>%
  mutate(tr_in_trial = tr_in_trial + 11)

df.1sub <- df.1sub %>%
  # insert 5 rows at the end of each run
  add_row(end_run1, .after = 240) %>%
  add_row(end_run2, .after = 485) %>%
  add_row(end_run3, .after = 730) %>%
  add_row(end_run4, .after = 975)

df.1sub <- df.1sub %>% 
  mutate(block = as.numeric(substr(run,6,6))) %>% 
  mutate(trial = trial + (block - 1) * 24) %>% 
  select(-block,-run) %>% 
  mutate(tr = row_number())

for (sub in unique(fmri3_pupil_data.for.regs$subject)){
  if(!is.na(sub)){
    df.pupil.for.regs_ <- left_join(df.1sub,fmri3_pupil_data.for.regs %>% 
                                      filter(subject == sub),
                                    by = c("trial","tr_in_trial")) %>% 
      relocate(subject)
    
    df.pupil.for.regs_$luminance[is.na(df.pupil.for.regs_$luminance)]  <- 0.002513428#insert approximate luminance for end-of-block screen
    
    # linear interpolation; rule=2 extrapolates leading/trailing NAs with first/last non-NA
    int <- zoo::na.approx(df.pupil.for.regs_$mean_pupil_nobc_scaled, rule = 2)
    df.pupil.for.regs_$mean_pupil_nobc_scaled_i <- NA
    df.pupil.for.regs_$mean_pupil_nobc_scaled_i[1:length(int)] <- int
    
    #lag pupil for arousal index
    df.pupil.for.regs_$mean_pupil_nobc_scaled_lag <- NA
    df.pupil.for.regs_$mean_pupil_nobc_scaled_lag <- c(
      lead(df.pupil.for.regs_$mean_pupil_nobc_scaled_i)[1:(nrow(df.pupil.for.regs_)-1)],
      df.pupil.for.regs_$mean_pupil_nobc_scaled_i[nrow(df.pupil.for.regs_)])
    
    int <- zoo::na.approx(df.pupil.for.regs_$pupil_x, rule = 2)
    df.pupil.for.regs_$pupil_x_i <- NA
    df.pupil.for.regs_$pupil_x_i[1:length(int)] <- int
    
    int <- zoo::na.approx(df.pupil.for.regs_$pupil_y, rule = 2)
    df.pupil.for.regs_$pupil_y_i <- NA
    df.pupil.for.regs_$pupil_y_i[1:length(int)] <- int
    
    int <- zoo::na.approx(df.pupil.for.regs_$blink, rule = 2)
    df.pupil.for.regs_$blink_i <- NA
    df.pupil.for.regs_$blink_i[1:length(int)] <- int
    
    int <- zoo::na.approx(df.pupil.for.regs_$sacc, rule = 2)
    df.pupil.for.regs_$sacc_i <- NA
    df.pupil.for.regs_$sacc_i[1:length(int)] <- int
    
    
    write_delim(data.frame(df.pupil.for.regs_$luminance),
                paste0("~/Desktop/midaffemo/analyses/pupil/glm_regs/",sub,"_luminance.1D"), col_names = F)
    write_delim(data.frame(df.pupil.for.regs_$pupil_x_i),
                paste0("~/Desktop/midaffemo/analyses/pupil/glm_regs/",sub,"_pupilx.1D"), col_names = F)
    write_delim(data.frame(df.pupil.for.regs_$pupil_y_i),
                paste0("~/Desktop/midaffemo/analyses/pupil/glm_regs/",sub,"_pupily.1D"), col_names = F)
    write_delim(data.frame(df.pupil.for.regs_$blink_i),
                paste0("~/Desktop/midaffemo/analyses/pupil/glm_regs/",sub,"_blink.1D"), col_names = F)
    write_delim(data.frame(df.pupil.for.regs_$sacc_i),
                paste0("~/Desktop/midaffemo/analyses/pupil/glm_regs/",sub,"_saccade.1D"), col_names = F)
    write_delim(data.frame(df.pupil.for.regs_$mean_pupil_nobc_scaled_i),
                paste0("~/Desktop/midaffemo/analyses/pupil/glm_regs/",sub,"_scaled.1D"), col_names = F)
    write_delim(data.frame(df.pupil.for.regs_$mean_pupil_nobc_scaled_lag),
                paste0("~/Desktop/midaffemo/analyses/pupil/glm_regs/",sub,"_scaled_lag.1D"), col_names = F)
  }
}



## for midaffemo2
fmri4_pupil_data<-read_csv("~/Desktop/midaffemo2/analyses/pupil/derivatives/pupillometry_baselineCorrected.csv")

fmri4_pupil_data.for.regs <- fmri4_pupil_data%>%
  mutate(tr_in_trial = floor((Time_sec+2)/2))  %>%
  filter(tr_in_trial > 0)%>%
  group_by(subject,block,trial,
           tr_in_trial)%>%
  summarise(pupil_x = mean(pupil_x,na.rm = T),
            pupil_y = mean(pupil_x,na.rm = T),
            blink = mean(blink,na.rm = T),
            sacc = mean(sacc,na.rm = T),
            mean_pupil_nobc_scaled= mean(pupilDiameter_scaled, na.rm = T))

df.1sub <- read_csv("~/Desktop/midaffemo2/data/timecourses/timecourses_b4_woOutliers_long.csv") %>%
  filter(subject == "ad241211", voi == "nacc",coreg == "ants") %>%
  mutate(trial_tr_noiti = ifelse(probe == 1, 9,8),
         trial_tr = trial_tr_noiti + iti/2) %>%
  rename(tr_in_trial = tr) %>% 
  select(trial,trial_tr,tr_in_trial,run) %>%
  filter(tr_in_trial <= trial_tr)

end_run1 <- df.1sub %>%
  subset(trial == 24) %>%
  .[1:5,] %>%
  mutate(tr_in_trial = tr_in_trial + 10)

# create the 5 TRs at the tail of run 2
end_run2 <- df.1sub %>%
  filter(run == "run-02")%>%
  subset(trial == 48) %>%
  .[1:5,] %>%
  mutate(tr_in_trial = tr_in_trial + 9)

df.1sub <- df.1sub %>%
  # insert 5 rows at the end of each run
  add_row(end_run1, .after = 240) %>%
  add_row(end_run2, .after = 485)

for (sub in unique(fmri4_pupil_data.for.regs$subject)){
  if(!is.na(sub)){
    df.pupil.for.regs_ <- left_join(df.1sub,fmri4_pupil_data.for.regs %>%
                                      filter(subject == sub),
                                    by = c("trial","tr_in_trial")) %>%
      relocate(subject)

    # linear interpolation; rule=2 extrapolates leading/trailing NAs with first/last non-NA
    int <- zoo::na.approx(df.pupil.for.regs_$mean_pupil_nobc_scaled, rule = 2)
    df.pupil.for.regs_$mean_pupil_nobc_scaled_i <- NA
    df.pupil.for.regs_$mean_pupil_nobc_scaled_i[1:length(int)] <- int

    #lag pupil for arousal index
    df.pupil.for.regs_$mean_pupil_nobc_scaled_lag <- NA
    df.pupil.for.regs_$mean_pupil_nobc_scaled_lag <- c(
      lead(df.pupil.for.regs_$mean_pupil_nobc_scaled_i)[1:(nrow(df.pupil.for.regs_)-1)],
      df.pupil.for.regs_$mean_pupil_nobc_scaled_i[nrow(df.pupil.for.regs_)])

    int <- zoo::na.approx(df.pupil.for.regs_$pupil_x, rule = 2)
    df.pupil.for.regs_$pupil_x_i <- NA
    df.pupil.for.regs_$pupil_x_i[1:length(int)] <- int

    int <- zoo::na.approx(df.pupil.for.regs_$pupil_y, rule = 2)
    df.pupil.for.regs_$pupil_y_i <- NA
    df.pupil.for.regs_$pupil_y_i[1:length(int)] <- int

    int <- zoo::na.approx(df.pupil.for.regs_$blink, rule = 2)
    df.pupil.for.regs_$blink_i <- NA
    df.pupil.for.regs_$blink_i[1:length(int)] <- int

    int <- zoo::na.approx(df.pupil.for.regs_$sacc, rule = 2)
    df.pupil.for.regs_$sacc_i <- NA
    df.pupil.for.regs_$sacc_i[1:length(int)] <- int

    write_delim(data.frame(df.pupil.for.regs_$pupil_x_i),
                paste0("~/Desktop/midaffemo2/analyses/pupil/glm_pupil_regs/",sub,"_pupilx.1D"), col_names = F)
    write_delim(data.frame(df.pupil.for.regs_$pupil_y_i),
                paste0("~/Desktop/midaffemo2/analyses/pupil/glm_pupil_regs/",sub,"_pupily.1D"), col_names = F)
    write_delim(data.frame(df.pupil.for.regs_$blink_i),
                paste0("~/Desktop/midaffemo2/analyses/pupil/glm_pupil_regs/",sub,"_blink.1D"), col_names = F)
    write_delim(data.frame(df.pupil.for.regs_$sacc_i),
                paste0("~/Desktop/midaffemo2/analyses/pupil/glm_pupil_regs/",sub,"_saccade.1D"), col_names = F)
    write_delim(data.frame(df.pupil.for.regs_$mean_pupil_nobc_scaled_i),
                paste0("~/Desktop/midaffemo2/analyses/pupil/glm_pupil_regs/",sub,"_scaled.1D"), col_names = F)
    write_delim(data.frame(df.pupil.for.regs_$mean_pupil_nobc_scaled_lag),
                paste0("~/Desktop/midaffemo2/analyses/pupil/glm_pupil_regs/",sub,"_scaled_lag.1D"), col_names = F)

  }
}
