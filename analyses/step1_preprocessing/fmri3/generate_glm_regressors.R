# get luminance time course for regression (downsample to TR = 2 sec)

fmri3_pupil_data.for.regs <- fmri3_pupil_data%>%
  mutate(tr = floor((Time_sec+2)/2))  %>% 
  filter(tr > 0)%>% 
  group_by(subject,block,trial,probe,rating_type,
           tr)%>%
  summarise(luminance= mean(luminance, na.rm = T),
            pupil_x = mean(pupil_x,na.rm = T),
            pupil_y = mean(pupil_x,na.rm = T),
            blink = mean(blink,na.rm = T),
            sacc = mean(sacc,na.rm = T),
            mean_pupil_nobc_scaled= mean(pupilDiameter_scaled, na.rm = T))

df.1sub <- read_csv("~/Desktop/midaffemo/analyses/timecourses/data/realsubj/MIDaffemo/timecourses_b0_woOutliers_long.csv") %>% 
  filter(subject == "nq240330", voi == "nacc",coreg == "ants") %>% 
  mutate(trial_tr_noiti = ifelse(probe == 1, 9,8),
         trial_tr = trial_tr_noiti + iti/2) %>% 
  select(trial,trial_tr,tr,run) %>% 
  filter(tr <= trial_tr)

end_run1 <- df.1sub %>%
  subset(trial == 24) %>%
  .[1:5,] %>%
  mutate(tr = tr + 10)

# create the 5 TRs at the tail of run 2
end_run2 <- df.1sub %>%
  filter(run == "run-02")%>%
  subset(trial == 24) %>%
  .[1:5,] %>%
  mutate(tr = tr + 9)

# create the 5 TRs at the tail of run 3
end_run3 <- df.1sub %>%
  filter(run == "run-03")%>%
  subset(trial == 24) %>%
  .[1:5,] %>%
  mutate(tr = tr + 10)

# create the 5 TRs at the tail of run 4
end_run4 <- df.1sub %>%
  filter(run == "run-04")%>%
  subset(trial == 24) %>%
  .[1:5,] %>%
  mutate(tr = tr + 11)

df.1sub <- df.1sub %>%
  # insert 5 rows at the end of each run
  add_row(end_run1, .after = 240) %>%
  add_row(end_run2, .after = 485) %>%
  add_row(end_run3, .after = 730) %>%
  add_row(end_run4, .after = 975)

df.1sub <- df.1sub %>% 
  mutate(block = as.numeric(substr(run,6,6))) %>% 
  mutate(trial = trial + (block - 1) * 24) %>% 
  select(-block,-run)

for (sub in unique(fmri3_pupil_data.for.regs$subject)){
  if(!is.na(sub)){
    df.pupil.for.regs_ <- left_join(df.1sub,fmri3_pupil_data.for.regs %>% 
                                      filter(subject == sub),
                                    by = c("trial","tr")) %>% 
      relocate(subject)
    
    df.pupil.for.regs_$luminance[is.na(df.pupil.for.regs_$luminance)]  <- 0.002513428#insert approximate luminance for this screen
    
    # # #linear interpolation
    int <- zoo::na.approx(df.pupil.for.regs_$mean_pupil_nobc_scaled)
    df.pupil.for.regs_$mean_pupil_nobc_scaled_i <- NA
    df.pupil.for.regs_$mean_pupil_nobc_scaled_i[1:length(int)] <-int
    df.pupil.for.regs_$mean_pupil_nobc_scaled_i[is.na(df.pupil.for.regs_$mean_pupil_nobc_scaled_i)] <- 0
    
    int <- zoo::na.approx(df.pupil.for.regs_$pupil_x)
    df.pupil.for.regs_$pupil_x_i <- NA
    df.pupil.for.regs_$pupil_x_i[1:length(int)] <-int
    df.pupil.for.regs_$pupil_x_i[is.na(df.pupil.for.regs_$pupil_x_i)] <- 0
    
    int <- zoo::na.approx(df.pupil.for.regs_$pupil_y)
    df.pupil.for.regs_$pupil_y_i <- NA
    df.pupil.for.regs_$pupil_y_i[1:length(int)] <-int
    df.pupil.for.regs_$pupil_y_i[is.na(df.pupil.for.regs_$pupil_y_i)] <- 0
    
    int <- zoo::na.approx(df.pupil.for.regs_$blink)
    df.pupil.for.regs_$blink_i <- NA
    df.pupil.for.regs_$blink_i[1:length(int)] <-int
    df.pupil.for.regs_$blink_i[is.na(df.pupil.for.regs_$blink_i)] <- 0
    
    int <- zoo::na.approx(df.pupil.for.regs_$sacc)
    df.pupil.for.regs_$sacc_i <- NA
    df.pupil.for.regs_$sacc_i[1:length(int)] <-int
    df.pupil.for.regs_$sacc_i[is.na(df.pupil.for.regs_$sacc_i)] <- 0
    
    
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
  }
}


## for midaffemo2
fmri4_pupil_data.for.regs <- fmri4_pupil_data%>%
  mutate(tr = floor((Time_sec+2)/2))  %>% 
  filter(tr > 0)%>% 
  group_by(subject,block,trial,
           tr)%>%
  summarise(pupil_x = mean(pupil_x,na.rm = T),
            pupil_y = mean(pupil_x,na.rm = T),
            blink = round(mean(blink,na.rm = T)),
            sacc = mean(sacc,na.rm = T),
            mean_pupil_nobc_scaled= mean(pupilDiameter_scaled, na.rm = T))

df.1sub <- read_csv("~/Desktop/midaffemo2/data/timecourses/timecourses_b4_woOutliers_long.csv") %>% 
  filter(subject == "ad241211", voi == "nacc",coreg == "ants") %>% 
  mutate(trial_tr_noiti = ifelse(probe == 1, 9,8),
         trial_tr = trial_tr_noiti + iti/2) %>% 
  select(trial,trial_tr,tr,run) %>% 
  filter(tr <= trial_tr)

end_run1 <- df.1sub %>%
  subset(trial == 24) %>%
  .[1:5,] %>%
  mutate(tr = tr + 10)

# create the 5 TRs at the tail of run 2
end_run2 <- df.1sub %>%
  filter(run == "run-02")%>%
  subset(trial == 48) %>%
  .[1:5,] %>%
  mutate(tr = tr + 9)

df.1sub <- df.1sub %>%
  # insert 5 rows at the end of each run
  add_row(end_run1, .after = 240) %>%
  add_row(end_run2, .after = 485)

for (sub in unique(fmri4_pupil_data.for.regs$subject)){
  if(!is.na(sub)){
    df.pupil.for.regs_ <- left_join(df.1sub,fmri4_pupil_data.for.regs %>% 
                                      filter(subject == sub),
                                    by = c("trial","tr")) %>% 
      relocate(subject)
  
    # # #linear interpolation
    int <- zoo::na.approx(df.pupil.for.regs_$mean_pupil_nobc_scaled)
    df.pupil.for.regs_$mean_pupil_nobc_scaled_i <- NA
    df.pupil.for.regs_$mean_pupil_nobc_scaled_i[1:length(int)] <-int
    df.pupil.for.regs_$mean_pupil_nobc_scaled_i[is.na(df.pupil.for.regs_$mean_pupil_nobc_scaled_i)] <- 0
    
    int <- zoo::na.approx(df.pupil.for.regs_$pupil_x)
    df.pupil.for.regs_$pupil_x_i <- NA
    df.pupil.for.regs_$pupil_x_i[1:length(int)] <-int
    df.pupil.for.regs_$pupil_x_i[is.na(df.pupil.for.regs_$pupil_x_i)] <- 0
    
    int <- zoo::na.approx(df.pupil.for.regs_$pupil_y)
    df.pupil.for.regs_$pupil_y_i <- NA
    df.pupil.for.regs_$pupil_y_i[1:length(int)] <-int
    df.pupil.for.regs_$pupil_y_i[is.na(df.pupil.for.regs_$pupil_y_i)] <- 0
    
    int <- zoo::na.approx(df.pupil.for.regs_$blink)
    df.pupil.for.regs_$blink_i <- NA
    df.pupil.for.regs_$blink_i[1:length(int)] <-int
    df.pupil.for.regs_$blink_i[is.na(df.pupil.for.regs_$blink_i)] <- 0
    
    int <- zoo::na.approx(df.pupil.for.regs_$sacc)
    df.pupil.for.regs_$sacc_i <- NA
    df.pupil.for.regs_$sacc_i[1:length(int)] <-int
    df.pupil.for.regs_$sacc_i[is.na(df.pupil.for.regs_$sacc_i)] <- 0
    
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
    
  }
}
