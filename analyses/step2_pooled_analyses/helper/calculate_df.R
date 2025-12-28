
## calculate df

vrmid1_data_antprobe <-  vrmid1_pupil_data%>%
  mutate(pupil_window = ifelse(sample_in_trial_t <= vrmid1_antprobe_out_window[2] & sample_in_trial_t > vrmid1_antprobe_out_window[1],
                               "out",
                               ifelse(sample_in_trial_t <= vrmid1_ant_window[2] & sample_in_trial_t > vrmid1_ant_window[1],
                                      "ant",NA)))%>%
  filter(probe == "anti", pupil_window %in% c("out","ant"))%>%
  ungroup()%>%
  group_by(probe,Subject,trial,reaction_outcome,trial_type,PosA_scaled,NegA_scaled,
           arousal,valence,pupil_window,last_ITI,
           arousal_scaled,valence_scaled,
           condition,reaction_time)%>%
  summarise(mean_pupil_bc = mean(pupil_Avg_bc_scaled, na.rm = T),
            mean_pupil_nobc = mean(pupil_Avg_scaled, na.rm = T),
            mean_pupil_baseline = mean(pupil_Avg_baseline_scaled, na.rm = T))%>%
  filter(abs(mean_pupil_bc) < 4,
         abs(mean_pupil_nobc) < 4,
         abs(mean_pupil_baseline) < 4)  %>% 
  pivot_wider(values_from = mean_pupil_bc:mean_pupil_baseline, names_from = pupil_window)%>%
  ungroup()%>%
  group_by(Subject)%>%
  mutate(log_rt = log(reaction_time),
         log_rt = as.numeric(scale(log_rt))) %>% 
  rowwise() %>% 
  mutate(trial_value = as.numeric(paste0(strsplit(as.character(trial_type), "\\$")[[1]][1],strsplit(as.character(trial_type), "\\$")[[1]][2]))) %>% 
  mutate(trial_outcome = ifelse(reaction_outcome == "Hit" & trial_value >= 0, trial_value,
                                ifelse(reaction_outcome == "Miss" & trial_value >= 0, 0,
                                       ifelse(reaction_outcome == "Hit" & trial_value < 0, 0,
                                              ifelse(reaction_outcome == "Miss" & trial_value < 0,trial_value))))) %>% 
  ungroup()

vrmid1_data_outprobe <-  vrmid1_pupil_data%>%
  mutate(pupil_window = ifelse(probe == "out" & sample_in_trial_t <= vrmid1_outprobe_out_window[2] & 
                                 sample_in_trial_t > vrmid1_outprobe_out_window[1],
                               "out",
                               ifelse(sample_in_trial_t <= vrmid1_ant_window[2] & 
                                        sample_in_trial_t > vrmid1_ant_window[1],
                                      "ant",NA)))%>%
  filter(probe == "out",pupil_window %in% c("out","ant"))%>% #probe == "out"
  ungroup()%>%
  group_by(probe,Subject,trial,trial_type,PosA_scaled,NegA_scaled,
           arousal,valence,pupil_window,last_ITI,
           arousal_scaled,valence_scaled,
           condition,reaction_outcome,reaction_time)%>%
  summarise(mean_pupil_bc = mean(pupil_Avg_bc_scaled, na.rm = T),
            mean_pupil_nobc = mean(pupil_Avg_scaled, na.rm = T),
            mean_pupil_baseline = mean(pupil_Avg_baseline_scaled, na.rm = T))%>%
  filter(abs(mean_pupil_bc) < 4,#this would remove trial 1 due to no baseline
         abs(mean_pupil_nobc) < 4, 
         abs(mean_pupil_baseline) < 4)  %>% 
  pivot_wider(values_from = mean_pupil_bc:mean_pupil_baseline, names_from = pupil_window)%>%
  ungroup()%>%
  group_by(Subject)%>%
  mutate(log_rt = log(reaction_time),
         log_rt = as.numeric(scale(log_rt)))%>% 
  rowwise() %>% 
  mutate(trial_value = as.numeric(paste0(strsplit(as.character(trial_type), "\\$")[[1]][1],strsplit(as.character(trial_type), "\\$")[[1]][2]))) %>% 
  mutate(trial_outcome = ifelse(reaction_outcome == "Hit" & trial_value >= 0, trial_value,
                                ifelse(reaction_outcome == "Miss" & trial_value >= 0, 0,
                                       ifelse(reaction_outcome == "Hit" & trial_value < 0, 0,
                                              ifelse(reaction_outcome == "Miss" & trial_value < 0,trial_value))))) %>% 
  ungroup()


fmri3_data_antprobe <-  fmri3_pupil_data%>%
  filter(!subject %in% c("tc240311","st240312")) %>% 
  mutate(pupil_window = ifelse(sample_in_trial_t <= fmri3_antprobe_out_window[2] & 
                                 sample_in_trial_t > fmri3_antprobe_out_window[1],
                               "out",
                               ifelse(sample_in_trial_t <= fmri3_ant_window[2] & 
                                        sample_in_trial_t > fmri3_ant_window[1],
                                      "ant",NA)))%>%
  filter(probe == 1, pupil_window %in% c("out","ant"))%>%
  ungroup()%>%
  group_by(subject) %>% 
  mutate(arousal_scaled = as.numeric(scale(arousal)),
         valence_scaled = as.numeric(scale(valence)))%>% 
  group_by(probe,subject,trial,hit,cue_value,emotion,emotion_cat,
           arousal,valence,pupil_window,rating_type,
           arousal_scaled,valence_scaled,cue_size,rt,iti,last_iti)%>%
  summarise(mean_pupil_bc_scaled = mean(pupilDiameter_bc_scaled, na.rm = T),
            mean_pupil_nobc_scaled= mean(pupilDiameter_scaled, na.rm = T),
            mean_pupil_baseline_scaled = mean(pupilDiameter_baseline_scaled, na.rm = T),
            mean_pupil_nobc_scaled_residual_conv = mean(residuals_conv,na.rm = T),
            luminance_mean = mean(luminance,na.rm = T),
            luminance_lag1_mean = mean(luminance_lag1,na.rm = T),
            luminance_convolved_mean = mean(convolved_s1,na.rm = T))%>%
  pivot_wider(values_from = c(mean_pupil_bc_scaled,mean_pupil_nobc_scaled,mean_pupil_nobc_scaled_residual_conv,luminance_mean,luminance_lag1_mean,luminance_convolved_mean), names_from = pupil_window)%>%
  ungroup()%>%
  group_by(subject)%>%
  mutate(log_rt = log(rt),
         log_rt = as.numeric(scale(log_rt)))%>% 
  rowwise() %>% 
  mutate(trial_value = as.numeric(paste0(strsplit(as.character(cue_value), "\\$")[[1]][1],strsplit(as.character(cue_value), "\\$")[[1]][2]))) %>% 
  mutate(trial_outcome = ifelse(hit == "Hit" & trial_value >= 0, trial_value,
                                ifelse(hit == "Miss" & trial_value >= 0, 0,
                                       ifelse(hit == "Hit" & trial_value < 0, 0,
                                              ifelse(hit == "Miss" & trial_value < 0,trial_value,NA))))) %>% 
  ungroup()

fmri3_data_outprobe <-  fmri3_pupil_data%>%
  filter(!subject %in% c("tc240311","st240312")) %>% 
  mutate(pupil_window = ifelse(probe == 2 & sample_in_trial_t <= fmri3_outprobe_out_window[2] 
                               & sample_in_trial_t > fmri3_outprobe_out_window[1],
                               "out",
                               ifelse(probe == 2 &  sample_in_trial_t <= fmri3_ant_window[2] 
                                      & sample_in_trial_t > fmri3_ant_window[1],
                                      "ant",NA)))%>%
  filter(probe == 2, pupil_window %in% c("out","ant"))%>%
  ungroup()%>%
  group_by(subject) %>% 
  mutate(arousal_scaled = as.numeric(scale(arousal)),
         valence_scaled = as.numeric(scale(valence)))%>% 
  group_by(probe,subject,trial,hit,cue_value,emotion,emotion_cat,
           arousal,valence,pupil_window,rating_type,
           arousal_scaled,valence_scaled,cue_size,rt,iti,last_iti)%>%
  summarise(mean_pupil_bc_scaled = mean(pupilDiameter_bc_scaled, na.rm = T),
            mean_pupil_nobc_scaled= mean(pupilDiameter_scaled, na.rm = T),
            mean_pupil_baseline_scaled = mean(pupilDiameter_baseline_scaled, na.rm = T),
            mean_pupil_nobc_scaled_residual_conv = mean(residuals_conv,na.rm = T),
            luminance_mean = mean(luminance,na.rm = T),
            luminance_lag1_mean = mean(luminance_lag1,na.rm = T),
            luminance_convolved_mean = mean(convolved_s1,na.rm = T))%>%
  pivot_wider(values_from = c(mean_pupil_bc_scaled,mean_pupil_nobc_scaled,mean_pupil_nobc_scaled_residual_conv,luminance_mean,luminance_lag1_mean,luminance_convolved_mean), names_from = pupil_window)%>%
  ungroup()%>%
  ungroup()%>%
  group_by(subject)%>%
  mutate(log_rt = log(rt),
         log_rt = as.numeric(scale(log_rt)))%>% 
  rowwise() %>% 
  mutate(trial_value = as.numeric(paste0(strsplit(as.character(cue_value), "\\$")[[1]][1],strsplit(as.character(cue_value), "\\$")[[1]][2]))) %>% 
  mutate(trial_outcome = ifelse(hit == "Hit" & trial_value >= 0, trial_value,
                                ifelse(hit == "Miss" & trial_value >= 0, 0,
                                       ifelse(hit == "Hit" & trial_value < 0, 0,
                                              ifelse(hit == "Miss" & trial_value < 0,trial_value,NA))))) %>%
  ungroup()



## combine two probe types
vrmid1_data_2probes <- rbind(vrmid1_data_antprobe %>% 
                               select(probe:reaction_time,arousal,valence,trial_value,trial_outcome,
                                      mean_pupil_nobc_ant,mean_pupil_bc_ant,mean_pupil_baseline_ant,last_ITI) %>% 
                               rename(mean_pupil_nobc=mean_pupil_nobc_ant,
                                      mean_pupil_bc = mean_pupil_bc_ant,
                                      mean_pupil_baseline = mean_pupil_baseline_ant),
                             vrmid1_data_outprobe %>% 
                               select(probe:reaction_time,arousal,valence,trial_value,trial_outcome,mean_pupil_nobc_out,mean_pupil_bc_out,
                                      mean_pupil_baseline_ant,last_ITI) %>% 
                               rename(mean_pupil_nobc=mean_pupil_nobc_out,
                                      mean_pupil_bc = mean_pupil_bc_out,
                                      mean_pupil_baseline = mean_pupil_baseline_ant))

fmri3_data_2probes <- rbind(fmri3_data_antprobe %>% 
                              select(probe:rt,trial_value,trial_outcome,
                                     mean_pupil_nobc_scaled_ant,mean_pupil_bc_scaled_ant,mean_pupil_baseline_scaled,last_iti) %>% 
                              rename(mean_pupil_nobc=mean_pupil_nobc_scaled_ant,
                                     mean_pupil_bc = mean_pupil_bc_scaled_ant,
                                     mean_pupil_baseline = mean_pupil_baseline_scaled),
                            fmri3_data_outprobe %>% 
                              select(probe:rt,trial_value,trial_outcome,
                                     mean_pupil_nobc_scaled_out,mean_pupil_bc_scaled_out,mean_pupil_baseline_scaled,last_iti) %>% 
                              rename(mean_pupil_nobc=mean_pupil_nobc_scaled_out,
                                     mean_pupil_bc = mean_pupil_bc_scaled_out,
                                     mean_pupil_baseline = mean_pupil_baseline_scaled))

write_csv(vrmid1_data_2probes,"../../data/vrmid1/derivatives/vrmid1_data_2probes.csv")
write_csv(fmri3_data_2probes,"../../data/fmri3/derivatives/fmri3_data_2probes.csv")
