library(dplyr)
library(tidyverse)
library(zoo)
library(chron)
library(lubridate)


bad_participants <- read.csv("~/Desktop/VRMID-analysis/mid-pupil/data/vrmid2/subjects_list/subject_to_drop_pupil.txt",
                             header = F)

sr = 120
# data w baseline correction
data_pupil <- read_csv("~/Desktop/VRMID-analysis/mid-pupil/data/vrmid2/derivatives/pupillometry.csv")%>%
  filter(!Subject %in% bad_participants,
         !is.na(Subject))%>%
  group_by(Subject, Time)%>%
  mutate(sample_in_sec = row_number())%>%
  ungroup()%>%
  group_by(Subject,trial)%>%
  mutate(sample_in_trial_n = row_number(),
         sample_in_trial_t = sample_in_trial_n/sr)

length(unique(data_pupil$Subject))
unique(data_pupil$emotion)

data_pupil%>%
  head(100)

# data_pupil$trial_type <- factor(data_pupil$trial_type,
#                                 levels = c("-$5","-$1","-$0",
#                                            "+$0","+$1","+$5"))

data_pupil <- data_pupil%>%
  filter(!is.na(Subject))

sort(unique(data_pupil$Time_sec))

data_pupil$Time_sec[1]
#baseline stuff
data_pupil_out <- {}
# get the last two seconds of previous ITI and paste it to this trial
for (i in 1:length(unique(data_pupil$Subject))){
  sub = unique(data_pupil$Subject)[i]
  if (!sub %in% bad_participants){
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
                 Time_sec == 18 + ITI_Duration)
        
        #if last second exists
        if (nrow(last_trial_iti_data) != 0){ 
          # last 1 second
          last_trial_iti_data <- last_trial_iti_data%>%
            filter(Time_in_trial_sec == max(last_trial_iti_data$Time_in_trial_sec))
          
          # change timing info
          last_trial_iti_data <- last_trial_iti_data%>%
            mutate(Time_sec = Time_sec -max(last_trial_iti_data$Time_sec),
                   Time_in_trial_sec = Time_in_trial_sec - max(last_trial_iti_data$Time_in_trial_sec),
                   sample_in_trial_n = sample_in_trial_n - max(last_trial_iti_data$sample_in_trial_n),
                   sample_in_trial_t = sample_in_trial_t - max(last_trial_iti_data$sample_in_trial_t),
                   current_stimulus = "prestim_baseline"
                   )%>%
            select(Subject,Time,Time_sec,Time_str,Time_in_trial_sec,current_stimulus,sample_in_sec,sample_in_trial_n:sample_in_trial_t,pupil_L:pupil_Avg)
          
          d <- merge(temp_trial_data,last_trial_iti_data,
                     all.x = T, all.y = T)%>%
            arrange(sample_in_trial_n)%>%
            fill(trial,trial_type:NegA_scaled,emotion, .direction = "downup")%>%
            relocate(trial,.after = Subject)
        }
        else {
          d <- temp_trial_data
        }
      }
      
      data_pupil_out <- rbind(data_pupil_out,d)
    } #trial loop
  }
}#subject loop

#get the mean value of the last second in a trial (ITI)
baseline <- data_pupil_out%>%
  group_by(Subject, trial)%>%
  filter(Time_sec == 0)%>%
  summarise_at(vars(pupil_L,pupil_R,pupil_Avg), ~mean(.x, na.rm = TRUE))%>%
  ungroup()

data_pupil_out$pupil_L_baseline <- NA
data_pupil_out$pupil_R_baseline <- NA
data_pupil_out$pupil_Avg_baseline <- NA

for (i in 1:nrow(baseline)){
  ind = which(data_pupil_out$Subject == baseline$Subject[i] &
                data_pupil_out$trial == baseline$trial[i])
  #baseline info
  if (!is_empty(ind)){
    #baseline is from last trial's ITI
    data_pupil_out$pupil_L_baseline[ind] = baseline$pupil_L[i]
    data_pupil_out$pupil_R_baseline[ind] = baseline$pupil_R[i]
    data_pupil_out$pupil_Avg_baseline[ind] = baseline$pupil_Avg[i]
  }
}

#baseline correction; first trial of each participant was excluded
data_pupil_out <- data_pupil_out%>%
  mutate(pupil_L_bc = pupil_L - pupil_L_baseline,
         pupil_R_bc = pupil_R - pupil_R_baseline,
         pupil_Avg_bc = pupil_Avg - pupil_Avg_baseline)

write_csv(data_pupil_out, "~/Desktop/VRMID-analysis/mid-pupil/data/vrmid2/derivatives/pupillometry_baselineCorrected.csv")

# remove(data_pupil)
