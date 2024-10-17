library(dplyr)
library(tidyverse)
library(zoo)
library(chron)
library(lubridate)

bad_participants <- read_delim("~/Desktop/VRMID-analysis/mid-pupil/data/fmri3/subjects_list/subject_to_drop_pupil.txt", col_names = F, delim = " ")
# bad_participants <- c('av0312b1','el0312b3','el0312b4','js0311b3','js0311b4','jy0308b2','sk0301b3','sk0301b4')

subjects <- as.data.frame(read_csv("~/Desktop/VRMID-analysis/mid-pupil/data/fmri3/subjects_list/subjects-midaffemo.txt", col_names = F))

subjects_midaffemo <- paste0(rep(subjects$X1,each = 4),"_",rep(c("b1","b2","b3","b4"),times = nrow(subjects)))

subjects <- subjects_midaffemo
df.beh <- {}
for (sub in subjects){
  df.beh.temp <- read_csv(paste0("~/Desktop/VRMID-analysis/mid-pupil/data/fmri3/raw/behavior/",sub,".csv"))%>%
    filter(TR == 1)
  if (sum(is.na(df.beh.temp$arating)) > 22){ #for subj that does not rate enough, do not calculate scale score
    df.beh.temp <- df.beh.temp %>% 
      mutate(arousal_scaled = as.numeric(scale(arating)),
             valence_scaled = as.numeric(scale(vrating)))%>%
      ungroup()
  } else {
    df.beh.temp <- df.beh.temp %>% 
      mutate(arousal_scaled = NA,
             valence_scaled = NA)%>%
      ungroup()
  }
  df.beh.temp <- df.beh.temp%>% 
    select(-TR)%>% 
    mutate(last_iti = lag(iti)) %>% 
    mutate(subject = strsplit(sub,"_")[[1]][1],
           block = strsplit(sub,"_")[[1]][2])
  
  if ("trial_tr" %in% names(df.beh.temp)){
    df.beh.temp <- df.beh.temp%>% 
      select(-trial_tr)
  }
  df.beh <- rbind(df.beh,df.beh.temp)
}

df.beh <- df.beh %>% 
  mutate(trial = as.numeric(trial)) %>% 
  rowwise() %>% 
  mutate(trial = ifelse(subject %in% c("cw240110","lr240201"),
                        (as.numeric(strsplit(block, "b")[[1]][2]) - 1) * 48 + trial, 
                        (as.numeric(strsplit(block, "b")[[1]][2]) - 1) * 24 + trial))%>% 
  relocate(subject,block,trial)

write_csv(df.beh,"../behavior/derivatives/beh.csv")

sr = 200
# data wo baseline correction
data_pupil <- read_csv("../pupil/derivatives/pupillometry.csv")

subjs_temp <- unique(data_pupil$subject)
subjs_temp
data_pupil$block <- NA
data_pupil$sub_id <- data_pupil$subject
for (sub in subjs_temp){
  print(sub)
  if (grepl("b",sub)){
    rowids <- which(data_pupil$subject == sub)
    # my.str <- strsplit(data_pupil$subject[rowids],"b")[[1]] #this doesn't work when participants' initials contain b!!!
    my.str1 <- substring(data_pupil$subject[rowids],1,6)
    my.str2 <- substring(data_pupil$subject[rowids],7,9)
    data_pupil$subject[rowids] <- my.str1
    data_pupil$subject[rowids] <- paste0(substring(my.str1,1,2),
                                         '24',substring(my.str1,3,6))
    data_pupil$block[rowids] <- my.str2
  }
}

unique(data_pupil$sub_id)

data_pupil <- data_pupil %>% 
  filter(!sub_id %in% bad_participants)

unique(data_pupil$subject)

unique(data_pupil$block)

length(unique(data_pupil$subject))

length(unique(data_pupil$sub_id))

data_pupil_ <- data_pupil%>%
  group_by(subject) %>% 
  mutate(times = (row_number()-1)/sr) %>% 
  filter(!subject %in% bad_participants,
         !is.na(subject))%>%
  ungroup() %>% 
  group_by(subject,block,trial)%>%
  mutate(sample_in_trial_n = row_number(),
         sample_in_trial_t = (sample_in_trial_n-1)/sr,
         Time_sec = floor(sample_in_trial_t))%>% 
  rowwise() %>% 
  mutate(trial = ifelse(subject %in% c("cw240110","lr240201"),
                        (as.numeric(strsplit(block, "b")[[1]][2]) - 1) * 48 + trial, 
                        (as.numeric(strsplit(block, "b")[[1]][2]) - 1) * 24 + trial))

data_pupil_ <- data_pupil_%>% 
  left_join(df.beh)

a = data_pupil_ %>% filter(subject == "pm240301")
max(a$sample_in_trial_t)

names(data_pupil_)

View(data_pupil_%>%
  tail(100))

sort(unique(data_pupil_$Time_sec)) #should be 0 to 21 (22 is rare)

data_pupil_$trial_duration = 16+as.numeric(data_pupil_$iti)

sort(unique(data_pupil_$trial_duration)) #16, 18, 20, 22

sort(unique(data_pupil_$Time_sec)) #16, 18, 20, 22

View(count(data_pupil_ %>% ungroup(),Time_sec))
data_pupil_$Time_sec[1] #should be 0

#baseline stuff
data_pupil_out <- {}
my.subjects <- unique(data_pupil_$subject)
my.subjects
# get the last two seconds of previous ITI and paste it to this trial
# for (i in 1:length(unique(data_pupil_$subject))){
for (sub in my.subjects){
  print(paste0("processing subject ", sub))
#   sub = unique(data_pupil_$subject)[i]
  if (!sub %in% bad_participants){
    temp_data <- data_pupil_%>%
      filter(subject == sub)%>%
      ungroup()
      
        for (j in 1:length(unique(temp_data$trial))){
          temp_trial_data <- temp_data%>%
            filter(trial == j)
          
          if (j == 1){
            d <- temp_trial_data
          } else{
            
            last_trial_iti_data <- temp_data%>%
              group_by(trial) %>% 
              mutate(max_time = max(Time_sec)) %>% 
              filter(trial == (j-1),
                     Time_sec == trial_duration-1) #time_sec start from 0 so should be t-1
            
            #if last second exists
            if (nrow(last_trial_iti_data) != 0){ 
              
              # change timing info
              last_trial_iti_data <- last_trial_iti_data%>%
                mutate(times = times - max(last_trial_iti_data$times),
                       Time_sec = Time_sec - max(last_trial_iti_data$Time_sec) - 1,
                       sample_in_trial_n = sample_in_trial_n - max(last_trial_iti_data$sample_in_trial_n),
                       sample_in_trial_t = sample_in_trial_t - max(last_trial_iti_data$sample_in_trial_t),
                       event = "baseline") %>% 
                ungroup() %>% 
                select(blink:Time_sec)
              
              d <- merge(temp_trial_data,last_trial_iti_data,
                         all.x = T, all.y = T)%>%
                arrange(sample_in_trial_n)%>%
                fill(subject,block,trial,trialonset:valence_scaled, .direction = "downup") %>% 
                relocate(trial)
            }
            else {
              d <- temp_trial_data
            }
          }
          
          data_pupil_out <- rbind(data_pupil_out,d)
        }
     #trial loop
  } #is bad subj?
}#subject loop


unique(data_pupil_out$subject)

#get the mean value of the last second in a trial (ITI)
baseline <- data_pupil_out%>%
  group_by(subject,block,trial)%>%
  filter(Time_sec == -1)%>%
  summarise_at(vars(pupilDiameter), ~mean(.x, na.rm = TRUE))%>%
  ungroup()

data_pupil_out$pupilDiameter_baseline <- NA

for (i in 1:nrow(baseline)){
  ind = which(data_pupil_out$subject == baseline$subject[i] &
                data_pupil_out$block == baseline$block[i] &
                data_pupil_out$trial == baseline$trial[i])
  #baseline info
  if (!is_empty(ind)){
    #baseline is from last trial's ITI
    data_pupil_out$pupilDiameter_baseline[ind] = baseline$pupilDiameter[i]
  }
}

#baseline correction; first trial of each participant was excluded
data_pupil_out <- data_pupil_out%>%
  mutate(pupilDiameter_bc = pupilDiameter - pupilDiameter_baseline)

write_csv(data_pupil_out, "~/Desktop/VRMID-analysis/mid-pupil/data/fmri3/derivatives/pupillometry_baselineCorrected.csv")
