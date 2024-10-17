library(dplyr)
library(tidyverse)
library(zoo)
library(chron)
library(lubridate)

#physio <- read.csv(paste0("./SENSORYDATASHEET_VRMID2.csv"))%>%
physio <- read_csv(paste0(getwd(),"/../data/study2/SENSORYDATASHEET_VRMID2.csv"))%>%
  filter(!is.na(LeftPDil))

#my.data <- read.csv("./BEHAVDATASHEET_VRMID2.csv")%>%
my.data <- read.csv(paste0(getwd(),"/../data/study2/BEHAVDATASHEET_VRMID2.csv"))%>%
  mutate(probe = ifelse(probe == "anti","out", 
                        ifelse(probe == "out","anti",NA)))%>%
  filter(subject %in% unique(physio$Subject))%>%
  rename(start_time = `start.time`,
         Subject = subject)%>%
  group_by(Subject)%>%
  mutate(start_time = times(start_time),
    end_time = times(lead(start_time,1)))%>%
  relocate(end_time, .after = start_time)

my.data$valence <- as.numeric(my.data$valence)
my.data$arousal <- as.numeric(my.data$arousal)

#converted valence and arousal to numeric
my.data <- my.data %>% group_by(Subject) %>% 
  mutate(mean_valence=mean(valence, na.rm=TRUE),
         mean_arousal=mean(arousal, na.rm=TRUE))

#the timing bit
physio$Time_str = times(physio$Time)

#some time points are wrong (from 12 pm to 1 am)
physio$Time_str <- times(ifelse(physio$Time_str < "7:00:00", times(physio$Time_str + "12:00:00"),
                                physio$Time_str))

joint_data <- {}
for (i in 1:length(unique(my.data$Subject))){
  print(i)
  curr_id = unique(my.data$Subject)[i]
  curr_data <- my.data[my.data$Subject == curr_id,]
  for (j in 1:nrow(curr_data)){
    nrowdata = curr_data[j,]
    nr = which(as.character(physio$Subject) == curr_id & 
                 times(physio$Time) >= times(nrowdata$start_time) &
                 times(physio$Time) < times(nrowdata$end_time))
    if (!is_empty(nr)){
      joint_data <- rbind(joint_data,
                          left_join(physio[nr,],nrowdata, by = "Subject"))
    }
  }
}
# data_pupil$Time_sec <- seconds(data_pupil$Time_in_trial_sec)
joint_data <- joint_data%>%
  group_by(Subject,trial)%>%
  mutate(Time_in_trial_sec = times(Time_str - start_time),
         Time_sec = as.numeric(seconds(Time_in_trial_sec)))%>%
  relocate(c(Time_str,Time_in_trial_sec,Time_sec),.after = Time)

joint_data <- joint_data%>%
  relocate(Time_sec,.after = Time)%>%
  ungroup()%>%
  mutate(current_stimulus = ifelse(Time_sec <= 4, "cue",
                                   ifelse(Time_sec %in% c(5,6),"fixation",
                                          ifelse((Time_sec %in% c(7,8) & probe == "out") |
                                                   (Time_sec %in% c(15,16) & probe == "anti"),"target",
                                                 ifelse((Time_sec %in% c(9,10) & probe == "out") |
                                                          (Time_sec %in% c(17,18) & probe == "anti"),"outcome",
                                                        ifelse((Time_sec %in% c(11:18) & probe == "out") |
                                                                 (Time_sec %in% c(7:14) & probe == "anti"),"self_report","ITI"))))),
         last_ITI = lag(ITI.Duration))%>%
  relocate(current_stimulus, .after = Time_sec)

joint_data <- joint_data %>% 
  ungroup()%>%
  mutate(sub_valence=valence - mean_valence,
         sub_arousal=arousal - mean_arousal)%>% 
  mutate(PosA = (sub_arousal+sub_valence)/sqrt(2), 
         NegA = (sub_arousal-sub_valence)/sqrt(2))%>%
  group_by(Subject)%>%
  mutate(PosA_scaled = as.numeric(scale(PosA)),
         NegA_scaled = as.numeric(scale(NegA)),
         arousal_scaled = as.numeric(scale(arousal)),
         valence_scaled = as.numeric(scale(valence)))

joint_data$condition <- factor(joint_data$condition, levels = c("large","receding","looming","small"))

joint_data$trial.type <- factor(joint_data$trial.type, levels = c("minus5","minus1","minus0","plus0","plus1","plus5"))

#get average pupil size (NA = closed)
joint_data <- joint_data%>%
  mutate(pupil_size = ifelse(LeftOpen == 1 & RightOpen == 1 & RightPDil != -1 &  LeftPDil != -1 &
                               RightPPos_x != -1 & LeftPPos_x != -1,
                             (RightPDil+LeftPDil)/2, NA))%>%
  relocate(Time_sec)%>%
  group_by(Subject)%>%
  mutate(pupil_size_scaled = as.numeric(scale(pupil_size)))

table(!is.na(joint_data$pupil_size_scaled))/nrow(joint_data)

joint_data$ITI.Duration <- as.factor(joint_data$ITI.Duration)
joint_data$last_ITI <- as.factor(joint_data$last_ITI)

write_csv(joint_data,paste0("/Users/rh/Desktop/VRMID-analysis/data/study2/per_second_data.csv"))