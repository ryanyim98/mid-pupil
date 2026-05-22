# Resolve paths under data/ whether getwd() is repo root, step2_pooled_analyses, or this helper/.
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

## vrmid1
vrmid1_pupil_data <- read_csv("~/Desktop/VRMID-analysis/mid-pupil/data/vrmid1/derivatives/pupillometry_lowpass_baselineCorrected.csv")%>%
  group_by(Subject)%>%
  mutate(arousal_scaled = as.numeric(scale(arousal)),
         valence_scaled = as.numeric(scale(valence)))%>%
  ungroup()

vrmid1_pupil_data <- vrmid1_pupil_data %>%
  rename_with(~ str_replace(., "_scale$", "_scaled"))

vrmid1_bad_participants <- read_delim(
  mid_pupil_resolve_data("vrmid1", "subjects_list", "subject_to_drop_pupil.txt"),
  col_names = F,
  delim = " "
)

vrmid1_pupil_data <- vrmid1_pupil_data %>% 
  filter(!Subject %in% c("cn221206","vx220916"), #manual bad participants
         !Subject %in% vrmid1_bad_participants) 

vrmid1_subjects <- unique(vrmid1_pupil_data$Subject)
vrmid1_subject_n<-length(vrmid1_subjects)
vrmid1_subject_n

ggplot(vrmid1_pupil_data %>% filter(Subject == "ay220916",trial < 5))+
  geom_line(aes(x = sample_in_trial_t, y = pupil_Avg_xyvar_scaled),color='red')+
  geom_line(aes(x = sample_in_trial_t, y = pupil_Avg_scaled))+
  facet_wrap(~trial)+
  labs(title = 1)

ggplot(vrmid1_pupil_data %>% filter(Subject == "ay220916",trial < 5))+
  geom_line(aes(x = sample_in_trial_t, y = AvgPPos_x_preproc),color='red')+
  geom_line(aes(x = sample_in_trial_t, y = AvgPPos_y_preproc))+
  facet_wrap(~trial)+
  labs(title = 1)


for (i in 1:vrmid1_subject_n){
  ggplot(vrmid1_pupil_data %>% filter(Subject == vrmid1_subjects[i]), aes(x = sample_in_trial_t, y = pupil_Avg_scaled))+
    geom_line()+
    facet_wrap(~trial)+
    labs(title = i)
  
  ggsave(paste0("~/Desktop/VRMID-analysis/mid-pupil/figures/vrmid1/examine_data/trial_figs/",vrmid1_subjects[i],".png"),dpi = 300, height = 7, width = 14)
}

vrmid1_pupil_data$trial_type <- factor(vrmid1_pupil_data$trial_type, levels = c("-$5","-$1","-$0","+$0","+$1","+$5"))

vrmid1_pupil_data <- vrmid1_pupil_data%>%
  mutate(sample_in_trial_t = round(sample_in_trial_t,4))

vrmid1_pupil_data<-vrmid1_pupil_data %>% 
  filter(abs(pupil_Avg_scaled)<5)


vrmid1_beh_data <- vrmid1_pupil_data %>% 
  ungroup() %>% 
  group_by(Subject,trial,probe,trial_type) %>% 
  slice(1) %>% 
  select(Subject,trial,probe,trial_type,reaction_outcome,reaction_time)


## fmri3
### ONLY if first time, run the code below to regress out the luminance component
# 
# fmri3_pupil_data <- read_csv("~/Desktop/VRMID-analysis/mid-pupil/data/fmri3/derivatives/pupillometry_baselineCorrected.csv") %>% 
#    mutate(cue_value = factor(trialtype, levels = 1:6, labels = c("-$0","-$1","-$5","+$0","+$1","+$5"))) %>% 
#    filter(!is.na(subject)) %>% 
#    ungroup() 
#   
#  fmri3_bad_participants <- read_delim("../../data/fmri3/subjects_list/subject_to_drop_pupil.txt", col_names = F, delim = " ")
#   
#  fmri3_pupil_data <- fmri3_pupil_data %>% 
#    filter(!sub_id %in% fmri3_bad_participants) %>% 
#    filter(!subject %in% c("sk240301","hc240324","js240311","jy240308","tc240311","el240312"))#these are manual exclusions based on data collection notes
# 
# # # ## run the code for binding luminance data to fmri3
# source("../step1_preprocessing/fmri3/bind_luminance_fmri3.R", local = TRUE)
# source("../step1_preprocessing/fmri3/convolve_luminance_fmri3.R", local = TRUE)

  
fmri3_pupil_data <- read_csv("~/Desktop/VRMID-analysis/mid-pupil/data/fmri3/derivatives/pupillometry_baselineCorrected_with_luminance_convolved.csv") %>% 
  mutate(cue_value = factor(trialtype, levels = 1:6, labels = c("-$0","-$1","-$5","+$0","+$1","+$5"))) %>% 
  filter(!is.na(subject)) %>% 
  ungroup() 

fmri3_subjects <- unique(fmri3_pupil_data$subject)
fmri3_subject_n<-length(fmri3_subjects)
fmri3_subject_n

fmri3_pupil_data$cue_value <- factor(fmri3_pupil_data$cue_value,
                                     levels = c("-$5","-$1","-$0","+$0","+$1","+$5"))

fmri3_pupil_data$emotion <- fmri3_pupil_data$erating

fmri3_pupil_data$emotion <- factor(fmri3_pupil_data$emotion, levels = 1:6, 
                                   labels = c("Anxious","Angry","Sad","Calm","Happy","Excited"))

fmri3_pupil_data$arousal <- fmri3_pupil_data$arating
fmri3_pupil_data$valence <- fmri3_pupil_data$vrating

fmri3_pupil_data$hit <- factor(fmri3_pupil_data$hit,levels = c(1,0), labels = c("Hit","Miss"))

fmri3_pupil_data <- fmri3_pupil_data %>% 
  mutate(emotion_cat = ifelse(emotion %in% c("Anxious","Angry","Excited"),"higharous",
                              ifelse(emotion %in% c("Sad","Calm","Happy"),"lowarous",NA)),
         emotion_cat2 = ifelse(emotion %in% c("Angry","Anxious","Sad"),
                               "pos_valence",
                               ifelse(emotion %in% c("Calm","Happy","Excited"),
                                      "neg_valence",NA))) %>% 
  ungroup()%>%
  group_by(subject) %>% 
  mutate(arousal_scaled = as.numeric(scale(arousal)),
         valence_scaled = as.numeric(scale(valence)))

fmri3_pupil_data<-fmri3_pupil_data %>% 
  filter(abs(pupilDiameter_scaled)<5) %>% 
  group_by(subject, trial, Time_sec) %>% 
  mutate(sample_in_sec = row_number())

for (i in 1:fmri3_subject_n){
  ggplot(fmri3_pupil_data %>% filter(subject == fmri3_subjects[i]), aes(x = sample_in_trial_t, y = pupilDiameter_scaled))+
    geom_line()+
    facet_wrap(~trial)+
    labs(title = i)
  
  ggsave(paste0("~/Desktop/VRMID-analysis/mid-pupil/figures/fmri3/examine_data/trial_figs/",fmri3_subjects[i],".png"),dpi = 300, height = 7, width = 14)
}


fmri3_pupil_data$cue_size <-factor(fmri3_pupil_data$cue_size, levels = c(1,2), labels = c("large","small"))

fmri3_pupil_data <- fmri3_pupil_data %>%
  filter(subject != "ag240330") #due to low hit rate
fmri3_subjects <- fmri3_subjects[fmri3_subjects!= "ag240330"]

fmri3_beh_data <- fmri3_pupil_data %>% 
  ungroup() %>% 
  group_by(subject,trial,probe,cue_value) %>% 
  slice(1) %>% 
  select(subject,trial,probe,rating_type,cue_value,binned_winpercent,hit,rt,
         vrating_t,arating_t,vrating,arating,scaledir,erating,erating_t)

fmri3_pupil_data$cue_size <- forcats::fct_recode(fmri3_pupil_data$cue_size,
                                                 small = "large",
                                                 large = "small")


## create time series of luminance-controlled pupil size
# 
# summary(l3<-lmer(scale(pupilDiameter_scaled) ~ scale(convolved_s1) + (1|subject),fmri3_pupil_data))
# summary(l4<-lmer(scale(pupilDiameter_scaled) ~ scale(convolved_s2) + (1|subject),fmri3_pupil_data))
summary(l5<-lmer(scale(pupilDiameter_scaled) ~ scale(convolved_s1)+scale(convolved_s2) + (1|subject),fmri3_pupil_data))

# summary(l3)$coefficient
# summary(l4)$coefficient
summary(l5)$coefficient
# anova(l3,l4,l5)
r.squaredGLMM(l5)

# # Extract rows used in the model
# model_data3 <- model.frame(l3)
# used_rows3 <- as.numeric(rownames(model_data3))
# # Create a residuals column in processed_data, initializing with NA
# fmri3_pupil_data$residuals_conv_system1 <- NA
# 
# residuals_l3 <- residuals(l3)
# # Assign residuals to the appropriate rows
# fmri3_pupil_data$residuals_conv_system1[used_rows3] <- residuals_l3

model_data5 <- model.frame(l5)
used_rows5 <- as.numeric(rownames(model_data5))
# Create a residuals column in processed_data, initializing with NA
fmri3_pupil_data$residuals_conv_system1 <- NA

residuals_l5 <- residuals(l5)
# Assign residuals to the appropriate rows
fmri3_pupil_data$residuals_conv[used_rows5] <- residuals_l5


## calculate normalized gaze position
fmri3_pupil_data <- fmri3_pupil_data %>%
  mutate(
  pupil_x_preproc_norm = pupil_x_preproc / (2560 - 1),  # -1 because 0-indexed (0 to 2559)
  pupil_y_preproc_norm = pupil_y_preproc / (1600 - 1))
