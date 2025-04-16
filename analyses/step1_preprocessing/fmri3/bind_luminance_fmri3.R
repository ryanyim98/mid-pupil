library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(magick)


folder_path <- "~/Desktop/midaffemo/analyses/pupil/luminance_calc_figs/" 
folder_subpath <- c("cue","feedback","others","rating")

luminance_results <- NULL

for (subfolder in folder_subpath){
  image_files <- list.files(paste0(folder_path,subfolder), pattern = "\\.(png)$", full.names = TRUE)
  
  #luminance <- array(NA, dim = c(length(image_files), 256, 160))
  for (i in seq_along(image_files)) {
    image <- image_read(image_files[i])
    image_name <- strsplit(image_files[i],"/")
    image_type = image_name[[1]][9]
    image_id = strsplit(image_name[[1]][10],"capturedImage_")[[1]][2]
    image_event = strsplit(image_id,".p")[[1]][1]
    
    image_rgb <- image_convert(image, colorspace = "RGB")
    
    image_data <- as.integer(image_data(image_rgb))
    
    R <- image_data[,,1]
    G <- image_data[,,2]
    B <- image_data[,,3]
    
    luminance <- 0.2126 * R + 0.7152 * G + 0.0722 * B
    
    average_luminance <- mean(luminance)
    
    luminance_results <- rbind(luminance_results,c(image_type, image_id,image_event, average_luminance))
  }
  
}

luminance_results <- as.data.frame(luminance_results)
names(luminance_results) <- c("type","img_id","img_class","luminance")
write_csv(luminance_results, "../../data/fmri3/derivatives/screen_luminance.csv")

fmri3_pupil_data <- fmri3_pupil_data %>% filter(!is.na(event)) %>%   # for now filter out NA values
  mutate(vrating_t = ifelse(vrating_t == -1,4,vrating_t),
         arating_t = ifelse(arating_t == -1,4,arating_t),
         erating_t = ifelse(erating_t == -1,8,erating_t))
subject_n <- unique(fmri3_pupil_data$subject)

luminance_values <- setNames(as.numeric(luminance_results$luminance), luminance_results$img_id)

tail <- ".png"

# Initialize a new column for luminance
fmri3_pupil_data$luminance <- NA

fmri3_pupil_data <- fmri3_pupil_data %>%
  group_by(subject,trial) %>%
  group_split() %>%
  lapply(function(trial_data) {
    # Extract trial-specific variables
    print(paste0("subject:", unique(trial_data$subject)))
    print(paste0("trial no:", unique(trial_data$trial)))
    
    target_seconds <- unique(trial_data$target_ms) # Convert target_ms to seconds
    reactiontime <- unique(trial_data$rt)                 # Reaction time in seconds
    if (reactiontime < 0) {    # missing reaction time
      reactiontime = 0
    }
    
    probe_type <- unique(trial_data$probe)  # probe type
    rating_type <- unique(trial_data$rating_type)  # probe type
    
    cue_size_str <- unique(ifelse(trial_data$cue_size == 1, "cue_small", "cue_large"))  # if cue is small or large
    
    # fetching the right luminance value depending on which amount of money is shown
    cue_money_type <- unique(case_when(
      trial_data$trialtype == 1 ~ "minus0",
      trial_data$trialtype == 2 ~ "minus1",
      trial_data$trialtype == 3 ~ "minus5",
      trial_data$trialtype == 4 ~ "plus0",
      trial_data$trialtype == 5 ~ "plus1",
      trial_data$trialtype == 6 ~ "plus5",
      TRUE ~ NA_character_ 
    ))
    
    cue_money_str = paste(cue_money_type, tail, sep = "")
    cue_type = paste(cue_size_str, cue_money_str, sep = "_")
    
    feedback_str <- unique(case_when(
      # Lose money trial (trialtype 1:3)
      trial_data$trialtype %in% 1:3 & trial_data$outcome_type == 1 ~ "feedback_0",  # Hit target
      trial_data$trialtype %in% 1:3 & trial_data$outcome_type != 1 ~ paste("feedback", cue_money_type, sep = "_"),  # Missed target
      
      # Gain money trial (trialtype 4:6)
      trial_data$trialtype %in% 4:6 & trial_data$outcome_type != 1 ~ "feedback_0",  # Missed target
      trial_data$trialtype %in% 4:6 & trial_data$outcome_type == 1 ~ paste("feedback", cue_money_type, sep = "_"),  # Hit target
    ))
    
    feedback_str <- paste(feedback_str, tail, sep = "")
    
    # # for how aroused, how do you feel, how positive do you feel
    # rating_str <- case_when(
    #   trial_data$rating_type == 1 ~ "arousal",
    #   trial_data$rating_type == 2 ~ "emotion",
    #   #trial_data$rating_type == 1 ~ "valence",
    # )
    # rating_str <- paste(rating_str, tail, sep = "")
    
    # Map luminance values based on sample_in_trial_t time intervals
    if (probe_type == 1) {
      if (rating_type == 1){ # valence and arousal ratings
        # affective probe after anticipatory fixation
        trial_data <- trial_data %>%
          mutate(
            luminance = case_when(
              sample_in_trial_t < 0 ~ luminance_values["iti.png"],  # all the seconds in baselining
              sample_in_trial_t >= 0 & sample_in_trial_t < 2 ~ luminance_values[cue_type],                         # Cue: 0–2 seconds
              sample_in_trial_t >= 2 & sample_in_trial_t < 4 ~ luminance_values["fix.png"],                   # Fixation: 2–4 seconds
              sample_in_trial_t >= 4 & sample_in_trial_t < 4+vrating_t ~ luminance_values['valence.png'],  # two ratings
              sample_in_trial_t >= 4+vrating_t & sample_in_trial_t < 8 ~ luminance_values['valence_yellow.png'],  # two ratings
              sample_in_trial_t >= 8 & sample_in_trial_t < 8+arating_t ~ luminance_values["arousal.png"],
              sample_in_trial_t >= 8+arating_t & sample_in_trial_t < 12 ~ luminance_values["arousal_yellow.png"],
              sample_in_trial_t >= 12 & sample_in_trial_t < 14+0.625 ~ luminance_values["fix.png"],
              sample_in_trial_t >= 14+0.625 & sample_in_trial_t < (14 + 0.625 + target_seconds) ~ luminance_values["target.png"],
              sample_in_trial_t >= (14 +0.625 + target_seconds) & sample_in_trial_t < 16 ~ luminance_values["allblack.png"],
              sample_in_trial_t >= 16 & sample_in_trial_t < 18 ~ luminance_values[feedback_str], 
              sample_in_trial_t >= 18 ~ luminance_values["iti.png"], # intertrial
              TRUE ~ NA_real_  # Default case (optional)
            )
          )
      } else if (rating_type == 2){ #emotion ratings
        # affective probe after anticipatory fixation
        trial_data <- trial_data %>%
          mutate(
            luminance = case_when(
              sample_in_trial_t < 0 ~ luminance_values["iti.png"],  # all the seconds in baselining
              sample_in_trial_t >= 0 & sample_in_trial_t < 2 ~ luminance_values[cue_type],                         # Cue: 0–2 seconds
              sample_in_trial_t >= 2 & sample_in_trial_t < 4 ~ luminance_values["fix.png"],                   # Fixation: 2–4 seconds
              sample_in_trial_t >= 4 & sample_in_trial_t < 4+erating_t ~ luminance_values['emotion.png'],  # two ratings
              sample_in_trial_t >= 4+erating_t & sample_in_trial_t < 12 ~ luminance_values['emotion_yellow.png'],  # two ratings
              sample_in_trial_t >= 12 & sample_in_trial_t < 14+0.625 ~ luminance_values["fix.png"],
              sample_in_trial_t >= 14+0.625 & sample_in_trial_t < (14 + 0.625 + target_seconds) ~ luminance_values["target.png"],
              sample_in_trial_t >= (14 +0.625 + target_seconds) & sample_in_trial_t < 16 ~ luminance_values["allblack.png"],
              sample_in_trial_t >= 16 & sample_in_trial_t < 18 ~ luminance_values[feedback_str], 
              sample_in_trial_t >= 18 ~ luminance_values["iti.png"], # intertrial
              TRUE ~ NA_real_  # Default case (optional)
            )
          )
      }
      
    } else if (probe_type == 2) {
      if (rating_type == 1){
        # affective probe after trial outcome feedback
        trial_data <- trial_data %>%
          mutate(
            luminance = case_when(
              sample_in_trial_t < 0 ~ luminance_values["iti.png"],  # all the seconds in baselining
              sample_in_trial_t >= 0 & sample_in_trial_t < 2 ~ luminance_values[cue_type],                         # Cue: 0–2 seconds
              sample_in_trial_t >= 2 & sample_in_trial_t < 4+0.625 ~ luminance_values["fix.png"],                   # Fixation: 2–4 seconds
              sample_in_trial_t >= 4+0.625 & sample_in_trial_t < (4+0.625 + target_seconds) ~ luminance_values["target.png"],  # Target: 4–(4 + target_seconds)
              sample_in_trial_t >= (4 +0.625+ target_seconds) & sample_in_trial_t < 6 ~ luminance_values["allblack.png"], # Blank screen
              sample_in_trial_t >= 6 & sample_in_trial_t < 8 ~ luminance_values[feedback_str],                   # Feedback: 6–8 seconds
              sample_in_trial_t >= 8 & sample_in_trial_t < 8+vrating_t ~ luminance_values['valence.png'],             # Rating screen: 8–16 seconds
              sample_in_trial_t >= 8+vrating_t & sample_in_trial_t < 12 ~ luminance_values['valence_yellow.png'],             # Rating screen: 8–16 seconds
              sample_in_trial_t >= 12 & sample_in_trial_t < 12+arating_t ~ luminance_values['arousal.png'],             # Rating screen: 8–16 seconds
              sample_in_trial_t >= 12+arating_t & sample_in_trial_t < 16 ~ luminance_values['arousal_yellow.png'],             # Rating screen: 8–16 seconds
              sample_in_trial_t >= 16  ~ luminance_values["iti.png"], # Intertrial: 16–(16 + iti_seconds)
              TRUE ~ NA_real_  # Default case (optional)
            )
          )
      } else if (rating_type == 2){
        # affective probe after trial outcome feedback
        trial_data <- trial_data %>%
          mutate(
            luminance = case_when(
              sample_in_trial_t < 0 ~ luminance_values["iti.png"],  # all the seconds in baselining
              sample_in_trial_t >= 0 & sample_in_trial_t < 2 ~ luminance_values[cue_type],                         # Cue: 0–2 seconds
              sample_in_trial_t >= 2 & sample_in_trial_t < 4+0.625 ~ luminance_values["fix.png"],                   # Fixation: 2–4 seconds
              sample_in_trial_t >= 4+0.625 & sample_in_trial_t < (4+0.625 + target_seconds) ~ luminance_values["target.png"],  # Target: 4–(4 + target_seconds)
              sample_in_trial_t >= (4 +0.625+ target_seconds) & sample_in_trial_t < 6 ~ luminance_values["allblack.png"], # Blank screen
              sample_in_trial_t >= 6 & sample_in_trial_t < 8 ~ luminance_values[feedback_str],                   # Feedback: 6–8 seconds
              sample_in_trial_t >= 8 & sample_in_trial_t < 8+erating_t ~ luminance_values['emotion.png'],             # Rating screen: 8–16 seconds
              sample_in_trial_t >= 8+erating_t & sample_in_trial_t < 16 ~ luminance_values['emotion_yellow.png'],             # Rating screen: 8–16 seconds
              sample_in_trial_t >= 16 ~ luminance_values["iti.png"], # Intertrial: 16–(16 + iti_seconds)
              TRUE ~ NA_real_  # Default case (optional)
            )
          )
      }
      
    }
    
    return(trial_data)
  }) %>%
  bind_rows()  # Combine back into a single data frame

# create lagged luminance
fmri3_pupil_data <- fmri3_pupil_data %>% 
  group_by(subject,trial) %>% 
  mutate(luminance_lag1 = ifelse(sample_in_trial_t <= 1,  luminance_values["iti.png"],
    lag(luminance,200))) #the luminance that matters is the luminance 1 sec ago or so
#write_csv(fmri3_pupil_data,"../../data/fmri3/derivatives/pupillometry_baselineCorrected_with_luminance.csv")
