library(readr)
library(dplyr)
library(factoextra)
library(tidyr)
library(ggplot2)
library(patchwork)

component_colors <- c(RColorBrewer::brewer.pal(10, "Paired"), "black")
names(component_colors) <- c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6","PC7","PC8","PC9","PC10","pupil")

my.data <- read_csv("~/Desktop/VRMID-analysis/vrmid-pupillometry/data/study1/pupil_for_pca.csv")
subject_n <- unique(my.data$Subject)

##### caculate all the pca together
## filter out outlier
my.data.mean <- my.data %>% 
  rename(subject = Subject) %>% 
  mutate(probe = ifelse(probe == "anti",1,
                        ifelse(probe == "out",2,NA)))


ant <- my.data.mean %>% 
  select(-Time,-current_stimulus) %>% 
  filter(sample_in_trial_t > 0,
         sample_in_trial_t <= 18) %>% 
  ungroup() %>% 
  filter(probe == 1)

counter = 1
for (i in 1:length(subject_n)){
  n_trials <- unique(ant[ant$subject == subject_n[i],'trial'])
  for (t in 1:length(n_trials$trial)){
    trialn = n_trials$trial[t]
    ant[ant$subject == subject_n[i] & ant$trial == trialn,'trial_id_unique'] <- counter
    counter = counter + 1
  }
}

ant <- ant%>% 
  ungroup() %>% 
  group_by(trial_id_unique) %>% 
  mutate(trial_nrow = max(row_number())) %>% 
  filter(trial_nrow > 120 * 17.8)

out <-my.data.mean %>% 
  select(-Time,-current_stimulus) %>% 
  filter(sample_in_trial_t > 0,
         sample_in_trial_t <= 18) %>% 
  filter(probe == 2)

counter = 1
for (i in 1:length(subject_n)){
  n_trials <- unique(out[out$subject == subject_n[i],'trial'])
  for (t in 1:length(n_trials$trial)){
    trialn = n_trials$trial[t]
    out[out$subject == subject_n[i]& out$trial == trialn,'trial_id_unique'] <- counter
    counter = counter + 1
  }
}


out <- out%>% 
  ungroup() %>% 
  group_by(trial_id_unique) %>% 
  mutate(trial_nrow = max(row_number())) %>% 
  filter(trial_nrow > 120 * 17.8)

##method 1
##caculate all the pca together, but for each subject's each trial
##method1

ant_mean <- ant %>% 
  group_by(sample_in_trial_t) %>% 
  summarise(pupilDiameter = mean(pupil_Avg_scale,na.rm = T)) %>% 
  ungroup()

ant_mean2 <- ant %>% 
  group_by(sample_in_trial_t,condition) %>% 
  summarise(pupilDiameter = mean(pupil_Avg_scale,na.rm = T)) %>% 
  ungroup()

out_mean <- out %>% 
  group_by(sample_in_trial_t) %>% 
  summarise(pupilDiameter = mean(pupil_Avg_scale,na.rm = T)) %>% 
  ungroup()

out_mean2 <- out %>% 
  group_by(sample_in_trial_t,condition) %>% 
  summarise(pupilDiameter = mean(pupil_Avg_scale,na.rm = T)) %>% 
  ungroup()

ggplot(ant_mean, aes(x = sample_in_trial_t, y = pupilDiameter))+
  geom_line()+
  ggplot(out_mean, aes(x = sample_in_trial_t, y = pupilDiameter))+
  geom_line()

ggplot(ant_mean2, aes(x = sample_in_trial_t, y = pupilDiameter, color = condition))+
  geom_line()+
  ggplot(out_mean2, aes(x = sample_in_trial_t, y = pupilDiameter, color = condition))+
  geom_line()

trial <- unique(ant$trial_id_unique);
trial_min <- array(NA, dim = c(length(trial)))
for (i in 1:length(trial)){
  var1<-ant[ant$trial_id_unique == trial[i], 'pupil_Avg_scale'];
  trial_min[i] <- length(var1$pupil_Avg_scale)
}

ant_pca_data  <- array(NA, dim = c(length(trial),min(trial_min))) 

trials_kept = c()
for (i in 1:length(trial)){
  var1<-ant[ant$trial_id_unique == trial[i], 'pupil_Avg_scale']
  ant_pca_data[i,]<-var1$pupil_Avg_scale[1:min(trial_min)]
  trials_kept = c(trials_kept,trial[i])
}

# Subset the matrix to only include rows without missing values
ant_pca_data_complete <- ant_pca_data[complete.cases(ant_pca_data), ]

ant_trial_included <- trials_kept[complete.cases(ant_pca_data)]

ant_pca <- prcomp(ant_pca_data_complete,center = TRUE,scale. = TRUE)

ncomp = dim(ant_pca$rotation)[2]

pca.loadings <- as_tibble(ant_pca$rotation[, 1:ncomp], rownames="time") %>% 
  mutate(time = as.numeric(time)*(1/120)) %>% 
  # mutate_at(vars(PC1:PC6), ~zoo::rollmean(., k = 100, fill = NA, align = "center"))%>% 
  left_join(ant_mean %>% rename(time= sample_in_trial_t,
                                pupil = pupilDiameter) %>% 
              mutate(pupil = pupil * 0.1), by = "time")
porigant<-
  pca.loadings %>%
  pivot_longer(cols = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6","PC7","PC8","PC9","PC10", "pupil"),
               names_to = "component",
               values_to = "loading") %>%
  ggplot(aes(x = time, y = loading, group = component)) +
  geom_line(aes(color = component, fill = component), linewidth = 1) +
  labs(x = "Time (s)", y = "Component loading", title = "Temporal PCA - anticipation probes") +
  scale_color_manual(values = component_colors) +
  scale_fill_manual(values = component_colors) +
  theme_bw() +
  geom_vline(xintercept = c(0, 4, 6, 14,16, 18), linetype = "dashed", color = "gray37") +
  annotate("text", x = 2, y = -0.05, label = "cue") +
  annotate("text", x = 5, y = -0.05, label = "fixation") +
  annotate("text", x = 10, y = -0.05, label = "rating") +
  annotate("text", x = 15, y = -0.05, label = "triang") +
  annotate("text", x = 17, y = -0.05, label = "out")

raw.loadings     <- ant_pca$rotation[,1:ncomp] %*% diag(ant_pca$sdev, ncomp, ncomp)
rotated.loadings <- varimax(raw.loadings)$loadings
# rotated.loadings2 <- psych::principal(ant_pca_data_complete,nfactors = 6,rotate = 'varimax')
inv.loadings     <- t(pracma::pinv(rotated.loadings))
loadings         <- rotated.loadings
scores_ant           <- scale(ant_pca_data_complete) %*% inv.loadings
subj.pred        <- as_tibble(scores_ant)

ss_loadings <- colSums(loadings^2)

var_explained <- as_tibble(ss_loadings/nrow(loadings)) %>% 
  mutate(percent_var_explained = value * 100) %>% 
  tibble::add_column(component = paste0("PC",1:ncomp))

var_explained <- var_explained %>% 
  arrange(desc(percent_var_explained)) %>% 
  mutate(component_renamed = paste("PC", row_number(), sep="")) %>% 
  rowwise() %>% 
  mutate(component_id = as.numeric(strsplit(component,"C")[[1]][2])) %>% 
  arrange(component_id)

pca.loadings.rot <- as_tibble(rotated.loadings[, 1:ncomp], rownames="time") %>% 
  mutate(time = as.numeric(time)*(1/120)) %>% 
  # mutate_at(vars(V1:V6), ~zoo::rollmean(., k = 100, fill = NA, align = "center"))%>% 
  left_join(ant_mean %>% rename(time= sample_in_trial_t,
                                pupil = pupilDiameter), by = "time")

names(pca.loadings.rot)[2:(ncomp+1)] <- var_explained$component_renamed #c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6","PC7","PC8","PC9","PC10")
#rename components based on variance explained
names(subj.pred) <- var_explained$component_renamed

sum(var_explained$percent_var_explained) #should be 100%

protant<-pca.loadings.rot %>%
  
  pivot_longer(cols = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6","PC7","PC8","PC9","PC10",
                        "PC11","PC12","PC13","PC14","PC15","pupil"),
               names_to = "component",
               values_to = "loading") %>%
  ggplot(aes(x = time, y = loading, group = component)) +
  geom_line(aes(color = component, fill = component), linewidth = 1) +
  labs(x = "Time (s)", y = "Component loading", title = "Temporal PCA - trial level - anticipation probes -varimax") +
  # scale_color_manual(values = component_colors) +
  # scale_fill_manual(values = component_colors) +
  theme_bw() +
  geom_vline(xintercept = c(0, 4, 6, 14,16, 18), linetype = "dashed", color = "gray37") +
  annotate("text", x = 2, y = -0.05, label = "cue") +
  annotate("text", x = 5, y = -0.05, label = "fixation") +
  annotate("text", x = 10, y = -0.05, label = "rating") +
  annotate("text", x = 15, y = -0.05, label = "triang") +
  annotate("text", x = 17, y = -0.05, label = "out")

pant <-porigant+protant

protant

#candidate_luminance_pc_ant <- c("PC1","PC2","PC5","PC6","PC8","PC10") #c("PC1","PC2","PC6","PC7","PC9","PC10") c("PC9")
#candidate_luminance_pc_ant <- c("PC5","PC10","PC12","PC15")
candidate_luminance_pc_ant <- c("PC1","PC4","PC9","PC10")

pca.loadings.rot %>%
  pivot_longer(cols = candidate_luminance_pc_ant,
               names_to = "component",
               values_to = "loading") %>%
  ggplot(aes(x = time, y = loading, group = component)) +
  geom_line(aes(color = component, fill = component), linewidth = 1) +
  labs(x = "Time (s)", y = "Component loading", title = "Temporal PCA - trial level - anticipation probes -varimax") +
  # scale_color_manual(values = component_colors) +
  # scale_fill_manual(values = component_colors) +
  theme_bw()  +
  geom_vline(xintercept = c(0, 4, 6, 14,16, 18), linetype = "dashed", color = "gray37") +
  annotate("text", x = 2, y = -0.05, label = "cue") +
  annotate("text", x = 5, y = -0.05, label = "fixation") +
  annotate("text", x = 10, y = -0.05, label = "rating") +
  annotate("text", x = 15, y = -0.05, label = "triang") +
  annotate("text", x = 17, y = -0.05, label = "out")


###### get 'cleaned' data
# Convert the candidate PCs to column indices
candidate_indices <- match(candidate_luminance_pc_ant, names(subj.pred))

# Zero out the scores of the candidate components
scores_ant[, candidate_indices] <- 0

sum(var_explained$percent_var_explained[!var_explained$component_renamed %in% candidate_luminance_pc_ant])

# Reconstruct the cleaned pupil data
cleaned_data <- scores_ant %*% t(rotated.loadings)

# Convert cleaned data to a tibble or data frame if needed
cleaned_data <- scale(cleaned_data, center = -ant_pca$center, scale = 1/ant_pca$scale)

cleaned_data <- as_tibble(cleaned_data)

#check the values
cor.test(as.numeric(cleaned_data[1,]),as.numeric(ant_pca_data[1,]))
mean(as.numeric(cleaned_data[1,])/ant_pca_data[1,])

#add other info
cleaned_data$trial_id_unique <- ant_trial_included
ant_info <- unique(ant %>% 
                     select(subject,trial:ITI_duration,arousal_scaled,valence_scaled,trial_id_unique))
cleaned_data <- left_join(cleaned_data,ant_info,by="trial_id_unique")

cleaned_data_long <- cleaned_data %>% 
  pivot_longer(V1:V2160, names_to = "sample_in_trial_n",names_prefix = "V", values_to = "cleaned_pupilDiameter") %>% 
  mutate(sample_in_trial_t = as.numeric(sample_in_trial_n) * (1/120))

cleaned_data_long_ant <- cleaned_data_long

cleaned_data_long_sum <- cleaned_data_long %>% 
  group_by(subject, sample_in_trial_t,condition) %>% 
  summarise(cleaned_pupilDiameter = mean(cleaned_pupilDiameter,na.rm = T)) %>% 
  group_by(sample_in_trial_t,condition) %>% 
  summarise(cleaned_pupilDiameter = mean(cleaned_pupilDiameter,na.rm = T)) 

pcleanant<-ggplot()+
  geom_line(data=cleaned_data_long_sum, aes(x = sample_in_trial_t, y = cleaned_pupilDiameter, color = condition),
            linewidth = 1)+
  ylim(-2.5,1.5)+
  labs(title = 'cleaned pupil',
       color = "cue luminance") +
  geom_vline(xintercept = c(0, 4, 6, 14,16, 18), linetype = "dashed", color = "gray37") +
  annotate("text", x = 2, y = -0.05, label = "cue") +
  annotate("text", x = 5, y = -0.05, label = "fixation") +
  annotate("text", x = 10, y = -0.05, label = "rating") +
  annotate("text", x = 15, y = -0.05, label = "triang") +
  annotate("text", x = 17, y = -0.05, label = "out")+
ggplot()+
  geom_line(data=ant_mean2, aes(x = sample_in_trial_t, y = pupilDiameter, color = as.factor(condition)),
            linewidth = 1)+
  ylim(-2.5,1.5)+
  labs(title = 'original pupil',
       color = "cue luminance") +
  geom_vline(xintercept = c(0, 4, 6, 14,16, 18), linetype = "dashed", color = "gray37") +
  annotate("text", x = 2, y = -0.05, label = "cue") +
  annotate("text", x = 5, y = -0.05, label = "fixation") +
  annotate("text", x = 10, y = -0.05, label = "rating") +
  annotate("text", x = 15, y = -0.05, label = "triang") +
  annotate("text", x = 17, y = -0.05, label = "out")
  #geom_line(data=ant_mean, aes(x = sample_in_trial_t, y = pupilDiameter),color = "blue")



### outcome

trial <- unique(out$trial_id_unique);
trial_min <- array(NA, dim = c(length(trial)))
for (i in 1:length(trial)){
  var1<-out[out$trial_id_unique == trial[i], 'pupil_Avg_scale'];
  trial_min[i] <- length(var1$pupil_Avg_scale)
}

out_pca_data  <- array(NA, dim = c(length(trial),min(trial_min))) 

trials_kept = c()
for (i in 1:length(trial)){
  var1<-out[out$trial_id_unique == trial[i], 'pupil_Avg_scale'];
  out_pca_data[i,]<-var1$pupil_Avg_scale[1:min(trial_min)];
  trials_kept = c(trials_kept,trial[i])
}

# Subset the matrix to only include rows without missing values
out_pca_data_complete <- out_pca_data[complete.cases(out_pca_data), ]
out_trial_included <- trials_kept[complete.cases(out_pca_data)]

out_pca <- prcomp(out_pca_data_complete,center = TRUE,scale. = TRUE)
summary(out_pca)

ncomp = dim(out_pca$rotation)[2]

pca.loadings <- as_tibble(out_pca$rotation[, 1:ncomp], rownames="time") %>%
  mutate(time = as.numeric(time)/120)%>% 
  # mutate_at(vars(PC1:PC6), ~zoo::rollmean(., k = 100, fill = NA, align = "left"))%>% 
  left_join(out_mean %>% rename(time= sample_in_trial_t,
                                pupil = pupilDiameter) %>% 
              mutate(pupil = pupil * 0.1), by = "time")

porigout<-pca.loadings %>%
  pivot_longer(cols = c("PC1", "PC2", "PC3", "PC4","PC5", "PC6","PC7","PC8","PC9","PC10","pupil"),  #, "PC5", "PC6"
               names_to = "component",
               values_to = "loading") %>%
  ggplot(aes(x = time,
             y = loading,
             color = component,
             fill = component)) +
  #stat_smooth(method="gam", se=F, alpha=0.3) +
  geom_line(linewidth = 1)+

  labs(x = "Time (s)",
       y = "Component loading",
       title =  "Temporal PCA - trial level - outcome probes") +
  scale_color_manual(values = component_colors) +
  scale_fill_manual(values = component_colors) +
  theme_bw()  +
  geom_vline(xintercept = c(0,4,6,8,10,18),
             linetype = "dashed",
             color = "gray37") +
  annotate("text",x = 2, y = -0.05,label = "cue")+
  annotate("text",x = 5, y = -0.05,label = "fixation")+
  annotate("text",x = 7, y = -0.05,label = "triang")+
  annotate("text",x = 9, y = -0.05,label = "out")+
  annotate("text",x = 14, y = -0.05,label = "rating")

raw.loadings     <- out_pca$rotation[,1:ncomp] %*% diag(out_pca$sdev, ncomp, ncomp)
rotated.loadings <- varimax(raw.loadings)$loadings
inv.loadings     <- t(pracma::pinv(rotated.loadings))
loadings         <- rotated.loadings
scores_out           <- scale(out_pca_data_complete) %*% inv.loadings
subj.pred        <- as_tibble(scores_out)
#names(subj.pred) <- paste0("PC",1:ncomp)

ss_loadings <- colSums(loadings^2)

var_explained <- as_tibble(ss_loadings/nrow(loadings)) %>% 
  mutate(percent_var_explained = value * 100) %>% 
  tibble::add_column(component = paste0("PC",1:ncomp))

var_explained <- var_explained %>% 
  arrange(desc(percent_var_explained)) %>% 
  mutate(component_renamed = paste("PC", row_number(), sep="")) %>% 
  rowwise() %>% 
  mutate(component_id = as.numeric(strsplit(component,"C")[[1]][2])) %>% 
  arrange(component_id)

pca.loadings.rot <- as_tibble(rotated.loadings[, 1:ncomp], rownames="time") %>% 
  mutate(time = as.numeric(time)*(1/120)) %>% 
  # mutate_at(vars(V1:V6), ~zoo::rollmean(., k = 100, fill = NA, align = "left"))%>% 
  left_join(out_mean %>% rename(time= sample_in_trial_t,
                                pupil = pupilDiameter), by = "time")

names(pca.loadings.rot)[2:(ncomp+1)] <- var_explained$component_renamed #c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6","PC7","PC8","PC9","PC10")
#rename components based on variance explained
names(subj.pred) <- var_explained$component_renamed

protout<-pca.loadings.rot %>%
  pivot_longer(cols = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6","PC7","PC8","PC9","PC10",
                        "PC11","PC12","PC13","PC14","PC15", "pupil"),
               names_to = "component",
               values_to = "loading") %>%
  ggplot(aes(x = time, y = loading, group = component)) +
  geom_line(aes(color = component, fill = component), linewidth = 1) +
  labs(x = "Time (s)", y = "Component loading", title = "- Temporal PCA - trial level - anticipation probes - varimax") +
  # scale_color_manual(values = component_colors) +
  # scale_fill_manual(values = component_colors) +
  theme_bw()   +
  labs(x = "Time (s)",
       y = "Component loading",
       title =  "Temporal PCA - trial level - outcome probes") +
  theme_bw()  +
  geom_vline(xintercept = c(0,4,6,8,10,18),
             linetype = "dashed",
             color = "gray37") +
  annotate("text",x = 2, y = -0.05,label = "cue")+
  annotate("text",x = 5, y = -0.05,label = "fixation")+
  annotate("text",x = 7, y = -0.05,label = "triang")+
  annotate("text",x = 9, y = -0.05,label = "out")+
  annotate("text",x = 14, y = -0.05,label = "rating")


pout <- porigout+protout
protout
###### get 'cleaned' data
candidate_luminance_pc_out<-c("PC1","PC2","PC9","PC10","PC11")
#candidate_luminance_pc_out <- c("PC2","PC3","PC4","PC5","PC6","PC8","PC9") #c("PC1","PC2","PC6","PC7","PC9","PC10") c("PC9")

pca.loadings.rot %>%
  pivot_longer(cols = candidate_luminance_pc_out,
               names_to = "component",
               values_to = "loading") %>%
  ggplot(aes(x = time, y = loading, group = component)) +
  geom_line(aes(color = component, fill = component), linewidth = 1) +
  labs(x = "Time (s)", y = "Component loading", title = "- Temporal PCA - trial level - anticipation probes - varimax") +
  theme_bw()   +
  labs(x = "Time (s)",
       y = "Component loading",
       title =  "Temporal PCA - trial level - outcome probes") +
  theme_bw() +
  geom_vline(xintercept = c(0,4,6,8,10,18),
             linetype = "dashed",
             color = "gray37") +
  annotate("text",x = 2, y = -0.05,label = "cue")+
  annotate("text",x = 5, y = -0.05,label = "fixation")+
  annotate("text",x = 7, y = -0.05,label = "triang")+
  annotate("text",x = 9, y = -0.05,label = "out")+
  annotate("text",x = 14, y = -0.05,label = "rating")


candidate_indices <- match(candidate_luminance_pc_out, names(subj.pred))

# Zero out the scores of the candidate components
scores_out[, candidate_indices] <- 0

sum(var_explained$percent_var_explained[!var_explained$component_renamed %in% candidate_luminance_pc_out])

# Reconstruct the cleaned pupil data
cleaned_data <- scores_out %*% t(rotated.loadings)

# Convert cleaned data to a tibble or data frame if needed
cleaned_data <- scale(cleaned_data, center = -out_pca$center, scale = 1/out_pca$scale)

cleaned_data <- as_tibble(cleaned_data)

cor.test(as.numeric(cleaned_data[1,]),as.numeric(out_pca_data[1,]))

#add other info
cleaned_data$trial_id_unique <- out_trial_included
out_info <- unique(out %>% 
                     select(subject,trial:ITI_duration,arousal_scaled,valence_scaled,trial_id_unique))
cleaned_data <- left_join(cleaned_data,out_info)

cleaned_data_long <- cleaned_data %>% 
  pivot_longer(V1:V2160, names_to = "sample_in_trial_n",names_prefix = "V", values_to = "cleaned_pupilDiameter") %>% 
  mutate(sample_in_trial_t = as.numeric(sample_in_trial_n) * (1/120))

cleaned_data_long_out <- cleaned_data_long

cleaned_data_long_sum <- cleaned_data_long %>% 
  group_by(subject, sample_in_trial_t,condition) %>% 
  summarise(cleaned_pupilDiameter = mean(cleaned_pupilDiameter,na.rm = T)) %>% 
  group_by(sample_in_trial_t,condition) %>% 
  summarise(cleaned_pupilDiameter = mean(cleaned_pupilDiameter,na.rm = T)) 

pcleanout<-ggplot()+
  geom_line(data=cleaned_data_long_sum, aes(x = sample_in_trial_t, y = cleaned_pupilDiameter, color = as.factor(condition)),
            linewidth = 1)+
  ylim(-2.5,1.5)+
  labs(title = 'cleaned pupil',
       color = "cue luminance")+
  geom_vline(xintercept = c(0,4,6,8,10,18),
             linetype = "dashed",
             color = "gray37") +
  annotate("text",x = 2, y = -0.05,label = "cue")+
  annotate("text",x = 5, y = -0.05,label = "fixation")+
  annotate("text",x = 7, y = -0.05,label = "triang")+
  annotate("text",x = 9, y = -0.05,label = "out")+
  annotate("text",x = 14, y = -0.05,label = "rating")+
  ggplot()+
  geom_line(data=out_mean2, aes(x = sample_in_trial_t, y = pupilDiameter, color = as.factor(condition)),
            linewidth = 1)+
  ylim(-2.5,1.5)+
  labs(title = 'original pupil',
       color = "cue luminance")+
  geom_vline(xintercept = c(0,4,6,8,10,18),
             linetype = "dashed",
             color = "gray37") +
  annotate("text",x = 2, y = -0.05,label = "cue")+
  annotate("text",x = 5, y = -0.05,label = "fixation")+
  annotate("text",x = 7, y = -0.05,label = "triang")+
  annotate("text",x = 9, y = -0.05,label = "out")+
  annotate("text",x = 14, y = -0.05,label = "rating")
#geom_line(data=ant_mean, aes(x = sample_in_trial_t, y = pupilDiameter),color = "blue")

pcleanant/pcleanout
ggsave("~/Desktop/VRMID-analysis/vrmid-pupillometry/figures/pca/study1_all_pupil_pca_cleaned.png",dpi = 300, height = 7, width = 14)

pant/pout
ggsave("~/Desktop/VRMID-analysis/vrmid-pupillometry/figures/pca/study1_all_pupil_pca.png",dpi = 300, height = 7, width = 14)


protant/protout
ggsave("~/Desktop/VRMID-analysis/vrmid-pupillometry/figures/pca/study1_all_pupil_pca_rot.png",dpi = 300, height = 10, width = 6)


all_clean_data <- rbind(cleaned_data_long_ant,cleaned_data_long_out)
write_csv(all_clean_data,"~/Desktop/VRMID-analysis/vrmid-pupillometry/data/study1/pupil_cleaned.csv")
