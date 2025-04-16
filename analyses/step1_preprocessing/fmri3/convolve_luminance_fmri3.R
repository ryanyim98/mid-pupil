#convolve luminance with a gamma function
k_s1=2.40
k_s2=3.24
c_s1=0.77
c_s2=0.43
theta_s1=0.29
theta_s2=0.18
fmri3_pupil_data.new = NULL
# gamma_pdf <- gamma_pdf[!is.na(gamma_pdf)]
# Perform the convolution grouped by subject and trial
for (sub in unique(fmri3_pupil_data$subject)){
  print(sub)
  my.df.temp <- fmri3_pupil_data %>% filter(subject == sub)
  for (t in unique(my.df.temp$trial)){
    my.df.temp.trial <- my.df.temp %>% 
      filter(trial == t)
    my.df.temp.trial$luminance_lag <- NULL
    my.df.temp.trial$luminance_lag <- c(rep(0,40),my.df.temp.trial$luminance)[1:nrow(my.df.temp.trial)] #lag by 200 ms
    tval <- my.df.temp.trial$sample_in_trial_t
    gamma_pdf_s1<-c_s1*dgamma(seq(0,20,0.005), shape = k_s1, scale = theta_s1)
    gamma_pdf_s1[is.infinite(gamma_pdf_s1)] <- 0
    gamma_pdf_s1 <- gamma_pdf_s1 / sum(gamma_pdf_s1)  # Normalize to sum to 1
    my.df.temp.trial$convolved_s1 <- convolve(my.df.temp.trial$luminance_lag, rev(gamma_pdf_s1), type = "open")[1:(nrow(my.df.temp.trial))]
    
  
    gamma_pdf_s2<-c_s2*dgamma(seq(0,20,0.005), shape = k_s2, scale = theta_s2)
    gamma_pdf_s2[is.infinite(gamma_pdf_s2)] <- 0
    gamma_pdf_s2 <- gamma_pdf_s2 / sum(gamma_pdf_s2)  # Normalize to sum to 1
    my.df.temp.trial <- my.df.temp.trial %>% 
      mutate(luminance_change = ifelse(trial == 1 & sample_in_trial_n == 2, luminance -  luminance_values["iti.png"],
                                       luminance - lag(luminance,1)),
             luminance_change = ifelse(is.na(luminance_change),0,
                                       ifelse(luminance_change > 0, 1,0)))
    my.df.temp.trial$luminance_change <- c(rep(0,40),my.df.temp.trial$luminance_change)[1:nrow(my.df.temp.trial)] #lag by 200 ms
    my.df.temp.trial$convolved_s2 <- convolve(my.df.temp.trial$luminance_change, rev(gamma_pdf_s2), type = "open")[1:(nrow(my.df.temp.trial))]
    # my.df.temp.trial$convolved_total = my.df.temp.trial$convolved_s1+200*my.df.temp.trial$convolved_s2
    # ggplot(my.df.temp.trial)+
    #   geom_line(aes(x = sample_in_trial_t, y = luminance), color = "red")+
    #   geom_line(aes(x = sample_in_trial_t, y = luminance_change))+
    #   geom_line(aes(x = sample_in_trial_t, y = convolved_s1), color = "red",size = 2)+
    #   geom_line(aes(x = sample_in_trial_t, y = convolved_s2*200),size = 2)
    fmri3_pupil_data.new <- rbind(fmri3_pupil_data.new,my.df.temp.trial)
  }
}
fmri3_pupil_data=fmri3_pupil_data.new
rm(fmri3_pupil_data.new)

write_csv(fmri3_pupil_data,'~/Desktop/VRMID-analysis/mid-pupil/data/fmri3/derivatives/pupillometry_baselineCorrected_with_luminance_convolved.csv')
