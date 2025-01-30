#convolve luminance with a gamma function
k_orig=0.276
k_new=2
theta=0.09
fmri3_pupil_data.new = NULL
# gamma_pdf <- gamma_pdf[!is.na(gamma_pdf)]
# Perform the convolution grouped by subject and trial
for (sub in unique(fmri3_pupil_data$subject)){
  print(sub)
  my.df.temp <- fmri3_pupil_data %>% filter(subject == sub)
  for (t in unique(my.df.temp$trial)){
    my.df.temp.trial <- my.df.temp %>% 
      filter(trial == t)
    tval <- my.df.temp.trial$sample_in_trial_t
    gamma_pdf_new <- dgamma(seq(0,20,0.005), shape = k_new, scale = theta) #dgamma(t_values, shape = k, scale = theta)
    gamma_pdf_orig<-dgamma(seq(0,20,0.005), shape = k_orig, scale = theta)
    gamma_pdf_orig[is.infinite(gamma_pdf_orig)] <- 0
    gamma_pdf_orig <- gamma_pdf_orig / sum(gamma_pdf_orig)  # Normalize to sum to 1
    gamma_pdf_new <- gamma_pdf_new / sum(gamma_pdf_new)  # Normalize to sum to 1
    my.df.temp.trial$convolved_orig <- convolve(my.df.temp.trial$luminance, rev(gamma_pdf_orig), type = "open")[1:(nrow(my.df.temp.trial))]
    my.df.temp.trial$convolved_new <- convolve(my.df.temp.trial$luminance,rev(gamma_pdf_new), type = "open")[1:(nrow(my.df.temp.trial))]
    
    # my.df.temp.trial$convolved_orig <- convolve(my.df.temp.trial$luminance, gamma_pdf_orig,type = "open")[1:(nrow(my.df.temp.trial))]
    # my.df.temp.trial$convolved_new <- convolve(my.df.temp.trial$luminance,gamma_pdf_new, type = "open")[1:(nrow(my.df.temp.trial))]
    fmri3_pupil_data.new <- rbind(fmri3_pupil_data.new,my.df.temp.trial)
  }
}
fmri3_pupil_data=fmri3_pupil_data.new
rm(fmri3_pupil_data.new)
write_csv(fmri3_pupil_data,'~/Desktop/VRMID-analysis/mid-pupil/data/fmri3/derivatives/pupillometry_baselineCorrected_with_luminance_convolved.csv')
