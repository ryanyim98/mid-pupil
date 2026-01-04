# The AFNI way

#pupil (alpha = .10, p < 0.001) NN2.2sided
3dClusterize -nosum -1Dformat -inset ../glm_maps/model33b_pupil_scaled_convolved_24z.nii.gz \
			-idat 1 -ithr 1 -NN 2 -clust_nvox 30 -2sided -3.227 3.227 -abs_table_data > tables/table_model33b_Z3.227.1D

			3dClusterize -nosum -1Dformat -inset ../glm_maps/model33c_pupil_scaled_convolved_24z.nii.gz \
						-idat 1 -ithr 1 -NN 2 -clust_nvox 30 -2sided -3.227 3.227 -abs_table_data > tables/table_model33c_Z3.227.1D
						
	#arousal (alpha = .10, p < 0.001) NN2.2sided
3dClusterize -nosum -1Dformat -inset ../glm_maps/model27c_param_arous_early_30z.nii.gz \
			-idat 1 -ithr 1 -NN 2 -clust_nvox 30 -2sided -3.227 3.227 -abs_table_data > tables/table_model27_Z3.227.1D
