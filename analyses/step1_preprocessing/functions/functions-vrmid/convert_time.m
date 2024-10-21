function [times_conv, pupil_position_conv, pupil_size_conv] =  convert_time(times,fmt,sr,pupil_position,pupil_size) 

secs = unique(times);
times_conv.timess = repelem(secs,sr); %120 Hz strictly
times_conv.time = repelem(secs,sr); %120 Hz strictly
pupil_position_conv = array2table(NaN([length(times_conv.timess) size(pupil_position,2)]),'VariableNames',pupil_position.Properties.VariableNames);
pupil_size_conv = array2table(NaN([length(times_conv.timess) size(pupil_size,2)]),'VariableNames',pupil_size.Properties.VariableNames);

times_conv.timess.Format = fmt;
times_conv.time.Format = "hh:mm:ss";
first_sec_ind = find(times == min(secs));
first_sec_time = times(times == min(secs),:);
this_second = secs(1);
sampling_duration = seconds(1/sr);

%if first second has less than sr data points, count backwards to
%recover time
if length(first_sec_time) <= sr
    for i = 1:length(first_sec_time)
            first_sec_time(i) = first_sec_time(i) + seconds(1) - sampling_duration *...
                (length(first_sec_time)-i+1);
    end
elseif length(first_sec_time) > sr
    for i = 1:sr
            first_sec_time(i) = first_sec_time(i) + seconds(1) - sampling_duration * ...
                (length(first_sec_time) - i+1);
    end
end

pupil_position_conv((121-length(first_sec_ind)):sr,:) =  pupil_position(first_sec_ind,:);
pupil_size_conv((121-length(first_sec_ind)):sr,:) =  pupil_size(first_sec_ind,:);

first_sec_time = [NaN(sr-length(first_sec_time),1); first_sec_time];
first_sec_time.Format = fmt;

times_conv.timess(1:length(first_sec_time)) = first_sec_time;

%move to the next sec
for j = 2:length(secs)
    this_sec_time = times(times == secs(j),:); 
    this_sec_ind_old = find(times == secs(j));
    this_sec_ind_new = find(times_conv.timess == secs(j));
    
    sl = length(this_sec_time); %can be any length from 0 to >sr
    
    if sl < sr %missing sample points
        this_sec_time = this_sec_time + sampling_duration * (1:sl)';
        %the rest is NaNs
        this_sec_time = [this_sec_time; NaN(sr-length(this_sec_time),1)];
        
        pupil_position_conv(this_sec_ind_new(1:sl),:) =  pupil_position(this_sec_ind_old,:);
        pupil_size_conv(this_sec_ind_new(1:sl),:) =  pupil_size(this_sec_ind_old,:);
        
    elseif sl >= sr 
        this_sec_time = this_sec_time(1:sr);
        this_sec_time = this_sec_time + sampling_duration * (0:(sr-1))';
        
        pupil_position_conv(this_sec_ind_new(1:sr),:) =  pupil_position(this_sec_ind_old(1:sr),:);
        pupil_size_conv(this_sec_ind_new(1:sr),:) =  pupil_size(this_sec_ind_old(1:sr),:);
    end
    this_sec_time.Format = fmt;
    times_conv.timess(this_sec_ind_new(1:sr),:) = this_sec_time;
end

end