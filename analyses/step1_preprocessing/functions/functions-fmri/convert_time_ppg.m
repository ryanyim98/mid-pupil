function [times_conv, ppg_conv] =  convert_time_ppg(times,fmt,sr,ppg) 

secs = unique(times);
times_conv.timess = repelem(secs,sr); %sr Hz strictly
times_conv.time = repelem(secs,sr); %sr Hz strictly
ppg_conv = array2table(NaN([length(times_conv.timess) size(ppg,2)]),'VariableNames',ppg.Properties.VariableNames);

times_conv.timess.Format = fmt;
times_conv.time.Format = "hh:mm:ss";
first_sec_ind = find(times == min(secs));
first_sec_time = times(times == min(secs),:);
this_second = secs(1);
sampling_duration = seconds(1/sr);

%if first second has less than 120 data points, count backwards to
%recover time
if length(first_sec_time) <= sr
    for i = 1:length(first_sec_time)
            first_sec_time(i) = first_sec_time(i) + seconds(1) - sampling_duration *...
                (length(first_sec_time)-i+1);
    end
elseif length(first_sec_time) > sr
    first_sec_ind = first_sec_ind(1:sr);
    for i = 1:sr
            first_sec_time(i) = first_sec_time(i) + seconds(1) - sampling_duration * ...
                (length(first_sec_time) - i+1);
    end
end

ppg_conv(((sr+1)-length(first_sec_ind)):sr,:) =  ppg(first_sec_ind,:);

first_sec_time = [NaN(sr-length(first_sec_time),1); first_sec_time];
first_sec_time.Format = fmt;

times_conv.timess(1:length(first_sec_time)) = first_sec_time;

%move to the next sec
sample_length = [];
for j = 2:length(secs)
    this_sec_time = times(times == secs(j),:); 
    this_sec_ind_old = find(times == secs(j));
    this_sec_ind_new = find(times_conv.time == secs(j));
    
    sl = length(this_sec_time); %can be any length from 0 to >sr
    sample_length = [sample_length,sl];
    
    if sl < sr %missing sample points
        this_sec_time = this_sec_time + sampling_duration * (1:sl)';
        %the rest is NaNs
        this_sec_time = [this_sec_time; NaN(sr-length(this_sec_time),1)];
        
        ppg_conv(this_sec_ind_new(1:sl),:) =  ppg(this_sec_ind_old,:);
        %disp("====sample point < sampling rate====");
        %break
        
    elseif sl >= sr 
        this_sec_time = this_sec_time(1:sr);
        this_sec_time = this_sec_time + sampling_duration * (0:(sr-1))';
        
        ppg_conv(this_sec_ind_new(1:sr),:) =  ppg(this_sec_ind_old(1:sr),:);
    end
    this_sec_time.Format = fmt;
    times_conv.timess(this_sec_ind_new(1:sr),:) = this_sec_time;
end

disp("sample length per sec: "+mean(sample_length));
end