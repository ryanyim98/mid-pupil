function valid_id = long_blink_mask(valid_id,missing_data) 
%create mask for long blinks
    ind_blinks = missing_data;
    if size(missing_data,2) ~= 3
        disp("Incorrect data dimensions");
    return;
    end
    
    if string(missing_data.Properties.VariableNames(1)) ~= 'start' | ...
            string(missing_data.Properties.VariableNames(2)) ~= 'end'|...
            string(missing_data.Properties.VariableNames(3)) ~= 'time'
        disp("Incorrect variable names");
    end

    ind_blinks = ind_blinks(ind_blinks.time > 10,:);
   
    if size(ind_blinks,1) > 0
       for i = 1:size(ind_blinks,1)
           valid_id(...
            ind_blinks.start(i):ind_blinks.end(i)) = 0;
       end
       disp("Percent removed due to long blinks: " + (1 - mean(valid_id)));
    end
   
%     if missing_data.start(1) == 1
%         valid_id(missing_data.start(1):missing_data.end(1)) = 0;
%         disp("Beginning of recording removed. ");
%     end
%     
%     if missing_data.end(size(missing_data,1)) == length(valid_id)
%         valid_id(missing_data.start(size(missing_data,1)):missing_data.end(size(missing_data,1))) = 0;
%         disp("End of recording removed. ");
%     end
    
end