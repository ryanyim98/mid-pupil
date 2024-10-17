    %Adapted from Kret et al.
    % Function for removing isolated sections of data

function valid_id = pupilIslandRemove(RAW,sr,valid_id)


    % 'Sample-islands' are clusters of samples that are temporally seperated
    % from other samples. The minimum distance used to consider
    % samples 'separated' is specified below:
    % Get settings:
    maxSep  = 20; %ms

    
    % When clusters are seperated, they need to have the minimum size
    % specified below. Sample islands smaller than this temporal width and
    % separated from other samples by the distance mentioned above are
    % marked as invalid.
    minIslandWidth = 50; %ms

    % Isolate the usable samples:
    validIndx = find(valid_id);
    mytime= (1:length(valid_id))*(1000/sr); %in ms
    tValidz   = mytime(validIndx)';

    % Return if there are not enough samples:
    if sum(valid_id)<3
        valid_id = valid_id;
        return
    end

    % find segments of missing data
    na_seg = find(isnan(RAW));
    na_ind = [];
    c = 1;
    
    for i = 1:length(na_seg)
        if na_seg(i) == 1 || i == 1
            na_ind(1,1) = na_seg(i);
        elseif na_seg(i) ~= na_seg(i-1) + 1 %increasing
                na_ind(c,2) = na_seg(i-1);
                c = c + 1;
                na_ind(c,1) = na_seg(i);
        end
        if i == length(na_seg)
            na_ind(c,2) = na_seg(i);
        end
    end

    na_ind(:,3) = (na_ind(:,2) - na_ind(:,1) + 1)/sr *1000; %duration of missing chunk in ms

    data_chunk_after = ([na_ind(:,1);NaN]-[NaN;na_ind(:,2)])/sr*1000;%%duration of data in ms
    na_ind(:,4) = data_chunk_after(2:end);

    %if last bit of data chunk is really small, consider the possibility of
    %being an island
    island_candidate_id = find(na_ind(:,3) > maxSep  & na_ind(:,4) <= minIslandWidth);

    for i = 1:length(island_candidate_id)
        thisid = island_candidate_id(i);
        startid = na_ind(thisid,2);
        endid = na_ind(thisid+1,1)-1;
        valid_id(startid:endid) = 0;
    end
       
end