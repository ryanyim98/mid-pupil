function valid_id = pupilSpeedfilter(RAW,sr,valid_id)
    curDilationSpeeds = diff(RAW)/(1/sr);

    % Generate a two column array with the back and forward dilation speeds:
    backFwdDilations = [[NaN;curDilationSpeeds] [curDilationSpeeds;NaN]];

    maxDilationSpeeds = max(abs(backFwdDilations),[],2);
    madMultiplier = 16; %16

    [med_d, mad, thresh] = madCalc(maxDilationSpeeds,madMultiplier);

    valid_id(maxDilationSpeeds > thresh) = 0;
    
end