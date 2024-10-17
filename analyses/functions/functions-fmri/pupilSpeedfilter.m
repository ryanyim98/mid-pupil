function RAW = pupilSpeedfilter(RAW,sr) %RAW = input signal, sr= sampling rate
    curDilationSpeeds = diff(RAW)/(1/sr);

    % Generate a two column array with the back and forward dilation speeds:
    backFwdDilations = [[NaN;curDilationSpeeds] [curDilationSpeeds;NaN]];

    maxDilationSpeeds = max(abs(backFwdDilations),[],2);
    madMultiplier = 16; %typically 16

    [med_d, mad, thresh] = madCalc(maxDilationSpeeds,madMultiplier);
    RAW(maxDilationSpeeds > thresh) = NaN;
    
%     valid_id(maxDilationSpeeds > thresh) = 0;
    
end