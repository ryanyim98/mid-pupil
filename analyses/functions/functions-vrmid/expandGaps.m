function [RAW_gapRemoved, exclude_ind] = expandGaps(RAW,sr)
% Function for removing samples around gaps, see 'Standard Settings'
% section of this m file.
%
%--------------------------------------------------------------------------

% Get settings:
minGap      = 75; %ms
maxGap      = 10000; %ms
backPadding = sr*0.2; %time points, ~ 50 ms
fwdPadding  = sr*0.05; %time points

% Blinks produce gaps in the data, the edges of these gaps may feature
% artifacts, as such, dilate gaps:

% Calculate the duration of each gap, and test whether it exceeds the
% thresholds:
if ~isnan(minGap) || ~isnan(maxGap)
    t = find(~isnan(RAW));% t is the time point in 120Hz data (each = 8.3 ms)
    gaps  = diff(t);
    isGapThatNeedsPadding = (gaps * 1000/sr) > minGap & (gaps * 1000/sr) < maxGap;
    isGapThatNeedsExcluding = (gaps * 1000/sr) >= maxGap;
    
    longBlinkStartTimes = t(...
        [isGapThatNeedsExcluding;false]);
    longBlinkEndTimes   = t(...
        [false;isGapThatNeedsExcluding]);
    blinkTime = [];
    
    if ~isempty(longBlinkStartTimes)
        disp("N = "+length(longBlinkStartTimes)+" long blinks detected");
        for i = 1:length(longBlinkStartTimes)
            blinkTime = [blinkTime ...
            longBlinkStartTimes(i):longBlinkEndTimes(i)];
        end
    else
        disp("No long blinks detected;");
    end
    
    blinkTime = blinkTime';
    exclude_ind = ones([length(RAW) 1]);
    exclude_ind(blinkTime) = 0;
    
    % Get the start and end times of each gap:
    gapStartTimes = t(...
        [isGapThatNeedsPadding;false]);
    gapEndTimes   = t(...
        [false;isGapThatNeedsPadding]);

    % Padd gaps that need padding:
    if backPadding>0 || fwdPadding>0
        
        % Detect samples around the gaps:
        isNearGap = any(bsxfun(@gt,t...
            ,gapStartTimes'-backPadding)...
            &bsxfun(@lt,t...
            ,gapEndTimes'+fwdPadding),2);
        
        %         printToConsole(4,'Edge Removal Filter: %i samples
        %         removed.\n'...
        %             ,sum(isNearGap&isValid_In(validIndx)));
        
        % Reject samples too near a gap:
        NearGaptime = t(isNearGap);
        RAW_gapRemoved = RAW;
        RAW_gapRemoved(NearGaptime) = NaN;
    end
    
end

end