%%
function [med_d, mad, thresh] = madCalc(d,n)
% madCalc calculates the rejection threshold using the mad method.
%
%--------------------------------------------------------------------------
% Calc the median:
med_d  = nanmedian(d);

% Calc the mad:
mad    = nanmedian(abs(d - med_d));

% Calc the threshold:
thresh = med_d + (n * mad);

end