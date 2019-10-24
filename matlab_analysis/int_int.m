function [cumSUM, SUM, intSUM]=int_int(signal, dt, intervals) 
%signal - the signal to integrate, dt- time between frames in signal,
%intervals - intervals to calculate integral cumSUM, format [start stop]xn
%output: cumSUM -cumulative sum of the integral, SUM - total sum of the
%integral, intSUM - intergral for individual intervals
if isempty(signal)
    cumSUM=nan;
    SUM=nan;
    intSUM=NaN(1, size(intervals,2));
elseif nansum(signal)==0
    cumSUM=nan;
    SUM=nan;
    intSUM=NaN(1, size(intervals,2));
elseif nansum(intervals)==0
    cumSUM=nan;
    SUM=nan;
    intSUM=NaN(1, size(intervals,2));
else
    cumSUM=cumsum(((signal(1:end-1)+signal(2:end))*dt)/2);
    SUM=cumSUM(length(cumSUM));
    for i=1:size(intervals,2)
        intSUM(i)=cumSUM(intervals(2,i))-cumSUM(intervals(1,i));
    end
end
