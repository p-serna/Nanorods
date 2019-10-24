function int_mean=mean_int(signal, intervals)

if isempty(signal)
    int_mean=NaN(1, size(intervals,2));
elseif nansum(signal)==0
    int_mean=NaN(1, size(intervals,2));
else
    int_mean=zeros(1, size(intervals,2));
    for i=1:size(intervals,2)
        int_mean(i)=mean(signal(intervals(1,i):intervals(2,i)));
    end
end