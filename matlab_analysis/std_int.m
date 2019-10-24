function int_std=std_int(signal, intervals)

if isempty(signal)
    int_std=NaN(1, size(intervals,2));
elseif nansum(signal)==0
    int_std=NaN(1, size(intervals,2));
else
    int_std=zeros(1, size(intervals,2));
    for i=1:size(intervals,2)
        int_std(i)=std(signal(intervals(1,i):intervals(2,i)));
    end
end