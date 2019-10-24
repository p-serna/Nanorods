function ind=boarders(crit) %takes sequence of zeros and ones (crit), gives boarders of uninterapted intervals of ones (ind)
k1=1;
k2=1;

if length(crit)>2
    if crit(1)==1
        start=1;
        k1=2;
        if  crit(2)==0
            stop(k2)=1;
            k2=k2+1;

        end
    end
    for i=2:length(crit)-1
        if crit(i)==1
            if crit(i-1)==0
                start(k1)=i;
                k1=k1+1;
            end
            if  crit(i+1)==0
                stop(k2)=i;
                k2=k2+1;

            end
        elseif crit(i)~=0 && crit(i)~=1
            error='not 0 or 1'
            return
        end
        
        
    end
    if crit(length(crit))==1
        stop(k2)=i;
        if crit(length(crit)-1)==0
            start(k1)=i;
        end
    elseif crit(i)~=0 && crit(i)~=1
        error='not 0 or 1'
        return
    end
%     if exist('start') && exist('stop')
%         ind=[start;stop];
%     else
%         ind=[nan;nan];
%     end
elseif length(crit)==2
    if crit(1)==1
        start=1;
        if  crit(2)==0
            stop=1;
        end
        if  crit(2)==1
            stop=2;
        end
    end
    if crit(1)==0
        if  crit(2)==1
            start=2;
            stop=2;
        end
    end
%     if exist('start') && exist('stop')
%         ind=[start;stop];
%     else
%         ind=[nan;nan];
%     end
elseif length(crit)==1
    if crit(1)==1
        start=1;
        stop=1;
    end
end

if exist('start') && exist('stop')
    ind=[start;stop];
else
    ind=[nan;nan];
end

        
            
        