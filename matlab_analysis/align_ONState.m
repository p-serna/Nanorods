function  [Signal_aligned, TimeOnState, IndOnState]=align_ONState(Signal,trh, normalize)
%input: 
%Signal
%trh - threshold for the legth of the onstate period
%normalize - flag for normalization (1 - no normalization, 1 - dF/F, 2-dRt)

%output:
%Signal_aligned - signal with removed blinking and aligned periods of the stimulation square wave
%TimeOnState - on state frames
%IndOnState indices of the on state frames

n=1;
NROI=size(Signal,2);
Nframes=size(Signal,1);

for roi=1:NROI 
    k=1;
    temp=reshape(Signal(:,roi),4,Nframes/4);
    if sum(sum(temp))==0
        temp(temp==0)=NaN;
    end
    temp_OnState=0;
    base=~isnan(mean(temp));
    m=0;
    
    crit=false(1,length(base));
    strt=false(2,length(base));
    for i=1:length(base)        
        if base(i)==1
            m=m+1;
        end
        if i>1 && base(i)==0 && base(i-1)==1
            temp_OnState(k)=m;          
            k=k+1;
            if m>trh-1
                crit(i-m:i-1)=true;
                strt(1,i-m)=true;
                strt(2,i-1)=true;
            end
        end
        if i==length(base) &&  base(i)==1 && base(i-1)==1
            temp_OnState(k)=m;
            k=k+1;
            if m>trh-1
                crit(i-m+1:i)=true;
                strt(1,i-m+1)=true;
                strt(2,i)=true;
            end
        end
        if base(i)==0
            m=0;
        end
    end
    
    if temp_OnState == Nframes/4
        TimeOnState{roi}=NaN;
    else
        TimeOnState{roi}=temp_OnState;
    end
    if nansum(nansum(temp))==0
         TimeOnState{roi}=NaN;
    end

    frames_start=find(strt(1,:));
    frames_stop=find(strt(2,:));
   
    temp=temp.*crit;
    temp(temp==0)=NaN;
    if normalize==1
        temp=temp-mean(temp(3:4,:));
    end
    
    if normalize==2
        temp=(temp-mean(temp(3:4,:)))./mean(temp(3:4,:));
    end
    
    temp=reshape(temp, 1, Nframes);
    temp(isnan(temp))=[];
    Signal_aligned{roi}=temp;
    IndOnState{roi}=[frames_start; frames_stop];
end

    

