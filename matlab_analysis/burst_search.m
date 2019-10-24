
%% User-defined variables
roi=23; %roi for the interval search
trh=1; %length of the minimum ON interval (in periods of the wave, x4 in frames)
normalizeF=2; %normalization flag
normalizeR=1;
wind=10; %size of moving mean window for score calculation (in periods)
%% Aligning and normalization of the signal (F and ratio)

[sig_aligned_wave, ~, ~]=align_ONState(Signal_wave_2ch_OnState(:,roi), trh, normalizeF);
[sig_aligned_contrl, ~, ~]=align_ONState(Signal_ctrl_2ch_OnState(:,roi), trh, normalizeF);
[ratio_aligned_wave, ~, ~]=align_ONState(Ratio_wave(:,roi), trh, normalizeR);
[ratio_aligned_contrl, ~, ~]=align_ONState(Ratio_ctrl(:,roi), trh, normalizeR);

sig_aligned_wave=sig_aligned_wave{1};
sig_dF_F_wave=reshape(sig_aligned_wave, 4, length(sig_aligned_wave)/4);
sig_dF_F_wave=mean(sig_dF_F_wave(1:2,:));

sig_aligned_contrl=sig_aligned_contrl{1};
sig_dF_F_contrl=reshape(sig_aligned_contrl, 4, length(sig_aligned_contrl)/4);
sig_dF_F_contrl=mean(sig_dF_F_contrl(1:2,:));

ratio_aligned_wave=ratio_aligned_wave{1};
sig_dRt_wave=reshape(ratio_aligned_wave, 4, length(ratio_aligned_wave)/4);
sig_dRt_wave=mean(sig_dRt_wave(1:2,:));

ratio_aligned_contrl=ratio_aligned_contrl{1};
sig_dRt_contrl=reshape(ratio_aligned_contrl, 4, length(ratio_aligned_contrl)/4);
sig_dRt_contrl=mean(sig_dRt_contrl(1:2,:));

%% Calculation of the burst score
%F
Nframes=size(sig_aligned_wave,2);
sig=reshape(sig_aligned_wave,4,Nframes/4);
sig=mean(sig(1:2,:));

score_F_wave=running_mean(sig',wind); %burst score wave

Nframes=size(sig_aligned_contrl,2);
sig=reshape(sig_aligned_contrl,4,Nframes/4);
sig=mean(sig(1:2,:));

score_F_contrl=running_mean(sig',wind); %burst score control

%Ratio
Nframes=size(ratio_aligned_wave,2);
sig=reshape(ratio_aligned_wave,4,Nframes/4);
sig=mean(sig(1:2,:));

score_R_wave=running_mean(sig',wind);

Nframes=size(ratio_aligned_contrl,2);
sig=reshape(ratio_aligned_contrl,4,Nframes/4);
sig=mean(sig(1:2,:));

score_R_contrl=running_mean(sig',wind);


clearvars sig Nframes ind trhF trhR normalizeF normalizeR trh

%% getting indices of positive and negatie score intervals

indices_F_wave_pos=boarders(score_F_wave>0);
indices_F_wave_neg=boarders(score_F_wave<0);
indices_F_contrl_pos=boarders(score_F_contrl>0);
indices_F_contrl_neg=boarders(score_F_contrl<0);

indices_R_wave_pos=boarders(score_R_wave>0);
indices_R_wave_neg=boarders(score_R_wave<0);
indices_R_contrl_pos=boarders(score_R_contrl>0);
indices_R_contrl_neg=boarders(score_R_contrl<0);



%% Calculating integrated score for positive and negative intervals and the sum of the integrated score for each ROI

clearvars SUM_F_wave_pos SUM_F_wave_neg SUM_F_contrl_pos SUM_F_wave_neg SUM_R_wave_pos SUM_R_wave_neg SUM_R_contrl_pos SUM_R_wave_neg
[~, ~,SUM_F_wave_pos]=int_int(score_F_wave, 1, indices_F_wave_pos);
[~, ~,SUM_F_wave_neg]=int_int(score_F_wave, 1, indices_F_wave_neg);
[~, ~,SUM_F_contrl_pos]=int_int(score_F_contrl, 1, indices_F_contrl_pos);
[~, ~,SUM_F_contrl_neg]=int_int(score_F_contrl, 1, indices_F_contrl_neg);
[~, ~,SUM_R_wave_pos]=int_int(score_R_wave, 1, indices_R_wave_pos);
[~, ~,SUM_R_wave_neg]=int_int(score_R_wave, 1, indices_R_wave_neg);
[~, ~,SUM_R_contrl_pos]=int_int(score_R_contrl, 1, indices_R_contrl_pos);
[~, ~,SUM_R_contrl_neg]=int_int(score_R_contrl, 1, indices_R_contrl_neg);
% 
%% select best intervals; intervals dimentions: start, stop (in cycles, x4 in frames), burst score, mean dF/F)
trhF=mean([SUM_F_contrl_neg, SUM_F_contrl_pos])+std([SUM_F_contrl_neg, SUM_F_contrl_pos])*3;%
trhR=mean([SUM_R_contrl_neg, SUM_R_contrl_pos])+std([SUM_R_contrl_neg, SUM_R_contrl_pos])*3;%

intervals_F_wave_selected=[indices_F_wave_pos(:,SUM_F_wave_pos>trhF), indices_F_wave_neg(:,SUM_F_wave_neg<trhF*(-1))];
intervals_F_wave_selected(3,:)=[SUM_F_wave_pos(SUM_F_wave_pos>trhF), SUM_F_wave_neg(SUM_F_wave_neg<trhF*(-1))];
intervals_F_wave_selected(4,:)=[mean_int(sig_dF_F_wave, indices_F_wave_pos(:,SUM_F_wave_pos>trhF)), mean_int(sig_dF_F_wave, indices_F_wave_neg(:,SUM_F_wave_neg<trhF*(-1)))];
intervals_F_wave_selected(5,:)=[std_int(sig_dF_F_wave, indices_F_wave_pos(:,SUM_F_wave_pos>trhF)), std_int(sig_dF_F_wave, indices_F_wave_neg(:,SUM_F_wave_neg<trhF*(-1)))];

intervals_F_contrl_selected=[indices_F_contrl_pos(:,SUM_F_contrl_pos>trhF), indices_F_contrl_neg(:,SUM_F_contrl_neg<trhF*(-1))];
intervals_F_contrl_selected(3,:)=[SUM_F_contrl_pos(SUM_F_contrl_pos>trhF), SUM_F_contrl_neg(SUM_F_contrl_neg<trhF*(-1))];
intervals_F_contrl_selected(4,:)=[mean_int(sig_dF_F_contrl, indices_F_contrl_pos(:,SUM_F_contrl_pos>trhF)), mean_int(sig_dF_F_contrl, indices_F_contrl_neg(:,SUM_F_contrl_neg<trhF*(-1)))];
intervals_F_contrl_selected(5,:)=[std_int(sig_dF_F_contrl, indices_F_contrl_pos(:,SUM_F_contrl_pos>trhF)), std_int(sig_dF_F_contrl, indices_F_contrl_neg(:,SUM_F_contrl_neg<trhF*(-1)))];

intervals_R_wave_selected=[indices_R_wave_pos(:,SUM_R_wave_pos>trhR), indices_R_wave_neg(:,SUM_R_wave_neg<trhR*(-1))];
intervals_R_wave_selected(3,:)=[SUM_R_wave_pos(SUM_R_wave_pos>trhR), SUM_R_wave_neg(SUM_R_wave_neg<trhR*(-1))];
intervals_R_wave_selected(4,:)=[mean_int(sig_dRt_wave, indices_R_wave_pos(:,SUM_R_wave_pos>trhR)), mean_int(sig_dRt_wave, indices_R_wave_neg(:,SUM_R_wave_neg<trhR*(-1)))];
intervals_R_wave_selected(5,:)=[std_int(sig_dRt_wave, indices_R_wave_pos(:,SUM_R_wave_pos>trhR)), std_int(sig_dRt_wave, indices_R_wave_neg(:,SUM_R_wave_neg<trhR*(-1)))];

intervals_R_contrl_selected=[indices_R_contrl_pos(:,SUM_R_contrl_pos>trhR), indices_R_contrl_neg(:,SUM_R_contrl_neg<trhR*(-1))];
intervals_R_contrl_selected(3,:)=[SUM_R_contrl_pos(SUM_R_contrl_pos>trhR), SUM_R_contrl_neg(SUM_R_contrl_neg<trhR*(-1))];
intervals_R_contrl_selected(4,:)=[mean_int(sig_dRt_contrl, indices_R_contrl_pos(:,SUM_R_contrl_pos>trhR)), mean_int(sig_dRt_contrl, indices_R_contrl_neg(:,SUM_R_contrl_neg<trhR*(-1)))];
intervals_R_contrl_selected(5,:)=[std_int(sig_dRt_contrl, indices_R_contrl_pos(:,SUM_R_contrl_pos>trhR)), std_int(sig_dRt_contrl, indices_R_contrl_neg(:,SUM_R_contrl_neg<trhR*(-1)))];

%% calculating mean dF/F and dRt for selected intervals
intervals_F_selected_stat(1)=nanmean(intervals_F_wave_selected(4,:));
intervals_F_selected_stat(2)=nanstd(intervals_F_wave_selected(4,:));
intervals_F_selected_stat(3)=nanmean(intervals_F_contrl_selected(4,:));
intervals_F_selected_stat(4)=nanstd(intervals_F_contrl_selected(4,:));
[~, intervals_F_selected_stat(5)]=ttest2(intervals_F_wave_selected(4,:),intervals_F_contrl_selected(4,:));

intervals_R_selected_stat(1)=nanmean(intervals_R_wave_selected(4,:));
intervals_R_selected_stat(2)=nanstd(intervals_R_wave_selected(4,:));
intervals_R_selected_stat(3)=nanmean(intervals_R_contrl_selected(4,:));
intervals_R_selected_stat(4)=nanstd(intervals_R_contrl_selected(4,:));
[~, intervals_R_selected_stat(5)]=ttest2(intervals_R_wave_selected(4,:),intervals_R_contrl_selected(4,:));

%% plotting dF/F trace with selected intervals

f_wave=figure;
f_wave.Position=[300, 150, 1800, 1000];

%plotting dF/F(t) and burst score
%WAVE
sig=reshape(sig_aligned_wave,4,length(sig_aligned_wave)/4);
sig=[mean(sig(1:2,:)); mean(sig(3:4,:))];
sig=sig(:);

ax1=subplot(2,1,1);
plot(sig(:),'.-', 'Color', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.3 0.3 0.3]) 
hold on

plot(1:2:length(score_F_wave)*2,score_F_wave, 'r-')
hold on
axis('tight');

%plotting dF/F(t) and burst score
%CONTROL
sig=reshape(sig_aligned_contrl,4,length(sig_aligned_contrl)/4);
sig=[mean(sig(1:2,:)); mean(sig(3:4,:))];
sig=sig(:);

ax2=subplot(2,1,2);
plot(sig(:),'.-', 'Color', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.3 0.3 0.3]) 
hold on

plot(1:2:length(score_F_contrl)*2,score_F_contrl, 'r-')
hold on
axis('tight');

%axis limits
ax1.XLim(1)=min(ax1.XLim(1), ax2.XLim(1));
ax2.XLim(1)=min(ax1.XLim(1), ax2.XLim(1));
ax1.XLim(2)=max(ax1.XLim(2), ax2.XLim(2));
ax2.XLim(2)=max(ax1.XLim(2), ax2.XLim(2));

ax1.YLim(1)=min(ax1.YLim(1), ax2.YLim(1));
ax2.YLim(1)=min(ax1.YLim(1), ax2.YLim(1));
ax1.YLim(2)=max(ax1.YLim(2), ax2.YLim(2));
ax2.YLim(2)=max(ax1.YLim(2), ax2.YLim(2));

%plotting intervals
%WAVE
for i=1:size(intervals_F_wave_selected,2)
    x1=intervals_F_wave_selected(1,i)*2-1;
    x2=intervals_F_wave_selected(2,i)*2;
    y1=ax1.YLim(1);
    y2=ax1.YLim(2);
    x = [x1, x2, x2, x1, x1];
    y = [y1, y1, y2, y2, y1];
    subplot(2,1,1)
    plot(x, y, 'r--');
    t=text(x1+(x2-x1)/2*0.1,y1+(y2-y1)*0.1,num2str(i));
    t.FontSize=14;
    t.Color='r';
    hold on
end
legend('dF/F(t)','burst score','selected intervals')

%plotting intervals
%CONTROL
for i=1:size(intervals_F_contrl_selected,2)
    x1=intervals_F_contrl_selected(1,i)*2-1;
    x2=intervals_F_contrl_selected(2,i)*2;
    y1=ax2.YLim(1);
    y2=ax2.YLim(2);
    x = [x1, x2, x2, x1, x1];
    y = [y1, y1, y2, y2, y1];
    subplot(2,1,2)
    plot(x, y, 'r--');
    t=text(x1+(x2-x1)/2*0.1,y1+(y2-y1)*0.1,num2str(i));
    t.FontSize=14;
    t.Color='r';  
    hold on
end
legend('dF/F(t)','burst score','selected intervals')

%plotting dashed lines at 0
subplot(2,1,1)
plot([ax2.XLim(1) ax2.XLim(2)], [0 0], 'k--')
xt=xticks;
xt={xt*0.02};
xticklabels(xt)
hold on

subplot(2,1,2)
plot([ax2.XLim(1) ax2.XLim(2)], [0 0], 'k--')
xt=xticks;
xt={xt*0.02};
xticklabels(xt)
hold on

%plotting color-coded F(t)
%WAVE
mn=min([score_F_wave' score_F_contrl']);
mx=max([score_F_wave' score_F_contrl']);
clr=64*(score_F_wave-mn+1/64)/(mx-mn+1/64);
clr(clr>64)=64;
clr(clr<1)=1;
map=parula;

sig=reshape(sig_aligned_wave,4,length(sig_aligned_wave)/4);
sig=[mean(sig(1:2,:)); mean(sig(3:4,:))];
sig=sig(:);

for i=1:length(score_F_wave)
    subplot(2,1,1)
    plot(2*i-1,sig(2*i-1), 'o', 'Color', [map(ceil(clr(i)),:)])
    hold on
end

title(['dF/F trace ROI', num2str(roi), ' with voltage modulation; int score, mean dF/F of the intervals is ', num2str(round(mean(intervals_F_wave_selected(4,:)),2))]);
set(gca,'FontSize',20)
xlabel('time(s)')
ylabel('dF/F')
xt=xticks;
xt={xt*0.02};
xticklabels(xt)

%colorbar properties
t=mn; %min(temp); %Burst score lower limit
m=mx; %max(temp); 
n=5; %number of thicks
c=colorbar;
c.Ticks = linspace(0, 1, n);
c.TickLabels = round(linspace(t, m, n),2);
c.FontSize=15;


%plotting color-coded F(t)
%CONTROL

clr=64*(score_F_contrl-mn+1/64)/(mx-mn+1/64);
clr(clr>64)=64;
clr(clr<1)=1;


sig=reshape(sig_aligned_contrl,4,length(sig_aligned_contrl)/4);
sig=[mean(sig(1:2,:)); mean(sig(3:4,:))];
sig=sig(:);

for i=1:length(score_F_contrl)
    subplot(2,1,2)
    plot(2*i-1,sig(2*i-1), 'o', 'Color', [map(ceil(clr(i)),:)])
    hold on
end

title(['dF/F trace ROI', num2str(roi), ' control; int score, mean dF/F of the intervals is ', num2str(round(mean(intervals_F_contrl_selected(4,:)),2))]);
set(gca,'FontSize',20)
xlabel('time(s)')
ylabel('dF/F')

%colorbar properties
t=mn; %min(temp); %Burst score lower limit
m=mx; %max(temp); 
n=5; %number of thicks
c=colorbar;
c.Ticks = linspace(0, 1, n);
c.TickLabels = round(linspace(t, m, n),2);
c.FontSize=15;

set(gcf,'renderer','Painters')

clearvars f c clr i m n map sig start stop t x x1 x2 y y1 y2 r j temp mn mx xt  ax1 ax2 

%% plotting dRt trace with selected intervals

f_wave=figure;
f_wave.Position=[300, 150, 1800, 1000];

%plotting dRt(t) and burst score
%WAVE
sig=reshape(ratio_aligned_wave,4,length(ratio_aligned_wave)/4);
sig=[mean(sig(1:2,:)); mean(sig(3:4,:))];
sig=sig(:);

ax1=subplot(2,1,1);
plot(sig(:),'.-', 'Color', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.3 0.3 0.3]) 
hold on

plot(1:2:length(score_R_wave)*2,score_R_wave, 'r-')
hold on
axis('tight');

%plotting dRt(t) and burst score
%CONTROL
sig=reshape(ratio_aligned_contrl,4,length(ratio_aligned_contrl)/4);
sig=[mean(sig(1:2,:)); mean(sig(3:4,:))];
sig=sig(:);

ax2=subplot(2,1,2);
plot(sig(:),'.-', 'Color', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.3 0.3 0.3]) 
hold on

plot(1:2:length(score_R_contrl)*2,score_R_contrl, 'r-')
hold on
axis('tight');

%axis limits
ax1.XLim(1)=min(ax1.XLim(1), ax2.XLim(1));
ax2.XLim(1)=min(ax1.XLim(1), ax2.XLim(1));
ax1.XLim(2)=max(ax1.XLim(2), ax2.XLim(2));
ax2.XLim(2)=max(ax1.XLim(2), ax2.XLim(2));

ax1.YLim(1)=min(ax1.YLim(1), ax2.YLim(1));
ax2.YLim(1)=min(ax1.YLim(1), ax2.YLim(1));
ax1.YLim(2)=max(ax1.YLim(2), ax2.YLim(2));
ax2.YLim(2)=max(ax1.YLim(2), ax2.YLim(2));

%plotting intervals
%WAVE
for i=1:size(intervals_R_wave_selected,2)
    x1=intervals_R_wave_selected(1,i)*2-1;
    x2=intervals_R_wave_selected(2,i)*2;
    y1=ax1.YLim(1);
    y2=ax1.YLim(2);
    x = [x1, x2, x2, x1, x1];
    y = [y1, y1, y2, y2, y1];
    subplot(2,1,1)
    plot(x, y, 'r--');
    t=text(x1+(x2-x1)/2*0.1,y1+(y2-y1)*0.1,num2str(i));
    t.FontSize=14;
    t.Color='r';
    hold on
end
subplot(2,1,1)
legend('dRt(t)','burst score','selected intervals')

%plotting intervals
%CONTROL
for i=1:size(intervals_R_contrl_selected,2)
    x1=intervals_R_contrl_selected(1,i)*2-1;
    x2=intervals_R_contrl_selected(2,i)*2;
    y1=ax2.YLim(1);
    y2=ax2.YLim(2);
    x = [x1, x2, x2, x1, x1];
    y = [y1, y1, y2, y2, y1];
    subplot(2,1,2)
    plot(x, y, 'r--');
    t=text(x1+(x2-x1)/2*0.1,y1+(y2-y1)*0.1,num2str(i));
    t.FontSize=14;
    t.Color='r';  
    hold on
end
subplot(2,1,2)
legend('dRt(t)','burst score','selected intervals')

%plotting dashed lines at 0
subplot(2,1,1)
plot([ax2.XLim(1) ax2.XLim(2)], [0 0], 'k--')
xt=xticks;
xt={xt*0.02};
xticklabels(xt)
hold on

subplot(2,1,2)
plot([ax2.XLim(1) ax2.XLim(2)], [0 0], 'k--')
xt=xticks;
xt={xt*0.02};
xticklabels(xt)
hold on

%plotting color-coded R(t)
%WAVE
mn=min([score_R_wave' score_R_contrl']);
mx=max([score_R_wave' score_R_contrl']);
clr=64*(score_R_wave-mn+1/64)/(mx-mn+1/64);
clr(clr>64)=64;
clr(clr<1)=1;
map=parula;

sig=reshape(ratio_aligned_wave,4,length(ratio_aligned_wave)/4);
sig=[mean(sig(1:2,:)); mean(sig(3:4,:))];
sig=sig(:);

for i=1:length(score_R_wave)
    subplot(2,1,1)
    plot(2*i-1,sig(2*i-1), 'o', 'Color', [map(ceil(clr(i)),:)])
    hold on
end

title(['dRt trace ROI', num2str(roi), ' with voltage modulation; int score, mean dRt of the intervals is ', num2str(round(mean(intervals_R_wave_selected(4,:)),2))]);
set(gca,'FontSize',20)
xlabel('time(s)')
ylabel('dRt')
xt=xticks;
xt={xt*0.02};
xticklabels(xt)

%colorbar properties
t=mn; %min(temp); %Burst score lower limit
m=mx; %max(temp); 
n=5; %number of thicks
c=colorbar;
c.Ticks = linspace(0, 1, n);
c.TickLabels = round(linspace(t, m, n),2);
c.FontSize=15;


%plotting color-coded R(t)
%CONTROL

clr=64*(score_R_contrl-mn+1/64)/(mx-mn+1/64);
clr(clr>64)=64;
clr(clr<1)=1;


sig=reshape(ratio_aligned_contrl,4,length(ratio_aligned_contrl)/4);
sig=[mean(sig(1:2,:)); mean(sig(3:4,:))];
sig=sig(:);

for i=1:length(score_R_contrl)
    subplot(2,1,2)
    plot(2*i-1,sig(2*i-1), 'o', 'Color', [map(ceil(clr(i)),:)])
    hold on
end

title(['dRt trace ROI', num2str(roi), ' control; int score, mean dRt of the intervals is ', num2str(round(mean(intervals_R_contrl_selected(4,:)),2))]);
set(gca,'FontSize',20)
xlabel('time(s)')
ylabel('dRt')

%colorbar properties
t=mn; %min(temp); %Burst score lower limit
m=mx; %max(temp); 
n=5; %number of thicks
c=colorbar;
c.Ticks = linspace(0, 1, n);
c.TickLabels = round(linspace(t, m, n),2);
c.FontSize=15;

set(gcf,'renderer','Painters')

clearvars f1 f2 c clr i m n map sig start stop t x x1 x2 y y1 y2 r j temp mn mx xt ax1 ax2 


%% plotting individual intervals dF
%WAVE
sig=reshape(sig_aligned_wave,4,length(sig_aligned_wave)/4);
sig=[mean(sig(1:2,:)); mean(sig(3:4,:))];
sig=sig(:);

for i=1:size(intervals_F_wave_selected,2)
    int=intervals_F_wave_selected(1,i)*2-1:intervals_F_wave_selected(2,i)*2;
    f_wave(i)=figure;
    ax_wave(i)=subplot(1,1,1);
    plot(int, sig(int),'.-', 'Color', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.3 0.3 0.3])
    axis ('tight')
    xlabel('time(s)')
    ylabel('dF/F')
    title(['ROI', num2str(roi), ' with v.m., interval ', num2str(i), ...
        ', mean dF/F=', num2str(round(intervals_F_wave_selected(4,i),2)), '\pm', ...
        num2str(round(intervals_F_wave_selected(5,i)/sqrt(sum(~isnan(intervals_F_wave_selected(1,i):intervals_F_wave_selected(2,i)))),2))])
    xt=xticks;
    xt={xt*0.02};
    xticklabels(xt)
    hold on
    for j=1:2:length(int)
        if sig(int(j))>0
            plot(int(j),sig(int(j)),'ro')
        else
            plot(int(j),sig(int(j)),'bo')
        end
        hold on
    end
     %legend('dF/F(t)', 'dF/F>0', 'dF/F<0')
end


clearvars mXLim mYLim i j int sig xt

    
%CONTROL
sig=reshape(sig_aligned_contrl,4,length(sig_aligned_contrl)/4);
sig=[mean(sig(1:2,:)); mean(sig(3:4,:))];
sig=sig(:);

for i=1:size(intervals_F_contrl_selected,2)
    int=intervals_F_contrl_selected(1,i)*2-1:intervals_F_contrl_selected(2,i)*2;
    f_contrl(i)=figure;
    ax_contrl(i)=subplot(1,1,1);
    plot(int, sig(int),'.-', 'Color', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.3 0.3 0.3])
    axis ('tight')
    xlabel('time(s)')
    ylabel('dF/F')
    title(['ROI', num2str(roi), ' control, interval ', num2str(i), ...
        ', mean dF/F=', num2str(round(intervals_F_contrl_selected(4,i),2)), '\pm', ...
        num2str(round(intervals_F_contrl_selected(5,i)/sqrt(sum(~isnan(intervals_F_contrl_selected(1,i):intervals_F_contrl_selected(2,i)))),2))])
    xt=xticks;
    xt={xt*0.02};
    xticklabels(xt)
    hold on
    for j=1:2:length(int)
        if sig(int(j))>0
            plot(int(j),sig(int(j)),'ro')
        else
            plot(int(j),sig(int(j)),'bo')
        end
        hold on
    end
     %legend('dF/F(t)', 'dF/F>0', 'dF/F<0')
end


mYLim(1)=ax_wave(1).YLim(1);
mYLim(2)=ax_wave(1).YLim(2);
for i=2:size(intervals_F_wave_selected,2)
    mYLim(1)=min(ax_wave(i).YLim(1),mYLim(1));
    mYLim(2)=max(ax_wave(i).YLim(2),mYLim(2));
end

for i=1:size(intervals_F_contrl_selected,2)
    mYLim(1)=min(ax_contrl(i).YLim(1),mYLim(1));
    mYLim(2)=max(ax_contrl(i).YLim(2),mYLim(2));
end


for i=1:size(intervals_F_wave_selected,2)
    ax_wave(i).YLim(1)=mYLim(1);
    ax_wave(i).YLim(2)=mYLim(2);
    plot(ax_wave(i), [ax_wave(i).XLim(1) ax_wave(i).XLim(2)], [0 0], 'k--')
end

for i=1:size(intervals_F_contrl_selected,2)
    ax_contrl(i).YLim(1)=mYLim(1);
    ax_contrl(i).YLim(2)=mYLim(2);
    plot(ax_contrl(i), [ax_contrl(i).XLim(1) ax_contrl(i).XLim(2)], [0 0], 'k--')
end 

clearvars f_wave f_contrl ax_wave ax_contrl mXLim mYLim i j int sig xt    
    
%% plotting individual intervals dRt
%WAVE
sig=reshape(ratio_aligned_wave,4,length(ratio_aligned_wave)/4);
sig=[mean(sig(1:2,:)); mean(sig(3:4,:))];
sig=sig(:);

for i=1:size(intervals_R_wave_selected,2)
    int=intervals_R_wave_selected(1,i)*2-1:intervals_R_wave_selected(2,i)*2;
    f_wave(i)=figure;
    ax_wave(i)=subplot(1,1,1);
    plot(int, sig(int),'.-', 'Color', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.3 0.3 0.3])
    axis ('tight')
    xlabel('time(s)')
    ylabel('dRt')
    title(['ROI', num2str(roi), ' with v.m., interval ', num2str(i), ...
        ', mean dRt=', num2str(round(intervals_R_wave_selected(4,i),2)), '\pm', ...
        num2str(round(intervals_R_wave_selected(5,i)/sqrt(sum(~isnan(intervals_R_wave_selected(1,i):intervals_R_wave_selected(2,i)))),2))])
    xt=xticks;
    xt={xt*0.02};
    xticklabels(xt)
    hold on
    for j=1:2:length(int)
        if sig(int(j))>0
            plot(int(j),sig(int(j)),'ro')
        else
            plot(int(j),sig(int(j)),'bo')
        end
        hold on
    end
     %legend('dRt(t)', 'dRt>0', 'dRt<0')
end


clearvars mXLim mYLim i j int sig xt

    
%CONTROL
sig=reshape(ratio_aligned_contrl,4,length(ratio_aligned_contrl)/4);
sig=[mean(sig(1:2,:)); mean(sig(3:4,:))];
sig=sig(:);

for i=1:size(intervals_R_contrl_selected,2)
    int=intervals_R_contrl_selected(1,i)*2-1:intervals_R_contrl_selected(2,i)*2;
    f_contrl(i)=figure;
    ax_contrl(i)=subplot(1,1,1);
    plot(int, sig(int),'.-', 'Color', [0.5 0.5 0.5], 'MarkerEdgeColor', [0.3 0.3 0.3])
    axis ('tight')
    xlabel('time(s)')
    ylabel('dRt')
    title(['ROI', num2str(roi), ' control, interval ', num2str(i), ...
        ', mean dRt=', num2str(round(intervals_R_contrl_selected(4,i),2)), '\pm', ...
        num2str(round(intervals_R_contrl_selected(5,i)/sqrt(sum(~isnan(intervals_R_contrl_selected(1,i):intervals_R_contrl_selected(2,i)))),2))])
    xt=xticks;
    xt={xt*0.02};
    xticklabels(xt)
    hold on
    for j=1:2:length(int)
        if sig(int(j))>0
            plot(int(j),sig(int(j)),'ro')
        else
            plot(int(j),sig(int(j)),'bo')
        end
        hold on
    end
     %legend('dRt(t)', 'dRt>0', 'dRt<0')
end


mYLim(1)=ax_wave(1).YLim(1);
mYLim(2)=ax_wave(1).YLim(2);
for i=2:size(intervals_R_wave_selected,2)
    mYLim(1)=min(ax_wave(i).YLim(1),mYLim(1));
    mYLim(2)=max(ax_wave(i).YLim(2),mYLim(2));
end

for i=1:size(intervals_R_contrl_selected,2)
    mYLim(1)=min(ax_contrl(i).YLim(1),mYLim(1));
    mYLim(2)=max(ax_contrl(i).YLim(2),mYLim(2));
end


for i=1:size(intervals_R_wave_selected,2)
    ax_wave(i).YLim(1)=mYLim(1);
    ax_wave(i).YLim(2)=mYLim(2);
    plot(ax_wave(i), [ax_wave(i).XLim(1) ax_wave(i).XLim(2)], [0 0], 'k--')
end

for i=1:size(intervals_R_contrl_selected,2)
    ax_contrl(i).YLim(1)=mYLim(1);
    ax_contrl(i).YLim(2)=mYLim(2);
    plot(ax_contrl(i), [ax_contrl(i).XLim(1) ax_contrl(i).XLim(2)], [0 0], 'k--')
end 


clearvars f_wave f_contrl ax_wave ax_contrl mXLim mYLim i j int sig xt   



%%  histograms of dF_F and dRt

f1=figure;
h1=histogram(sig_dF_F_wave, 20);
title(['dF/F of roi', num2str(roi)])
hold on
h2=histogram (sig_dF_F_contrl, h1.BinEdges);
legend ('with v.m.', 'control')
h1.DisplayStyle=('stairs');
h2.DisplayStyle=('stairs');
h1.EdgeColor='red';
h2.EdgeColor='blue';
xlabel('dF/F')
ylabel('frequency')

f2=figure;
plot(h1.BinLimits(1)+h1.BinWidth/2:h1.BinWidth:h1.BinLimits(2), cumsum(h1.BinCounts)/sum(h1.BinCounts), 'r-');
title(['dF/F cummulative probability of roi', num2str(roi)])
hold on
plot(h2.BinLimits(1)+h2.BinWidth/2:h2.BinWidth:h1.BinLimits(2), cumsum(h2.BinCounts)/sum(h2.BinCounts), 'b-');
legend ('with v.m.', 'control')
xlabel('dF/F')
ylabel('probability')

f3=figure;
h1=histogram(sig_dRt_wave, 20);
title(['dRt of roi', num2str(roi)])
hold on
h2=histogram (sig_dRt_contrl, h1.BinEdges);
legend ('with v.m.', 'control')
h1.DisplayStyle=('stairs');
h2.DisplayStyle=('stairs');
h1.EdgeColor='red';
h2.EdgeColor='blue';
xlabel('dRt')
ylabel('frequency')

f4=figure;
plot(h1.BinLimits(1)+h1.BinWidth/2:h1.BinWidth:h1.BinLimits(2), cumsum(h1.BinCounts)/sum(h1.BinCounts), 'r-');
title(['dRt cummulative probability of roi', num2str(roi)])
hold on
plot(h2.BinLimits(1)+h2.BinWidth/2:h2.BinWidth:h1.BinLimits(2), cumsum(h2.BinCounts)/sum(h2.BinCounts), 'b-');
legend ('with v.m.', 'control')
xlabel('dRt')
ylabel('probability')

clearvars h1 h2 f1 f2
%% plotting dF/F and dR for selected intervals
sig=[intervals_F_selected_stat(1),intervals_F_selected_stat(3)];
sig(isnan(sig))=0;
stdev=[intervals_F_selected_stat(2),intervals_F_selected_stat(4)];
stdev(isnan(stdev))=0;
temp1=intervals_F_wave_selected(4,:);
temp2=intervals_F_contrl_selected(4,:);


f1=figure('position',[100,200,600,1000]);subplot(2,1,1)
b=bar(sig);
hold on
errorbar([1 2],[mean(temp1), mean(temp2)], [std(temp1)/sqrt(length(temp1)),std(temp2)/sqrt(length(temp2))], '.');

title(['mean dF/F \pm ste of selected intervals, ROI', num2str(roi), ', p=', num2str(round((intervals_F_selected_stat(5)),3))])
ylabel('dF/F')
xticklabels({'with v.m.', 'control'})

subplot(2,1,2)
scatter(repmat(1, length(temp1),1) ,temp1, '.', 'MarkerEdgeColor', [0.7 0.7 0.7])
hold on
e=errorbar(1,mean(temp1), std(temp1)/sqrt(length(temp1)));
e.Marker='o';
e.MarkerFaceColor=[0.7 0.7 0.7];
e.MarkerEdgeColor=[0 0 0];
e.Color=[0 0 0];
hold on
scatter(repmat(2, length(temp2),1) ,temp2, '.', 'MarkerEdgeColor', [0.3 0.3 0.3])
hold on
e=errorbar(2,mean(temp2), std(temp2)/sqrt(length(temp2)));
e.Marker='o';
e.MarkerFaceColor=[0.3 0.3 0.3];
e.MarkerEdgeColor=[0 0 0];
e.Color=[0 0 0];
hold on
plot([0.5 2.5], [0 0], 'k--')
xlim([0.5 2.5])
l=ylim;
xticks([1 2])
xticklabels({'with v.m.','control'})
%yticks(l(1):(l(2)-l(1))/10:l(2))
title(['mean dF/F \pm ste of selected intervals, ROI', num2str(roi), ', p=', num2str(round((intervals_F_selected_stat(5)),3))])
ylabel('dF/F')
text(1.1, mean(temp1), ['dF/F=',num2str(round(mean(temp1),2)), '\pm', num2str(round(std(temp1)/sqrt(length(temp1)),2))])
text(2.1, mean(temp2), ['dF/F=',num2str(round(mean(temp2),2)), '\pm', num2str(round(std(temp2)/sqrt(length(temp2)),2))])

clearvars temp temp1 temp2 e b l
 

sig=[intervals_R_selected_stat(1),intervals_R_selected_stat(3)];
sig(isnan(sig))=0;
stdev=[intervals_R_selected_stat(2),intervals_R_selected_stat(4)];
stdev(isnan(stdev))=0;
temp1=intervals_R_wave_selected(4,:);
temp2=intervals_R_contrl_selected(4,:);

f2=figure('position',[100,200,600,1000]);
subplot(2,1,1)
b=bar(sig);
hold on
errorbar([1 2],[mean(temp1), mean(temp2)], [std(temp1)/sqrt(length(temp1)),std(temp2)/sqrt(length(temp2))], '.');
title(['mean dRt \pm ste of selected intervals, ROI', num2str(roi), ', p=', num2str(round((intervals_R_selected_stat(5)),3))])
ylabel('dRt')
xticklabels({'with v.m.', 'control'})

subplot(2,1,2)
scatter(repmat(1, length(temp1),1) ,temp1, '.', 'MarkerEdgeColor', [0.7 0.7 0.7])
hold on
e=errorbar(1,mean(temp1), std(temp1)/sqrt(length(temp1)));
e.Marker='o';
e.MarkerFaceColor=[0.7 0.7 0.7];
e.MarkerEdgeColor=[0 0 0];
e.Color=[0 0 0];
hold on
scatter(repmat(2, length(temp2),1) ,temp2, '.', 'MarkerEdgeColor', [0.3 0.3 0.3])
hold on
e=errorbar(2,mean(temp2), std(temp2)/sqrt(length(temp2)));
e.Marker='o';
e.MarkerFaceColor=[0.3 0.3 0.3];
e.MarkerEdgeColor=[0 0 0];
e.Color=[0 0 0];
hold on
plot([0.5 2.5], [0 0], 'k--')
xlim([0.5 2.5])
l=ylim;
xticks([1 2])
xticklabels({'with v.m.','control'})
%yticks(l(1):(l(2)-l(1))/10:l(2))
title(['mean dRt \pm ste of selected intervals, ROI', num2str(roi), ', p=', num2str(round((intervals_R_selected_stat(5)),3))])
ylabel('dRt')
text(1.1, mean(temp1), ['dRt=',num2str(round(mean(temp1),3)), '\pm', num2str(round(std(temp1)/sqrt(length(temp1)),3))])
text(2.1, mean(temp2), ['dRt=',num2str(round(mean(temp2),3)), '\pm', num2str(round(std(temp2)/sqrt(length(temp2)),3))])

clearvars temp temp1 temp2 e b l f1 f2


 
