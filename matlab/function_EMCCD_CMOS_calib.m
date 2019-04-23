function [ax ay bx by]=EMCCD_CMOS_calib(folder)
%folder='C:\Users\ludwig\Documents\MATLAB\18_11_29_pd23_11_div6_WIS_NR-BeRST\EMCCD-CMOS calib'; 
files = dir(fullfile(folder, '*.tif'));
files_csv = dir(fullfile(folder, '*.csv'));
%importing calib images
k=0;
n=1;
m=1;
for i=1:size(files) %real traces
    file = files(i).name;
    k=strfind(file, 'CMOS');
    if k~=0
        CMOS(:,:,n)=imread(fullfile(folder,file));
        n=n+1;
    else
        EMCCD(:,:,m)=imread(fullfile(folder,file));
        m=m+1;
    end 
end

CMOS=CMOS(:,1:size(CMOS,2)/2,:);
k=0;
n=1;
m=1;
for i=1:size(files_csv) %real traces
    file = files_csv(i).name;
    k=strfind(file, 'CMOS');
    if k~=0
        CMOS_pts{n}=ReadROIs(fullfile(folder,file));
        n=n+1;
    else
        EMCCD_pts{m}=ReadROIs(fullfile(folder,file));
        m=m+1;
    end 
end
 clearvars i file folder files k temp myvars k n m files_csv

% %Plotting images
% temp=squeeze(CMOS(1,:,:));
% img=imadjust(temp);
% figure;
% imshow(img)
% 
% clearvars img temp


% define ROIs
ROIsize=9;

for i=1:size(CMOS_pts,2)
    CMOS_pts_center{i}=CenterROIs(CMOS(:,:,i),CMOS_pts{i},ROIsize);
    EMCCD_pts_center{i}=CenterROIs(EMCCD(:,:,i),EMCCD_pts{i},ROIsize);
end
clearvars i CMOS_pts EMCCD_pts

%% Plot SumMoive with ROIs
% image=1;
% RB=0;
% norm=0;
% PlotROIs(CMOS(:,:,image),[CMOS_pts_center{image}],ROIsize,RB, norm)
% PlotROIs(EMCCD(:,:,image),[EMCCD_pts_center{image}],ROIsize,RB, norm)
% clearvars RB norm image

%%
ptsC=vertcat(CMOS_pts_center{:});
ptsE=vertcat(EMCCD_pts_center{:});
Xc=ptsC(:,1);
Yc=ptsC(:,2);
Xe=ptsE(:,1);
Ye=ptsE(:,2);

fitobject = fit(Xc,Xe,'poly1');
ax=fitobject.p1;
bx=fitobject.p2;

fitobject = fit(Yc,Ye,'poly1');
ay=fitobject.p1;
by=fitobject.p2;