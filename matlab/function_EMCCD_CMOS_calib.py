import os
import warning

def EMCCD_CMOS_calib(folder = './'):

    #Extracting images and csv files
    files_tif = []
    files_csv = []
    filenames = []
    for file in os.listdir(folder):
        if file.endswith('.tif'):
            files_tif.append(file)
            filenames.append(file.split('_')[0])
        elif file.endswith('.csv'):
            files_csv.append(file)
    
    filenames = list(set(filenames))
    filenames.sort()
    files_tif.sort()
    files_csv.sort()
    files = {name: [] for name in filenames}
        
    # Let us subdivide files with name of the cell
    for file in files_tif:
        filename = file.split('_')[0]
        files[filename].append(file)
    
    
    for i,name in enumerate(filenames)): #real traces
        fileCMOS = [file for file in files[name] if file.find('CMOS')>0]
        if len(fileCMOS) != 1:
            raise 'Sorry, there are '+str(len(fileCMOS))+' files with CMOS in their name, for '+name+'.'
        fileCMOS = fileCMOS[0]
        
        fileEMCCD = [file for file in files[name] if file.find('EMCCD')>0]
        if  len(fileEMCCD) != 1:
            if len(fileEMCCD)==0:
                raise 'We stop, no files with EMCCD in their name, for '+name+'.'
            else:
                warning.warn('There are '+str(len(fileEMCCD))+' files with EMCCD in their title, we take '+fileEMCCD[0])
        fileEMCCD = fileEMCCD[0]
        
        # Change these lines to appropiate reading of tif files:
        CMOS = imread(folder+fileCMOS)
        EMCCD= imread(folder+fileEMCCD)
        
        CMOS = CMOS[:,:CMOS.shape[1]//2]


#function [ax ay bx by]=EMCCD_CMOS_calib(folder)

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
