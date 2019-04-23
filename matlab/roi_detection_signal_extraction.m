%% --- user defined variables --- %
root='C:\Users\ludwig\Documents\MATLAB';
folders={'\19_02_05_pd1_02_div4_NR_BeRST'};
folders_wave='\image files\wave';
folders_contrl='\image files\control';
count=1;
for i=1:length(folders)
    files_wave{i} = dir(fullfile(root, folders{i}, folders_wave, '*.tif'));
    files_contrl{i} = dir(fullfile(root, folders{i}, folders_contrl, '*.tif'));
end

ROIsize=5;
calib_file='calib.tif';
clearvars i
%% extract some info from the files; it is assumed that tacq and nframes are identical for all recordings
for f=1:length(folders)
    imcalib=imread(fullfile(root, folders{f},calib_file));
    pts_calib = ReadROIs(fullfile(root, folders{f},[calib_file(1:end-3) 'csv']));
    pts_calib = CenterROIs(imcalib,pts_calib,9);
    dx(f)=round(mean(pts_calib(2:2:end,1)-pts_calib(1:2:end,1)));
    dy(f)=round(mean(pts_calib(2:2:end,2)-pts_calib(1:2:end,2)));
end

temp=files_wave{1};
info = imfinfo(fullfile(root, folders{1}, folders_wave, temp(1).name));
info = info(1, :);
NFramesStr = regexp(info.ImageDescription, 'images=(\d*)', 'tokens');
Nframes = str2double(NFramesStr{1});
TacqStr = regexp(info.ImageDescription, 'finterval=(\d?\D?\d*)','tokens');
Tacq= str2double(TacqStr{1});



% Generating stimulation trace
Stim=repmat([0 0 -60 -60],1,Nframes/4);

clearvars imcalib pts_calib TacqStr NFramesStr info f


%% read movie files, determine ROIs, select blinking traces, threshold
for f=1:length(folders)
    f_wave=files_wave{f};
    f_contrl=files_contrl{f};
    for ff=1:length(f_wave)
        %% reading wave movie
        folder_wave=fullfile(root, folders{f}, folders_wave);
        file_wave=f_wave(ff).name;
        
        Movie=ReadLongMovie(fullfile(folder_wave,file_wave), Nframes);
        SumMovie_wave=sum(Movie,3);     

        %% Generating ROIs within CMOS mask
        load(fullfile(root, folders{f}, 'mat files', [file_wave(1:5), '_mask.mat']), 'Mask_CMOS')
        
        clearvars pts
        ROIsize=5;
        n=1;
        SE = strel('square',5);
        J = imdilate(SumMovie_wave,SE);
        Mask_CMOS(Mask_CMOS==0)=nan;
        for x=2*ROIsize+1:size(Mask_CMOS,1)-2*ROIsize
            for y=2*ROIsize:size(Mask_CMOS,2)-2*ROIsize
                if ~isnan(Mask_CMOS(x,y))&& J(x, y)==SumMovie_wave(x, y)
                    pts(n,:)=[y x];
                    n=n+1;
                end
            end
        end
        
        
        %%  check if ROI are inside the Movie
        x=pts(:,1)+dx(f);
        y=pts(:,2)+dy(f);
        
        a=x+floor(ROIsize/2)<size(SumMovie_wave,2);
        x=x.*a;
        y=y.*a;
        
        a=y+floor(ROIsize/2)<size(SumMovie_wave,1);
        x=x.*a;
        y=y.*a;
        
        a=y-floor(ROIsize/2)>0;
        x=x.*a;
        y=y.*a;
        
        x(x==0)=[];
        y(y==0)=[];
        pts=[x-dx(f) y-dy(f)];
        
        
        clearvars idx x y n s J SE a
        
        %% define ROIs
        pts_blue=pts(:,1:2);
        pts_blue_center=CenterROIs(SumMovie_wave,pts_blue,ROIsize);
        
        % transtlate blue ROIs to red side of the camera and re-center
        pts_red = pts_blue_center+repmat([dx(f) dy(f)],size(pts_blue_center,1),1);
        pts_red_center = pts_red;%CenterROIs(SumMovie_wave,pts_red,ROIsize);
        
        ROIs_Movie = MakeROIMoive(Movie,[pts_blue_center;pts_red_center],ROIsize);
        
        NROI_auto=size(ROIs_Movie,4)/2;
        
        clearvars pts_blue pts_red;
        
        %% Determining ROI signal
        % Creation of a matrix with all ROI data. Dimensions of Linear_ROIs_Moive: (ROI Pixels,Frame,ROI#)
        Linear_ROIs_Movie=reshape(ROIs_Movie,ROIsize^2,1,Nframes,NROI_auto*2);
        Linear_ROIs_Movie=permute(Linear_ROIs_Movie,[1,3,4,2]);
        
        % Sorting ROI Pixels in Linear_ROIs_Moive from min to max
        for j=1:NROI_auto*2 % loop over all ROI
            for i=1:Nframes % loop over all frames
                Linear_ROIs_Movie(:,i,j)=sortrows(Linear_ROIs_Movie(:,i,j));
            end
        end
        
        Linear_ROIs_Movie_Back=Linear_ROIs_Movie(1:ROIsize,:,:);
        Linear_ROIs_Movie_Signal=Linear_ROIs_Movie(ROIsize+1:end,:,:);
        
        
        Signal_auto=squeeze(sum(Linear_ROIs_Movie_Signal)-mean(Linear_ROIs_Movie_Back)*(ROIsize^2-ROIsize))/Tacq;
        
        Signal_auto_2ch=Signal_auto(:,1:NROI_auto)+Signal_auto(:,NROI_auto+1:end);
        clearvars  Linear_ROIs_Movie_Back Linear_ROIs_Movie_Signal j i
        
        %% classifier
        
        temp1=mean(Signal_auto_2ch);
        temp2=mean((Signal_auto_2ch-temp1).^2);
        temp3=mean((Signal_auto_2ch-temp1).^3)./temp2.^(3.0/2.0);
        temp4=mean((Signal_auto_2ch-temp1).^4)./temp2.^(4.0/2.0);
        classifier=(temp3.^2+1)./temp4;
        [~, idx] = sort(classifier, 'descend');
        
        clearvars temp1 temp2 temp3 temp4
        % selecting ROIs
        
        clearvars ROI_selected
        ROI_selected=find(classifier>0.45);
        
        if sum(ROI_selected)==0
            disp('no ROI selected')
            return;
        end
        
        NROI=size(ROI_selected,2);
        Signal_wave=Signal_auto(:,[ROI_selected NROI_auto+ROI_selected]);
        
     
       
        pts_blue=pts_blue_center(ROI_selected,:);
        pts_red=pts_red_center(ROI_selected,:);
        
        %%  reading control file
          
        folder_contrl=fullfile(root, folders{f}, folders_contrl);
        file_contrl=f_contrl(ff).name;
        Movie=ReadLongMovie(fullfile(folder_contrl,file_contrl), Nframes);
        
        SumMovie_contrl=sum(Movie,3);
        pts_blue_center=CenterROIs(SumMovie_contrl,pts_blue,ROIsize);
        pts_red_center=CenterROIs(SumMovie_contrl,pts_red,ROIsize);
        
        ROIs_Movie = MakeROIMoive(Movie,[pts_blue_center;pts_red_center],ROIsize);
        
        %% Determining ROI signal in control movie
        % Creation of a matrix with all ROI data. Dimensions of Linear_ROIs_Moive: (ROI Pixels,Frame,ROI#)
        Linear_ROIs_Movie=reshape(ROIs_Movie,ROIsize^2,1,Nframes,NROI*2);
        Linear_ROIs_Movie=permute(Linear_ROIs_Movie,[1,3,4,2]);
        
        % Sorting ROI Pixels in Linear_ROIs_Moive from min to max
        for j=1:NROI*2 % loop over all ROI
            for i=1:Nframes % loop over all frames
                Linear_ROIs_Movie(:,i,j)=sortrows(Linear_ROIs_Movie(:,i,j));
            end
        end
        
        Linear_ROIs_Movie_Back=Linear_ROIs_Movie(1:ROIsize,:,:);
        Linear_ROIs_Movie_Signal=Linear_ROIs_Movie(ROIsize+1:end,:,:);
        
        
        Signal_contrl=squeeze(sum(Linear_ROIs_Movie_Signal)-mean(Linear_ROIs_Movie_Back)*(ROIsize^2-ROIsize))/Tacq;
        Signal_contrl_2ch=Signal_contrl(:,1:NROI)+Signal_contrl(:,NROI+1:end);
        clearvars  Linear_ROIs_Movie_Back Linear_ROIs_Movie_Signal j i
        
         %% classifier control
        
        temp1=mean(Signal_contrl_2ch);
        temp2=mean((Signal_contrl_2ch-temp1).^2);
        temp3=mean((Signal_contrl_2ch-temp1).^3)./temp2.^(3.0/2.0);
        temp4=mean((Signal_contrl_2ch-temp1).^4)./temp2.^(4.0/2.0);
        classifier=(temp3.^2+1)./temp4;
        [~, idx] = sort(classifier, 'descend');
        
        clearvars temp1 temp2 temp3 temp4
        % selecting ROIs
        
        clearvars ROI_selected_contrl
        ROI_contrl=find(classifier>0.45);
        
        if sum(ROI_contrl)==0
            disp('no ROI selected')
            return;
        end
      
        clearvars j button x y m n limy limx temp fitobject gof STD
        
       
%% 
       Exp_Data(count).folder=folder_wave;
       Exp_Data(count).file_wave=file_wave;
       Exp_Data(count).file_contrl=file_contrl;
       Exp_Data(count).Tacq=Tacq;
       Exp_Data(count).Nframes=Nframes;
       Exp_Data(count).Mask_CMOS=Mask_CMOS;
       Exp_Data(count).SumMovie_wave=SumMovie_wave;
       Exp_Data(count).SumMovie_contrl=SumMovie_contrl;
       Exp_Data(count).pts_blue=pts_blue;
       Exp_Data(count).pts_red=pts_red;
       Exp_Data(count).NROI=NROI;
       Exp_Data(count).ROI_contrl=ROI_contrl;
       Exp_Data(count).Signal_wave=Signal_wave;
       Exp_Data(count).Signal_contrl=Signal_contrl;
       
       count=count+1;
       
       clearvars -except root dx dy f ff f_contrl f_wave FileData files_contrl files_wave firstrun folders folders_contrl folders_wave Nframes ...
           Nframes Stim root ROIsize Tacq Exp_Data count
   
    end
end

save(fullfile(root, 'summary.mat'))