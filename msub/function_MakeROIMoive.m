function ROIs_Movie=MakeROIMoive(Movie,pts,ROIsize)

NROI=size(pts,1);
ROIs=zeros(2,ROIsize,NROI);
ROIs_Movie=zeros(ROIsize,ROIsize,size(Movie,3),NROI);
for i=1:NROI
    % dimensions of ROIs_Movie: (X,Y,Frame,ROI#)
    % won't work for 1x1
    ROIs(1,:,i)=pts(i,1)-(ROIsize-1)/2 : pts(i,1)+(ROIsize-1)/2;
    ROIs(2,:,i)=pts(i,2)-(ROIsize-1)/2 : pts(i,2)+(ROIsize-1)/2;
    ROIs_Movie(:,:,:,i)=Movie(ROIs(2,:,i),ROIs(1,:,i),:);
end