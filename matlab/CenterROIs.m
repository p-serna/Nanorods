function pts_center=CenterROIs(SumMovie,pts,ROIsize)

NROI=size(pts,1);
ROIs=zeros(2,ROIsize,NROI);
SumROIMovie=zeros(ROIsize,ROIsize,NROI);
for i=1:NROI
    % dimensions of ROIs_Movie: (X,Y,Frame,ROI#)
    % won't work for 1x1
    ROIs(1,:,i)=pts(i,1)-(ROIsize-1)/2 : pts(i,1)+(ROIsize-1)/2;
    ROIs(2,:,i)=pts(i,2)-(ROIsize-1)/2 : pts(i,2)+(ROIsize-1)/2;
    SumROIMovie(:,:,i)=SumMovie(ROIs(2,:,i),ROIs(1,:,i));
end

center=ROIsize-floor(ROIsize/2);
pts_center=zeros(size(pts));
for i=1:NROI
    [Ix, Iy]=find(SumROIMovie(:,:,i)==max(max(SumROIMovie(:,:,i))),1);
    pts_center(i,1)=pts(i,1)+(Iy-center);
    pts_center(i,2)=pts(i,2)+(Ix-center);
end