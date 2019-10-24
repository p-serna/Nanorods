function PlotROIs(SumMovie,pts,ROIsize,RB,norm, start, f)
% RB is a flag to number red side and blue side with the same number
% norm is a flag to use normalized image 
% [pts_blue ; pts_red]
% SumMovie=sum(Movie,3);
% f - figure handle

if nargin > 5
  s = start;
else
  s = 1;
end 

if nargin > 6
  fh = f;
else
  fh = figure;
end 

if norm img=SumMovie./max(max(SumMovie));
    img2=imadjust(img);%[0.0 0.1],[]
else img2=SumMovie;
end

imshow(img2)
NROI=size(pts,1);
j=1;
for i=1:NROI
    rectangle('Position',[pts(i,1)-(ROIsize/2) pts(i,2)-(ROIsize/2) ROIsize ROIsize],'EdgeColor','red','LineWidth',1)
    if RB
        if j==NROI/2+1
            j=1;
        end
    end
    text(pts(i,1)+(ROIsize/2),pts(i,2)+(ROIsize/2),num2str(s-1+j),'Color','red','FontSize',10,'FontWeight','bold');
    j=j+1;
end
