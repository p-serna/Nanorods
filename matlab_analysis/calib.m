function Mask_CMOS=calib(Mask, ax, ay, bx, by)

Mask_CMOS=zeros(size(Mask,1), size(Mask,2));
for Xc=1:size(Mask_CMOS,2)
    for Yc=1:size(Mask_CMOS,1)
        Xe=round(ax*Xc+bx);
        Ye=round(ay*Yc+by);
        if Xe>0 && Xe<size(Mask,2)+1 && Ye>0 && Ye<size(Mask,1)
            Mask_CMOS(Yc,Xc)=Mask(Ye, Xe);
        end
    end
end
Mask_CMOS(isnan(Mask_CMOS))=0;
SE = strel('square',2);
Mask_CMOS = imdilate(Mask_CMOS,SE);
Mask_CMOS(Mask_CMOS==0)=NaN;