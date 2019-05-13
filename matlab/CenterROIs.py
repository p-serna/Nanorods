def centerROIs(img,pts,ROIsize = 5):
    rs2 = ROIsize//2
    ptsn = pts*1
    for i,p in enumerate(pts):
        x,y = p
        roi = img[(y-rs2):(y+rs2+1),(x-rs2):(x+rs2+1)]
        idx = roi.argmax()
        pxy, pxx = idx//5, idx%5
        ptsn[i,:] = [x+pxx-2,y+pxy-2]
    return(ptsn)
