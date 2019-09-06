import sys
sys.path.append("/users/bssn/serna/GitIBENS/Nanorods/")
from ROIS_analysis_all_speeding import fitROI

%pylab
def ftvaugmentr(img0,magn=2):
    img = img0.transpose()
    sh = img.shape
    sh2 = array([sh[0],sh[1]])*(magn-1)
    sh2 = ((sh2[0])//2+1,(sh2[0])//2+1)
    fftim  = fftshift(fft2(img,axes=(0,1)),axes=(0,1))
    fftim  =  pad(fftim,(sh2[0],sh2[1]),'constant') 
    imgn = real(ifft2(ifftshift(fftim,axes=(0,1)),axes=(0,1)))
    imgn = imgn[:,:,sh2[0]:-sh2[0]]
    return(imgn.transpose())
    
rois = load("roi_moviewave.npy")     
roirs = rois.reshape(rois.shape[0],6000,5,5)  
#roia = map(ftvaugment,roirs[:4,:,:,:])


for i,roi in enumerate(roirs):
    roia = ftvaugmentr(roi,magn=2)
    roia = roia.reshape(6000,121)
    popts = fitROI(roia)
    save("sptrack2/posh_"+str(i).zfill(4)+".npy",popts)
    print(" ROI: ",i)

for i,roi in enumerate(roirs):
    roia = roi.reshape(6000,25)
    popts = fitROI(roia)
    save("sptrack/posh_"+str(i).zfill(4)+".npy",popts)
    print(" 2-ROI: ",i)
