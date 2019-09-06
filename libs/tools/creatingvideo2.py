import sys, os
sys.path.append("/users/bssn/serna/GitIBENS/Nanorods/")
from tools.creatingvideo_sub import ROIanimation
from matplotlib
wdirt = "selpeaks/"
files = os.listdir(wdirt)
dfiles = []
for f in files:
    if f[:11] == "rois_random":
        dfiles.append(wdirt+f)


dfiles.sort()

k = 0
for fname in dfiles:
    roi = load(fname)
    sh22 = roi.shape[2]//2
    for i in range(sh22):
        roip = roi[:,:,i]+roi[:,:,i+sh22]
        St,ani,fig  = ROIanimation(roip.reshape(6000,5,5))
        k += 1
        ani.save("movies/movie_roi"+str(k).zfill(4)+".mp4",fps =  100) 
        close()
