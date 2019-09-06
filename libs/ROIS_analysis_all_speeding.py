from numpy import *
from matplotlib.pylab import *
import os
from scipy import signal
#from core_libs import *
from scipy.fftpack import fft
import scipy.optimize as opt
import cythLocalization as loc
import time

#%matplotlib inline

def fitROI(roi,n0 = 100,name=""):
    ''' roi: video of x*y: shape=( # frames, # pixels (x*y))
        n0: to estimate sx and sy from n0 frames selected from 600 brightest frames
        name: just  a silly thing.
        It returns an array (nframes ,8) with columns:
        Amplitude, sx, sy, theta, background, x, y, error flag
        
        error flag is 1 if there was an error in that fit, 0 otherwise

    '''
    # shape: nframes, pixels (x*y)
    sh = roi.shape
    testm = roi.transpose()
    
    sel = (-testm.max(axis=0)).argsort()[:600]
    Amax = testm[:,sel[0]].max()
    sel2 = permutation(sel[n0:])[:n0]
    
    testm2 = testm/Amax
    
    arga = argmax(testm2,axis=0)
    argi = argmin(testm2,axis=0)
    
    testm2 = transpose(row_stack((testm2,arga,argi)))
    #poptssel = array(loc.fitvideoP(testm2[sel2,:])) 
    poptssel2 = array(loc.fitvideoPfxth(testm2[sel2,:])) 
    sel = (poptssel2[:,1]>0)*(poptssel2[:,1]<2)*(poptssel2[:,2]>0)*(poptssel2[:,2]<2)
    if sum(sel) <= 2:
        print("Screwed with this ROI "+name)
        print("No correct sigmas for any of selected pre fits")
        print("Attempting full fit")
        popts = array(loc.fitvideoP(testm2))
    else:
        poptssel2 = poptssel2[sel,:]
        s1, s2 = poptssel2.mean(axis=0)[1:3]
        testm3 = row_stack((testm2.transpose(),s1**2*ones(sh[0]),s2**2*ones(sh[0]))).transpose()
        popts = array(loc.fitvideoPfxsh(testm3)) 

        popts[:,:2] = popts[:,:2]*Amax

        # Rearranging columns to have same format as before:
        # Amplitude, sx, sy, theta,offset,x,y , error
        popts = column_stack((popts[:,:1],popts[:,4:6],0*popts[:,:1],popts[:,1:2],popts[:,2:4],popts[:,-1]))

    # This is the full fit - much slower!! 
    # ~ popts = array(loc.fitvideoP(testm2)) 
    return(popts)
       
def main():

    if len(sys.argv)>1:
        wdir = sys.argv[1]
        if wdir[-1] != '/':
            wdir = wdir+'/'
    else:
        wdir = "./"
        print("No argument, we do it for current folder")

    outputdir = wdir+'sptrack/'
    if not os.path.isdir(outputdir):
        try:
            os.system("mkdir "+outputdir)
        except ValueError:
            print("I cannot create the folder for the output!")
            raise SystemExit(0)

    dirt = wdir

    basedir = dirt
    files = os.listdir(basedir)

    dfiles = []
    for f in files:
        if f[:5]=='roi_s' and f[-3:] == 'npy': dfiles.append(basedir+f)
     
    dfiles.sort()

    for filepath in dfiles:
        print("Starting with "+filepath)
        namesp = (filepath.split(sep=".")[0]).split(sep="/")[-1]
        rois = load(filepath)
        sh = shape(rois)
        if len(sh)>2 :
            print("This is a different one: "+filepath)
        rois = reshape(rois,(1,sh[0],sh[1]))


        for i in range(rois.shape[0]):
            starttime = time.time()
            testm = transpose(rois[i,:,:])

          
            
            sel = (-testm.max(axis=0)).argsort()[:600]
            Amax = testm[:,sel[0]].max()
            sel2 = permutation(sel[100:])[:100]
            
            testm2 = testm/Amax
            
            arga = argmax(testm2,axis=0)
            argi = argmin(testm2,axis=0)
            
            testm2 = transpose(row_stack((testm2,arga,argi)))
            #poptssel = array(loc.fitvideoP(testm2[sel2,:])) 
            poptssel2 = array(loc.fitvideoPfxth(testm2[sel2,:])) 
            sel = (poptssel2[:,1]>0)*(poptssel2[:,1]<2)*(poptssel2[:,2]>0)*(poptssel2[:,2]<2)
            if sum(sel) <= 2:
                print("Screwed with file "+namesp)
                print("No correct sigmas for any of selected pre fits")
                print("Attempting full fit")
                popts = array(loc.fitvideoP(testm2))
            else:
                poptssel2 = poptssel2[sel,:]
                s1, s2 = poptssel2.mean(axis=0)[1:3]
                testm3 = row_stack((testm2.transpose(),s1**2*ones(sh[0]),s2**2*ones(sh[0]))).transpose()
                popts = array(loc.fitvideoPfxsh(testm3)) 
        
                popts[:,:2] = popts[:,:2]*Amax
        
                # Rearranging columns to have same format as before:
                # Amplitude, sx, sy, theta,offset,x,y , error
                popts = column_stack((popts[:,:1],popts[:,4:6],0*popts[:,:1],popts[:,1:2],popts[:,2:4],popts[:,-1]))

            # This is the full fit - much slower!! 
            # ~ popts = array(loc.fitvideoP(testm2)) 
             
            endtime = time.time()
            print("It took ",endtime-starttime,"seconds")
            save(outputdir+"posh_"+namesp+".npy",popts)
            
            #poptroi[i] = popts
        print("File "+filepath+" finished")


if __name__ == '__main__':
    main()
