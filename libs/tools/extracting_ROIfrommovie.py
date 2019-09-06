from numpy import zeros,save,row_stack,ones,array,column_stack,log
from numpy.random import permutation
from sub.subs import readBigtifFile,selROI
import cythLocalization as loc


def extractROIS(cfile,pos,ROIsize=5,output_folder=None,output_names=True):
    # ----------------------------------------------------------------
    # We read the file.     

    movie = readBigtifFile(cfile)
    tduration, mheight,mwidth = movie.shape
   
    rois = zeros((tduration,ROIsize*ROIsize,len(pos)))
    
    if type(ROIsize) is int:
        if ROIsize%2 == 0:
            print("Not odd size for ROIsize")
            return
            
        for i,cl in enumerate(pos):
            roi = selROI(movie,cl[0],cl[1],t=0,size=ROIsize//2,xlim=(0,mwidth),ylim=(0,mheight))
            roi = roi.reshape(tduration,ROIsize*ROIsize)
            rois[:,:,i] = roi
            if output_folder is not None:
                try:
                    save(output_folder+"roi"+str(i).zfill(4)+".npy",roi)
                except:
                    print("output_folder is not an address?")
    else:
        print("ROIsize is not integer")
        return   
         
    return(rois)


def fitRois(rois,nsel0 = 600, nsel1 = 100,sx = 2, sy = 2):
    '''rois: (tduration, roisize linearised( (5,5)=25), nrois)
    nsel0, nsel1 : for initial selection to make fits faster
    sx,sy : limits for size of sigma of distribution (in pixels)
    '''
    sh = rois.shape
    fullp = zeros((sh[2],sh[0],8))
    for i in range(sh[2]):
        testm = rois[:,:,i].transpose()
        sel = (-testm.max(axis=0)).argsort()[:nsel0]
        Amax = testm[:,sel[0]].max()
        sel2 = permutation(sel[nsel1:])[:nsel1]

        testm2 = testm/Amax

        arga = testm2.argmax(axis=0)
        argi = testm2.argmin(axis=0)

        testm2 = row_stack((testm2,arga,argi)).transpose()
        poptssel2 = array(loc.fitvideoPfxth(testm2[sel2,:])) 
        sel = (log(poptssel2[:,0])<15)*(poptssel2[:,1]>0)*(poptssel2[:,1]<sx)*(poptssel2[:,2]>0)*(poptssel2[:,2]<sy)

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
        
        fullp[i,:,:] = popts
    return(fullp)

def projectmovie(cfile):
    # ----------------------------------------------------------------
    # We read the file.     

    movie = readBigtifFile(cfile)
    tduration, mheight,mwidth = movie.shape
   
    return(movie.mean(axis=0))
