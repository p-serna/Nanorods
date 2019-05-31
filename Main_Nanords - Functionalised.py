# coding: utf-8
# Now some text explaining stuff

# ## General requirements

from matplotlib.pylab import *
from numpy import *
from Nanorods.sub.subs import *
from matlab.function_EMCCD_CMOS_calib import *
import PIL
import PIL.ExifTags as Exiftags
from matlab.tifmethods import *
from matlab.function_EMCCD_CMOS_calib import EMCCD_CMOS_calib
from matlab.calib import extractpeaks2d0,afintransf, CMOS_2fields_calib
from scipy.optimize import minimize


def stitchedtrace(x,framesxperiod = 4, index = False):
    size = x.shape[0]
    if size % framesxperiod != 0:
        print('Not divisible by %i!' % framesxperiod)
        return x
    xb = x.reshape(size//framesxperiod,framesxperiod)
    idx = arange(size).reshape(size//framesxperiod,framesxperiod)
    sel = isnan(xb.sum(axis=1))
    xb[sel,:] = nan
    # only those periods that are not nans go to new index
    idx = idx[~sel,:]
    xb = xb[~sel,:]
   
    if index:
        return(xb.flatten(),idx.flatten())
    return(xb.flatten())


# In[ ]:


def FFT_score(signal, Tacq, length_trh=100, lint=20):
    NROI=signal.shape[1]
    score_FFT=zeros(NROI)

    for roi in range(NROI):
        sig=stitchedtrace(signal[:,roi])
        Nframes=sig.shape[0]
        if Nframes>length_trh:
            Y = fft(sig)
            P2 = abs(Y/Nframes)
            P1 = P2[arange(0,Nframes//2+1)]
            #P1[1:-1] = 2*P1[1:-1]
            #print(Y.shape,P2.shape,P1.shape)
            f = (1/Tacq)*(arange(0,Nframes//2+1))/Nframes
            
            ind_start = argmin(abs(f-0.99/(4*Tacq)))
            ind_stop = argmin(abs(f-1.01/(4*Tacq)))
            inx_start = argmin(abs(f-(1/(4*Tacq)-lint/2)))
            inx_stop = argmin(abs(f-(1/(4*Tacq)+lint/2)))

            mP=mean(P1[inx_start:inx_stop])
            sP=std(P1[inx_start:inx_stop])
            temp1=zeros(ind_stop-ind_start)
            temp1=(P1[ind_start:ind_stop]-mP)/sP
            #print(mP, sP, temp1[0:2])
            score_FFT[roi]=max(temp1)
        else:
            score_FFT[roi]=nan
    return score_FFT

def movavgy (img, wind):
    sizeX, sizeY = img.shape
    imga=img*0
    wind2=wind//2
    for y in range(sizeY):
        for x in range(wind2,sizeX-wind2):
            imga[x,y]=nanmean(img[x-wind2:x+wind2+1,y])
    return imga
def movavgx (img, wind):
    return movavgy(img.transpose(), wind).transpose()

def calib(Mask, ax, ay, bx, by):
    Mask_CMOS=zeros(Mask.shape);
    for Xc in range(Mask_CMOS.shape[1]):
        for Yc in range(Mask_CMOS.shape[0]):
            Xe=int(ax*Xc+bx);
            Ye=int(ay*Yc+by);
            if Xe>-1 and Xe<Mask.shape[1] and Ye>-1 and Ye<Mask.shape[0]:
                Mask_CMOS[Yc,Xc]=Mask[Ye, Xe];
    Mask_CMOS[isnan(Mask_CMOS)]=0
    
    #kernel=2
    #kernelm = cv2.UMat(ones((kernel,kernel),float32)/float(kernel)**2)
    #dilated = cv2.dilate(Mask_CMOS,kernel=kernelm)
    return Mask_CMOS #array(dilated.get())

def MakeROIMovie(Movie,pts,ROIsize):

    NROI=shape(pts)[0]
    ROIs=zeros((2,ROIsize,NROI))
    ROIs_Movie=zeros((shape(Movie)[0],ROIsize,ROIsize,NROI))
    for i in range(NROI):
        # dimensions of ROIs_Movie: (X,Y,Frame,ROI#)
        # won't work for 1x1
        r2=(ROIsize-1)//2
        x0,xf = pts[i,0]-r2 , pts[i,0]+r2+1
        y0,yf = pts[i,1]-r2 , pts[i,1]+r2+1
        ROIs_Movie[:,:,:,i]=Movie[:, y0:yf , x0:xf]
    return ROIs_Movie


def extractLinearROIs(movie, pts_all, ROIsize):
    NROI=shape(pts_all)[0]
    Nframes, sizeY, sizeX = movie.shape
    ROIs_Movie = MakeROIMovie(movie,array(pts_all,dtype=int),ROIsize)
    # we reshape it to have ROIS as linear vectors
    Linear_ROIs_Movie=reshape(ROIs_Movie,(Nframes, ROIsize**2,NROI))
    
    # Sort pxs each frame
    Linear_ROIs_Movie.sort(axis=1)
    # We decide the background and the signal
    Linear_ROIs_Movie_Back=Linear_ROIs_Movie[:,0:ROIsize,:]
    Linear_ROIs_Movie_Signal=Linear_ROIs_Movie[:,ROIsize:,:]
    
    # Final signal and sum of the 2 channels
    Signal=sum(Linear_ROIs_Movie_Signal-mean(Linear_ROIs_Movie_Back, axis=1, keepdims=True), axis=1)/Tacq   
    Signal_2ch=Signal[:,0:NROI//2]+Signal[:,NROI//2:]
    
    return(Signal, Signal_2ch)

def classifier0(Signal_2ch, threshold = 0.45):
    temp1=mean(Signal_2ch, axis=0)
    temp2=mean((Signal_2ch-temp1)**2, axis=0)
    temp3=mean((Signal_2ch-temp1)**3, axis=0)/temp2**(3.0/2.0)
    temp4=mean((Signal_2ch-temp1)**4, axis=0)/temp2**(4.0/2.0)
    classifier=(temp3**2+1)/temp4

    idx = (-classifier).argsort()
    # selecting ROIs

    ROI_selected=arange(len(classifier))[classifier>threshold]

    if sum(ROI_selected)==0:
        print('no ROI selected')
    return(classifier, ROI_selected)

def getSignalsep(Signal,ROI_selected,NROI = Signal.shape[0]//2):
    Signalsep=Signal[:,ROI_selected]
    Signalsep=column_stack((Signalsep, Signal[:,NROI//2+ROI_selected])) 
    return(Signalsep)

def gaussd(x,par): return(exp(-(x-par[0])**2/2/par[1]**2)/sqrt(2*pi*par[1]**2) )
def dblgaussd(x,par): return(par[0]**2*gaussd(x,par[1:3])+(1.0-par[0]**2)*gaussd(x,par[3:]))

def dblgausfit(x,y,wy=1.0,par0=array([sqrt(0.5),-1,.5,1,.5])):
    def minf(par): return( sum(wy*(dblgaussd(x,par)-y)**2)/sum(wy))
    minx  = minimize(minf,par0)
    return minx


# In[ ]:


def get_treshold(Signal, NROI):
    threshold=zeros(NROI)
    for j in range(NROI):
        A, edges=histogram(Signal[:,j],40, density=True)
        xData, yData = edges[0:-1], A
        sig_mean=mean(Signal[:,j]);
        sig_std=std(Signal[:,j]);
        seed=array([1/2, sig_mean-sig_std, sig_std/2, sig_mean+sig_std, sig_std/2]);

        fitresult = dblgausfit(xData, yData, par0=seed);
        _, lm, ls, hm, hs=fitresult.x;
        if lm>hm:
            tm, ts = lm, ls
            lm, ls = hm, hs
            hm, hs = tm, ts
        if lm+ls*2<hm:
            threshold[j]=lm+ls*2
        else:
            threshold[j]=mean(lm, hm);
    return threshold



def main():
    windows = False

    if windows:
        folder = r"C:\Users\ludwig\data\NR sample files\18_11_29_pd23_11_div6_WIS_NR-BeRST\image files"
        fname = r"\wave\cell4_7.tif"
        fwave = r"\wave\cell4_7.tif"
        fctrl = r"\control\cell4_8.tif"
        file_Berst='\cell4_BeRST.tif'
        foldercalib = "\EMCCD-CMOS calib\\"
    else:
        folder = '/mnt/data/Anastasia/sample folder with NR data/image files/'
        fname = 'wave/cell1_1.tif'
        fwave = 'wave/cell1_1.tif'
        fctrl = 'control/cell1_2.tif'
        file_Berst='cell1_BeRST.tif'
        foldercalib = "EMCCD-CMOS calib/"

    #Only to visualize it:
    #fname = folder+file_Berst
    #movie = readBigTifFile(folder+fname)
    #movie_sum = movie.sum(axis=0) # sum(movie,axis=0)

    #visualization(movie)


    # ##     Mask generation and three fields alignment.
    # 
    # ### cMOS camera fields alignment

    f=3 #averaging factor for the mask
    ax, ay, bx, by = EMCCD_CMOS_calib(folder+foldercalib)

    print('Calibration parameters',ax, ay,bx, by)

    # ### EMCCD mask generation
    # ### Transposition of mask two CMOS images
    # 

    info = readtifInfo(folder+file_Berst)
    NframesStr = [d for d in info[270].split('\n') if d.find('frames')>=0][0].split('=')[-1]
    Nframes = int(NframesStr)
    Stim=tile([-60, 0],Nframes//2);

    movie = readBigTifFile(folder+file_Berst)
    movie_sum = movie.sum(axis=0)
    SizeX,SizeY = movie_sum.shape


    # Get Coefficients
    CorCoef_wave = zeros((SizeX, SizeY))
    for x in range(SizeX):
        for y in range(SizeY):
            temp = corrcoef(movie[:,x,y],Stim)
            CorCoef_wave[x,y]=temp[1,0];

    #imshow(CorCoef_wave, cmap='hot')
    #colorbar()

    thld=mean(CorCoef_wave)+2.2*std(CorCoef_wave)
    Mask_EMCCD=CorCoef_wave*1
    Mask_EMCCD[Mask_EMCCD<thld]=0;
    Mask_EMCCD[Mask_EMCCD>=thld]=1;
    #imshow(Mask_EMCCD);




    temp1=movavgx (Mask_EMCCD,5);
    temp2=movavgy (temp1,5);
    #imshow(temp1)
    #figure()
    #imshow(temp2)
    #figure()
    temp=temp2.flatten()
    #hist(temp[temp>0],51)
    temp3=Mask_EMCCD*1;
    temp3[temp2<=0.15]=nan;
    temp3[temp3==0]=nan;
    #figure()
    #imshow(temp3)



    f=3
    temp1=movavgx(temp3,f);
    Mask_EMCCD_exp=movavgy (temp1,f);
    #imshow(Mask_EMCCD_exp)


    #imshow(Mask_EMCCD_exp)
    #ylim(512,-10)




    Mask_CMOS=calib(Mask_EMCCD_exp, ax, ay, bx, by)
    #imshow(Mask_CMOS)


    # ##   Local maxima detection and ROI detection
    # ### Local maxima detection
    # ### Blinking coefficient and selection
    # ### Neural-net?
    # 



    # In[ ]:


    v = CMOS_2fields_calib(folder+foldercalib)


    # In[ ]:


    # ~ print(v)


    # In[ ]:


    calibim = readtifImage(folder+foldercalib+'calib_CMOS_2.tif')


    # In[ ]:


    ptsB = extractpeaks2d0(calibim[:,:512])


    # In[ ]:


    ptsR2 = 1.0*ptsB
    ptsR2[:,0] += 512+v[0]
    ptsR2[:,1] += v[1]
    # ~ visualization(calibim,row_stack((ptsB[:100,:],ptsR2[:100,:])))


    # In[ ]:


    dx, dy = 512+v[0],v[1]


    # In[ ]:


    #folder = r"C:\Users\ludwig\data\NR sample files\18_11_29_pd23_11_div6_WIS_NR-BeRST\image files"
    fname = fwave
    movie = readBigTifFile(folder+fname)
    movie_sum = movie.sum(axis=0) # sum(movie,axis=0)


    # In[ ]:


    info = readtifInfo(folder+fname)
    NframesStr = [d for d in info[270].split('\n') if d.find('frames')>=0][0].split('=')[-1]
    Nframes = int(NframesStr)


    # In[ ]:


    TacqStr = [d for d in info[270].split('\n') if d.find('finterval')>=0][0].split('=')[-1]
    Tacq = float(TacqStr)


    # In[ ]:


    # Generating stimulation trace
    Stim=tile([0, 0, -60, -60],Nframes//4);


    # In[ ]:


    # ~ visualization(movie_sum)


    # In[ ]:


    kernel = 5
    ROIsize=5
    km = ones((kernel,kernel),float32)#/float(kernel)**2
    J = zeros(Mask_CMOS.shape)
    rs2 = kernel//2
    mx,my = Mask_CMOS.shape

    for x in range(0,mx-kernel):
        for y in range(0,my-kernel):
            temp = km*movie_sum[x:(x+kernel),y:(y+kernel)]
            J[x+rs2,y+rs2] = temp.max()


    # In[ ]:


    mx,my = Mask_CMOS.shape
    pts = [] 
    for x in range(2*ROIsize+1,mx-2*ROIsize):
        for y in range(2*ROIsize,my-2*ROIsize):
            if Mask_CMOS[x,y]>0 and J[x, y]==movie_sum[x, y]:
                pts.append([y, x])


    pts = array(pts)


    # In[ ]:


    # ~ visualization(movie_sum,pts,figsize=(14,14))


    # In[ ]:


    x=pts[:,0]+dx
    y=pts[:,1]+dy


    # In[ ]:


    sel = x+floor(ROIsize/2)<movie_sum.shape[1]
    sel = sel*(y+floor(ROIsize/2)<movie_sum.shape[0])
    sel = sel*(y-floor(ROIsize/2)>0)
    x = x[sel]
    y = y[sel]


    # In[ ]:


    ptsr=column_stack((x,y))
    ptsb=pts*1
    pts_all=row_stack((ptsb,ptsr))
    # ~ visualization(movie_sum,pts_all,figsize=(14,14))


    # ##   Processing ROIs:
    # ### Linearization of ROI
    # ### Blinking substraction
    # ### Stitching cycles
    # ### Mobility
    # 
    # 



    NROI=shape(pts_all)[0]


    # In[ ]:





    # In[ ]:


    Signal, Signal_2ch = extractLinearROIs(movie,pts_all, ROIsize)
    # ~ temp = randint(Signal_2ch.shape[1])

    # ~ plot(Signal_2ch[:,temp])


    # In[ ]:


    # In[ ]:


    classifier, ROI_selected = classifier0(Signal_2ch, threshold = 0.45)


    # In[ ]:


    # ~ temp = randint(Signal_2ch.shape[1])
    # ~ #temp = 1
    # ~ plot(Signal_2ch[:,temp])
    # ~ print(temp,':',classifier[temp])
    # ~ print(ROI_selected)


    # In[ ]:




    # In[ ]:


    Signal_wave = getSignalsep(Signal,ROI_selected,NROI = NROI)

    pts_blue=pts_all[ROI_selected,:]
    pts_red=pts_all[ROI_selected+NROI//2,:]
    pts=row_stack((pts_blue, pts_red))


    # In[ ]:


    # Same for the control
    fname = fctrl
    movie = readBigTifFile(folder+fname)
    # Do we use same points: or we get new ones?
    Signalct, Signal_2ch_ctrl = extractLinearROIs(movie,pts_all, ROIsize)
    classifierct, ROI_selectedct = classifier0(Signal_2ch_ctrl, threshold = 0.45)


    # In[ ]:


    # ~ print(ROI_selectedct, ROI_selected)


    # In[ ]:


    ROI_sel=[]
    for i in ROI_selected:
        if i in ROI_selectedct:
            ROI_sel.append(i)
    # ~ print(ROI_sel)
    ROI_sel=array(ROI_sel)


    # In[ ]:


    Signal_wave = getSignalsep(Signal,ROI_sel,NROI = NROI)
    Signal_ctrl = getSignalsep(Signalct,ROI_sel,NROI = NROI)
    NROI=shape(Signal_wave)[1]//2
    Signal_2ch_wave = Signal_wave[:,0:NROI] + Signal_wave[:,NROI:]
    Signal_2ch_ctrl = Signal_ctrl[:,0:NROI] + Signal_ctrl[:,NROI:]


    # ~ temp=randint(ROI_sel.shape[0])
    # ~ figure(figsize=(18,6))
    # ~ plot(arange(6000)*0.01,Signal_wave[:,temp])
    # ~ plot(arange(6000)*0.01,Signal_wave[:,temp+ROI_sel.shape[0]], alpha=.7)

    # ~ plot(arange(6000)*0.01+60,Signal_ctrl[:,temp])
    # ~ plot(arange(6000)*0.01+60,Signal_ctrl[:,temp+ROI_sel.shape[0]], alpha=.7)


    # In[ ]:


    # In[ ]:



    threshold_2ch_wave=get_treshold(Signal_2ch_wave, NROI)
    threshold_2ch_ctrl=get_treshold(Signal_2ch_ctrl, NROI)


    # In[ ]:


    # ~ temp=randint(27)
    # ~ figure(figsize=(18,6))
    # ~ plot(arange(6000)*0.01,Signal_2ch_wave[:,temp])
    # ~ plot(arange(6000)*0.01, arange(6000)*0+threshold_2ch_wave[temp])


    # In[ ]:


    OnState_wave=Signal_2ch_wave<=threshold_2ch_wave
    # ~ print(OnState_wave.sum(axis=0))


    # In[ ]:


    # ~ print(sum(Signal_2ch_wave[:,1]<threshold_2ch_wave[1]))


    # In[ ]:


    OnState_wave=Signal_2ch_wave<=threshold_2ch_wave
    Signal_2ch_OnState_wave=Signal_2ch_wave*1
    Signal_2ch_OnState_wave[OnState_wave]=nan;

    OnState_wave_2=column_stack((OnState_wave,OnState_wave))

    Signal_OnState_wave=Signal_wave*1;
    Signal_OnState_wave[OnState_wave_2]=nan;


    # In[ ]:


    # ~ figure(figsize=(18,6))
    # ~ plot(arange(6000)*0.01,Signal_OnState_wave[:,temp])
    # ~ plot(arange(6000)*0.01,Signal_OnState_wave[:,temp+ROI_sel.shape[0]])


    # In[ ]:


    OnState_ctrl=Signal_2ch_ctrl<=threshold_2ch_ctrl
    Signal_2ch_OnState_ctrl=Signal_2ch_ctrl*1
    Signal_2ch_OnState_ctrl[OnState_ctrl]=nan;

    OnState_ctrl_2=column_stack((OnState_wave,OnState_ctrl))

    Signal_OnState_ctrl=Signal_ctrl*1;
    Signal_OnState_ctrl[OnState_ctrl_2]=nan;


    # In[ ]:


    Ratio_wave=zeros((Nframes,NROI))
    Ratio_wave=Signal_OnState_wave[:,0:NROI]/Signal_2ch_OnState_wave
    Ratio_ctrl=zeros((Nframes,NROI))
    Ratio_ctrl=Signal_OnState_ctrl[:,0:NROI]/Signal_2ch_OnState_ctrl


    # In[ ]:


    split_wave=nanmean(Ratio_wave, axis=0)
    Sig_FFT_wave=zeros((Nframes,NROI))
    ind=arange(split_wave.shape[0])[split_wave>0.5]
    Sig_FFT_wave[:,ind]=Signal_OnState_wave[:,ind]
    ind=arange(split_wave.shape[0])[split_wave<=0.5]
    Sig_FFT_wave[:,ind]=Signal_OnState_wave[:,NROI+ind]

    split_ctrl=nanmean(Ratio_ctrl, axis=0)
    Sig_FFT_ctrl=zeros((Nframes,NROI));
    ind=arange(split_ctrl.shape[0])[split_ctrl>0.5]
    Sig_FFT_ctrl[:,ind]=Signal_OnState_ctrl[:,ind]
    ind=arange(split_ctrl.shape[0])[split_ctrl<=0.5]
    Sig_FFT_ctrl[:,ind]=Signal_OnState_ctrl[:,NROI+ind]




    score_FFT=FFT_score(Sig_FFT_wave,Tacq)

    return(score_FFT)


if __name__ == '__main__':
    main()
    
