from numpy import *
from matplotlib.pylab import *
import os
from scipy import signal
from scipy.fftpack import fft
import scipy.optimize as opt
import time



wdir = "./"
fs = 1
fs = 200

if len(sys.argv)>4:
    try:
        fs = float(sys.argv[1])
        bt = float(sys.argv[2])
        bw = float(sys.argv[3])
        wdir = sys.argv[4]
        if wdir[-1] != '/':
            wdir = wdir+'/'
        print("Ok, we are using files in folder ",wdir)
        print("where the sampling freq is %.1f Hz and the target freq %.2f Hz with a bandwidth of %.2f Hz"%(fs,bt,bw))
    except:
        print("Something went wrong with the args: (sampling frequency, band target, band width directory)")
else:
    print("Not enough arguments, we do it for current folder with a sampling 1 Hz")

outputdir = wdir+'spect/'
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

putative = []
putative2 = []

scores = zeros((len(dfiles),3,2))
for i,filepath in enumerate(dfiles):
    namesp = (filepath.split(sep=".")[0]).split(sep="/")[-1]
    rois = load(filepath)
    sh = shape(rois)
    for j in range(2):
        xtt = rois*1.0
        xtt = transpose(xtt)-j*xtt.min(axis=1)
        xtt = xtt.mean(axis=0)

        S = xtt
        Sshape = S.shape
        nt = namesp
        namefigure = nt

        St = S*1.0

        f, t, Sxx = signal.spectrogram(St, fs,nperseg=100)
        sel = abs(f-bt)<bw
        sel0 = abs(f-bt)<min(bw/2.0,1.0)
        bands = Sxx[sel0,:];
        wbands = Sxx[sel,:];
       
        mab = max(wbands.flatten())
        meb = mean(wbands.flatten())
        stb = std(wbands.flatten())

        aSxmax = arange(Sxx.shape[1])[(bands.max(axis=0)-meb)/stb>5]
        if aSxmax.shape[0]>0:
            print(namesp,"more than 5 sigmas in spectrogram at times",t[aSxmax])
            if aSxmax.shape[0]>1:
                dd = aSxmax[1:]-aSxmax[:-1]
                de = concatenate(([1],1.0*(dd==1)))
                nintervals= sum((de[1:]-de[:-1])==1)
                if nintervals>0:
                    print("#-------------- With ",nintervals," long intervals!")
                    if j ==0:
                        putative.append(i)
                    if j ==1:
                        putative2.append(i)
        
        # Number of samplepoints
        N = St.shape[0]
        # sample spacing
        T = 1.0/fs
        x = np.linspace(0.0, N*T, N)
        y = St
        yf = fft(y)
        xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
        yyf =  2.0/N * np.abs(yf[8:N//2])
        # ~ plot(xf[8:],yyf,'.-',alpha=0.3)
        sel = (abs(xf[8:]-bt)<bw/2.0)*(abs(xf[8:]-bt)>0.1)
        nbw = array([mean(yyf[sel]),std(yyf[sel])])
        sel0 = (abs(xf[8:]-bt)<0.1)
        ybw = yyf[sel0]
        xhi = xf[8:][sel0][argmax(ybw)]
        bbw = ybw
        
        pbw = array([max(bbw),mean(bbw)])
        score = concatenate(((pbw-nbw[0])/nbw[1],[xhi]))
        
        if(score[0]>5): print(namesp,"more than 5 sigmas for freq=",xhi)
        # ~ if(score[0]>3): print(namesp,"more than 3 sigmas for freq=",xhi)

        scores[i,:,j] = score
        
save(outputdir+"scores.npy",scores)
save(outputdir+"putative.npy",array(putative+putative2))

sel = (scores[:,0,1]>3.0)+(scores[:,0,0]>3.0)

# ~ sel = argsort(-scores[:,0,1])
print(putative)
for i in arange(len(dfiles))[sel]:
    filepath  = dfiles[i]
    namesp = (filepath.split(sep=".")[0]).split(sep="/")[-1]
    rois = load(filepath)
    sh = shape(rois)
    xtt = rois*1.0
    xtt = transpose(xtt)-xtt.min(axis=1)
    xtt = xtt.mean(axis=0)

    S = xtt
    Sshape = S.shape
    nt = namesp
    namefigure = nt
    St = S*1.0
    
    f, t, Sxx = signal.spectrogram(St, fs,nperseg=100)
    sel = abs(f-bt)<bw
    sel0 = abs(f-bt)<min(bw/2.0,1.0)
    bands = Sxx[sel0,:];
    wbands = Sxx[sel,:];
    
    mab = max(wbands.flatten())
    meb = mean(wbands.flatten())
    stb = std(wbands.flatten())
    
    fig = figure("Spectrogram")
    imt = Sxx*1.0
    imt[Sxx>meb+2*stb] = meb+2*stb
    imt[Sxx<meb-1*stb] = meb-1*stb
    pcolormesh(t,f,imt,cmap="inferno")
    arrow(-2,bt,1,0,head_width=0.5)
    xlim(-2,60)
    ylabel("Frequency (Hz)")
    xlabel("Time")
    savefig(outputdir+"sp100ms_"+namesp+".png")
    title(namesp)
    # ~ fig.show()
    #i = input("...")
    # ~ fig.clear()
    close("Spectrogram")

    fig = figure("Fourier Transform")    
    N = St.shape[0]
    # sample spacing
    T = 1.0/fs
    x = np.linspace(0.0, N*T, N)
    y = St
    yf = fft(y)
    xf = np.linspace(0.0, 1.0/(2.0*T), N//2)
    yyf =  2.0/N * np.abs(yf[8:N//2])
    plot(xf[8:],yyf,'.-',alpha=0.5)
    sel = (abs(xf[8:]-bt)<bw/2.0)
    my = max(yyf[sel])
    arrow(bt,1.3*my,0,-my*0.1,head_width=1,head_length=my*0.1)
    xlabel("f(Hz)")
    ylabel("|FT(F)|")
    savefig(outputdir+"FT_"+namesp+".png")
    close("Fourier Transform")
    
