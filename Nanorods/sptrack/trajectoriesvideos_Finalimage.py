from numpy import *
from matplotlib.pylab import *
import h5py,os
from scipy import signal
#from core_libs import *
from scipy.fftpack import fft
import scipy.optimize as opt
from scipy.optimize import minimize
import pickle
from matplotlib.collections import LineCollection
from scipy.fftpack import fft2,fftshift,ifftshift,ifft2

def gaussd(x,par): return(exp(-(x-par[0])**2/2/par[1]**2)/sqrt(2*pi*par[1]**2) )
def dblgaussd(x,par): return(par[0]**2*gaussd(x,par[1:3])+(1.0-par[0]**2)*gaussd(x,par[3:]))

def dblgausfit(x,y,wy=1.0,par0=array([sqrt(0.5),-1,.5,1,.5])):
    def minf(par): return( sum(wy*(dblgaussd(x,par)-y)**2)/sum(wy))
    minx  = minimize(minf,par0)
    return minx
    
    
def running_mean(x, N):
    cumsumt = cumsum(insert(x, 0, 0)) 
    return (cumsumt[N:] - cumsumt[:-N]) / float(N)

if len(sys.argv)>1:
    wdir = sys.argv[1]
    if wdir[-1] != '/':
        wdir = wdir+'/'
    if len(sys.argv)>2:
        driftcorrected = sys.argv[2]=='1'
    else:
        driftcorrected = False
else:
    wdir = "./"
    driftcorrected = False
    print("No argument, we do it for current folder")

if driftcorrected:
    outputdir = wdir+"moviesDC/"
    print("We are correcting the drift for dir ",wdir)
else:
    outputdir = wdir+"movies0/"
    print("No drift is corrected for dir ",wdir)
    
if not os.path.isdir(outputdir):
    try:
        os.system("mkdir "+outputdir)
    except ValueError:
        print("I cannot create the folder for the output in:",outputdir)
        raise SystemExit(0)
        
        
        
dirt = wdir
basedir = dirt
files = os.listdir(basedir)



dfiles = []

if driftcorrected:
    for f in files:
        if f[-4:]=='.npy'and f[:4]=='posh' and f[-6:-4]=='DC': dfiles.append(basedir+f)
else:
    for f in files:
        if f[-4:]=='.npy'and f[:4]=='posh' and f[-6:-4]!='DC': dfiles.append(basedir+f)
    
dfiles.sort()

namesm = dfiles

cmapgnu = get_cmap("gnuplot")

#ths = load("thresholds.npy")
data = []

selROI = arange(len(dfiles))
try:
    ths = load(wdir+"thresholds.npy")
    if len(ths)<len(namesm):
        ths = concatenate((ths,zeros(len(namesm)-len(ths))))
except:
    ths = zeros(len(namesm))

#fig = figure()
#ion()

nokarray = zeros(len(namesm))
# ~ for i,name in enumerate(namesm):
    # ~ ax = fig.add_subplot(111)
    # ~ nameNR = name.split("_")[-2]
    # ~ popts = load(name)
    # ~ posx = popts[:,5]
    # ~ posy = popts[:,6]
    # ~ amp = popts[:,0]
    # ~ th = ths[i]
    # ~ th2 = 13.0
    # ~ sel = (popts[:,-1]==0)*(posx>0.1)*(posy>0.1)*(posx<4.9)*(posy<4.9)*(amp>exp(th))*(amp<exp(th2))
    # ~ y = amp[sel]
    # ~ avs[i,:] = [mean(y),std(y)]
    # ~ y = (y-mean(y))/std(y)
    # ~ h = histogram(y,51)
    # ~ yd = h[0]/sum(h[0])
    # ~ wy = sqrt(h[0])/sum(h[0])
    # ~ xd = (h[1][1:]+h[1][:-1])/2.0
    # ~ yd = yd/(xd[1]-xd[0])
    # ~ plot(xd,yd,'.-')
    # ~ gfit = dblgausfit(xd,yd,wy=1.0/wy**2,par0=array([0.5,-1,.5,1,.5]))
    # ~ x = linspace(min(xd),max(xd),101)
    # ~ plot(x,dblgaussd(x,gfit.x),'r--')
    # ~ thr = gfit.x[1]+2*gfit.x[2]
    # ~ if thr>gfit.x[3]-gfit.x[4]:
        # ~ thr = min(0,gfit.x[3]-gfit.x[4])
    # ~ ths[i] = thr
    # ~ vlines(thr,0,max(yd),linestyles='--')
    # ~ fig.show()
    # ~ key = input("Is it OK?")
    # ~ try:
        # ~ ths[i] = float(key)
        # ~ nokarray[i] = 1 
    # ~ except:
        # ~ pass
    
    # ~ fig.clear()
   
nths = ths#*avs[:,1]+avs[:,0]
    
for i,name in enumerate(namesm):
    nameNR = name.split("_")[-2]
    popts = load(name)
    posx = popts[:,5]
    posy = popts[:,6]
    amp = popts[:,0]
    if nths[i]>0:
        th = log(nths[i])
    else:
        th = 0
    th2 = 13.0
    sel = (popts[:,-1]==0)*(posx>0.1)*(posy>0.1)*(posx<4.9)*(posy<4.9)*(amp>exp(th))*(amp<exp(th2))
    
    
    m1,m2 = (mean(amp[sel]),std(amp[sel]))
    amb = (amp-m1)/m2
    am = (amp[sel]-m1)/m2
    h = histogram(am,arange(min(am),max(am),0.2))
    hd = h[0]/sum(h[0])/0.2
    hx = (h[1][1:]+h[1][:-1])/2.0
    dgfit = dblgausfit(hx,hd,par0=array([0.5,-1,.5,1,.5]))
    dgt = concatenate(([m1,m2],dgfit.x))
    s1 = max(dgt[6],dgt[4])
    yt = abs((dgt[3]-dgt[5])/s1)
    
    sel2 = sel*(amb>min(dgt[5]+2*dgt[6],dgt[3]+2*dgt[4]))
    if sel2.sum()<500:
        sel2 = sel*(amb>max(dgt[5]-2*dgt[6],dgt[3]-2*dgt[4]))
    if sel2.sum()>500:
        sel = sel2
    
    
    ts = arange(popts.shape[0])[sel]

    nframes = popts.shape[0]
    if nframes == 12000:
        tsr = ts*5e-3
        ttf = nframes*5e-3
    else:
        ttf = nframes*10e-3
        tsr = ts*10e-3
    
    
    popts= popts.transpose()
        
    mex = mean(popts[5,sel]); sdx = std(popts[5,sel])
    xlims =[ floor(10*(mex-2*sdx)*.98)/10.0,ceil(10*(mex+2*sdx)*1.02)/10.0]
    mey = mean(popts[6,sel]); sdy = std(popts[6,sel])

    ylims =[ floor(10*(mey-2*sdy)*.98)/10.0,ceil(10.0*(mey+2*sdy)*1.02)/10.0]
    
    npixelx = floor(xlims[1]-xlims[0])+1
    npixely = floor(ylims[1]-ylims[0])+1
    npixell = max(npixelx,npixely)
    npixelx = npixell
    npixely = npixell
    xwdthleft = npixelx-(xlims[1]-xlims[0])
    ywdthleft = npixely-(ylims[1]-ylims[0])
    xlims[0] = xlims[0]-xwdthleft/2.0
    xlims[1] = xlims[1]+xwdthleft/2.0
    ylims[0] = ylims[0]-ywdthleft/2.0
    ylims[1] = ylims[1]+ywdthleft/2.0
    
    scdiv = 325/20.0

    nf = 50
    ns = 50
    k = sum(sel)
    sel2 = range(0,sum(sel))

    figure(1)
    xx,yy = (popts[5,sel][sel2],popts[6,sel][sel2])
    xm = running_mean(xx,ns)
    ym = running_mean(yy,ns)
    
    fig = figure(1)
    ax = fig.add_subplot(111)
    sel2 = range(0,k)
    sc = scatter(popts[5,sel][sel2],popts[6,sel][sel2],c=tsr[sel2],cmap=cmapgnu,alpha=0.5)
    xlim(xlims[0],xlims[1])
    ylim(ylims[0],ylims[1])
    axis('off')
    plot([xlims[1]*0.99-1.0/scdiv,xlims[1]*0.99],[ylims[0]*1.01,ylims[0]*1.01],'k-',linewidth=2.0)
    text(xlims[1]*0.99-.70/scdiv,ylims[0]*1.013,"20nm")
    sc.set_clim(0,12000*5e-3)
    axc = colorbar()
    axc.set_label("t(s)")
    xlabel("pixel")
    ylabel("pixel")
    subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
        hspace = 0, wspace = 0)
    margins(0,0)
    gca().xaxis.set_major_locator(NullLocator())
    gca().yaxis.set_major_locator(NullLocator())
    fig.patch.set_facecolor((.4,.4,.4))
    plot(xm,ym,'-',c='C0',alpha=0.7)
    savefig(outputdir+"final"+str(i).zfill(3)+".png",bbox_inches='tight',pad_inches = 0,facecolor=((0.4,0.4,0.4)))
    close(1)


   
    
