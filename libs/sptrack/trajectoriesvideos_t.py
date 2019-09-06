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
def dblgaussd(x,par): return(par[0]*gaussd(x,par[1:3])+(1.0-par[0])*gaussd(x,par[3:]))

def dblgausfit(x,y,wy=1.0,par0=array([0.5,-1,.5,1,.5])):
    def minf(par): return( sum(wy*(dblgaussd(x,par)-y)**2)/sum(wy))
    minx  = minimize(minf,par0)
    
    return minx
def running_mean(x, N):
    cumsumt = cumsum(insert(x, 0, 0)) 
    return (cumsumt[N:] - cumsumt[:-N]) / float(N)

dirs = ["./sptrack/"]

for dirt in dirs:
    basedir = dirt
    files = os.listdir(basedir)

    dfiles = []
    for f in files:
        if f[-4:]=='.npy'and f[:4]=='pots': dfiles.append(basedir+f)

    arrays = {}
    namesf = []
    namesm = []

    print("Data files:")
    for filepath in dfiles:
        f = load(filepath)
        name = filepath.split(sep=".")[1]
        namesf.append(name)
        namesm.append(name)
        arrays[name] = array(f)

namesm.sort()
cmapgnu = get_cmap("gnuplot")

#ths = load("thresholds.npy")
data = []
print(namesm)

selROI = arange(len(namesm))
try:
    ths = load("sptrack/thresholds.npy")
    if len(ths)<len(namesm):
        ths = concatenate((ths,zeros(len(namesm)-len(ths))))
except:
    ths = zeros(len(namesm))

for i in arange(len(namesm)):
    name = namesm[i]
    popts = arrays[name]
    th = ths[i]
    sel = (popts[-1,:]==0)*(popts[1,:]>0.32)*(popts[2,:]>0.32)*(
        popts[5,:]>0.1)*(popts[6,:]>0.1)*(popts[5,:]<4.9)*(popts[6,:]<4.9)*(popts[0,:]>exp(th))

    a = popts[:,sel]
    ts = arange(12000)[sel]

    tsr = ts*5e-3

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

    for k in range(0,sum(sel),20):   
        sel2 = range(0,k)

        xx,yy = (popts[5,sel][sel2],popts[6,sel][sel2])
        xm = running_mean(xx,20)
        ym = running_mean(yy,20)
        
        fig = figure(1)
        ax = fig.add_subplot(111)
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
        ax.set_facecolor((.4,.4,.4))
        plot(xm,ym,'-',c='C0')
        savefig("movies/im"+str(k).zfill(5)+".png",bbox_inches='tight',pad_inches = 0,facecolor=((0.4,0.4,0.4)))
        close(1)
    figure(1)
    xx,yy = (popts[5,sel][sel2],popts[6,sel][sel2])
    xm = running_mean(xx,20)
    ym = running_mean(yy,20)
    
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
    plot(xm,ym,'-',c='C0')
    savefig("movies/im"+str(0).zfill(5)+".png",bbox_inches='tight',pad_inches = 0,facecolor=((0.4,0.4,0.4)))
    savefig("movies/final"+str(i).zfill(3)+".png",bbox_inches='tight',pad_inches = 0,facecolor=((0.4,0.4,0.4)))
    close(1)


    os.system("cat movies/im*.png | ffmpeg -framerate 11 -i - movies/vid"+str(i).zfill(3)+".mp4")
    os.system("rm movies/im*.png")
    
    
