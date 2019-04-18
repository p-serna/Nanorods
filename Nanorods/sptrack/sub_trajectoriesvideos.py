from numpy import *
from matplotlib.pylab import *
from scipy.optimize import minimize
from matplotlib.collections import LineCollection

def gaussd(x,par): return(exp(-(x-par[0])**2/2/par[1]**2)/sqrt(2*pi*par[1]**2) )
def dblgaussd(x,par): return(par[0]**2*gaussd(x,par[1:3])+(1.0-par[0]**2)*gaussd(x,par[3:]))

def dblgausfit(x,y,wy=1.0,par0=array([sqrt(0.5),-1,.5,1,.5])):
    def minf(par): return( sum(wy*(dblgaussd(x,par)-y)**2)/sum(wy))
    minx  = minimize(minf,par0)
    return minx
    
    
def running_mean(x, N):
    cumsumt = cumsum(insert(x, 0, 0)) 
    return (cumsumt[N:] - cumsumt[:-N]) / float(N)

    
def plotNR(popts, th = 0, th2 = 15., ns = 50, fs = 10e-3, cmap = get_cmap("gnuplot")):
    posx = popts[:,5]
    posy = popts[:,6]
    amp = popts[:,0]
    
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
    
    if sel.sum()>500:
    
        ts = arange(popts.shape[0])[sel]

        nframes = popts.shape[0]
        ttf = nframes*fs
        tsr = ts*fs
        
        
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

        k = sum(sel)
        sel2 = range(0,sum(sel))
        xx,yy = (popts[5,sel][sel2],popts[6,sel][sel2])
        xm = running_mean(xx,ns)
        ym = running_mean(yy,ns)
        
        fig = figure(1)
        ax = fig.add_subplot(111)
        sel2 = range(0,k)
        sc = scatter(popts[5,sel][sel2],popts[6,sel][sel2],c=tsr[sel2],cmap=cmap,alpha=0.5)
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
        #savefig(outputdir+"final"+str(i).zfill(3)+".png",bbox_inches='tight',pad_inches = 0,facecolor=((0.4,0.4,0.4)))
        #close(1)
    return(fig)

   
    
