from scipy.optimize import minimize,least_squares
from numpy import *
from matplotlib.pylab import *

  

def gaussd(x,par): return(exp(-(x-par[0])**2/2/par[1]**2)/sqrt(2*pi*par[1]**2) )
def dblgaussd(x,par): return(par[0]*gaussd(x,par[1:3])+(1.0-par[0])*gaussd(x,par[3:]))

def dblgausfit(x,y,par0=array([0.5,-1,.5,1,.5])):
    def minf(par): return( sum((dblgaussd(x,par)-y)**2))
    minx  = minimize(minf,par0)
    
    return minx

def dblgausfit_tseries(S,bins=41,fname="Distribution",par0=array([0.8,-.5,.5,2,2]),normalization=True):
    h0 = histogram(S,bins=bins)
    hy = h0[0]/sum(h0[0])
    hx = h0[1]
    
    if normalization:
        mS = mean(S)
        sS = sqrt(var(S))
    else:
        mS = 0.0
        sS = 1.0
            
    hxt = ((hx[1:]+hx[:-1])/2.-mS)/sS
    hy = hy/(hxt[1]-hxt[0])
    
    figure(fname)
    plot(hxt,hy,'-o',label="Data")
    
    def minf(par): return( sum((dblgaussd(hxt,par)-hy)**2))
    minx  = minimize(minf,par0)
    plot(x,dblgaussd(x,minx.x),'--',label="Fit double gaussian")
    
    return minx,mS,sS


def movingfit(S,bins=20,nt=2000,fname="Distribution",par0=array([0.8,-.5,.5,2,2]),normalization=True,tol=1e0):
    h0 = histogram(S,bins=bins)
    hy = h0[0]/sum(h0[0])
    hx = h0[1]
    
    if normalization:
        mS = mean(S)
        sS = sqrt(var(S))
    else:
        mS = 0.0
        sS = 1.0
            
    hxt = ((hx[1:]+hx[:-1])/2.-mS)/sS
    hy = hy/(hxt[1]-hxt[0])
    
    def minf(par,hxt,hy): return( sum((dblgaussd(hxt,par)-hy)**2))
    minx  = minimize(minf,par0,args=(hxt,hy))
    
    par0 = minx.x
    part = par0
    para = part
    pars = []
    #max(floor(nt/4),1)
    for i in range(0,nt//2,1):
        pars.append(par0)
    for i in range(0,len(S)-nt,1):
        St = S[i:(i+nt)]
        h0 = histogram(St,bins=bins)
        
        hx = h0[1]
        hxt = ((hx[1:]+hx[:-1])/2.-mS)/sS
        hy = h0[0]/sum(h0[0])/(hxt[1]-hxt[0])
        minxt = minimize(minf,part,args=(hxt,hy))
        part = minxt.x
        if sqrt(sum((part-para)**2)/sum(para**2))>tol:
            minxt = minimize(minf,par0,args=(hxt,hy))
            part = minxt.x
            if sqrt(sum((part-para)**2)/sum(para**2))>tol:
                part = par0
        if abs(part[1]-par0[1])>par0[2] or abs(part[3]-par0[3])>par0[4]:
            part = par0
        pars.append(part)
    for i in range(0,nt//2,1):
        pars.append(par0)
            
    return array(pars)




def gaussian_fits(Sn,bins=linspace(-4,4,101)):
    Sshape = Sn.shape
    hs = []
    pars = []
    xb = bins
    xh = (xb[1:]+xb[:-1])/2.0
    dx = xb[1]-xb[0]
    for i in arange(Sshape[0]):
        ht = histogram(Sn[i,:],xb)[0]
        hs.append(ht)
        yh = ht/sum(ht)/dx
        plot(xh,yh,label=i)
        minx = dblgausfit(xh,yh)
        print(minx.x)
        pars.append(minx.x)
    return pars



def gaussian_fit(Sn,bins=linspace(-4,4,101),verbose=False,par0=array([0.5,-1,.5,1,.5])):
    xb = bins
    xh = (xb[1:]+xb[:-1])/2.0
    dx = xb[1]-xb[0]
    ht = histogram(Sn,xb)[0]
    yh = ht/sum(ht)/dx
    if verbose: plot(xh,yh,label=i)
    minx = dblgausfit(xh,yh,par0=par0)
    if verbose: print(minx.x)
    return minx.x
    
def isbimodal(par):
    p,m1,s1,m2,s2 = par
    d = sqrt((m1-m2)**2/s1/s2/4.0)
    A1 = (d<=1)
    if d>1:
        A2 = abs(log(1-p)-log(p))>=2*log(d-sqrt(d**2-1))+2*d*sqrt(d**2-1)
    else:
        A2 = False
    return not(A1 or A2)

def get_intervals_fixedD(Sn,par,probl = 0.1):
    threshold = par[1]+par[2]*2.0
    prob = sum(Sn>threshold)/len(Sn)
    intervals = []
    if prob>probl:
        Stn = array(Sn)
        sel = Stn<threshold
        Stn[sel] = 0.0
        csel = cumsum(1-sel)
        selt = 0
        ica = 0
        on = True
        for ic,c in enumerate(csel):
            if c == selt:
                if not on:
                    intervals.append([ica,ic])
                    ica = ic
                on = True
            else:
                if on:
                    intervals.append([ica,ic])
                    ica = ic
                    selt = c
                selt = c
                on = False
        return intervals
    else:
        print("Not enough points in the second gaussian!")
        return nan

def get_intervals_fixedth(Sn,pth,probl = 0.1):
    threshold = pth
    prob = sum(Sn>threshold)/len(Sn)
    intervals = []
    if prob>probl:
        Stn = array(Sn)
        sel = Stn<threshold
        Stn[sel] = 0.0
        csel = cumsum(1-sel)
        selt = 0
        ica = 0
        on = True
        for ic,c in enumerate(csel):
            if c == selt:
                if not on:
                    intervals.append([ica,ic])
                    ica = ic
                on = True
            else:
                if on:
                    intervals.append([ica,ic])
                    ica = ic
                    selt = c
                selt = c
                on = False
        return intervals
    else:
        print("Not enough points in the second gaussian!")
        return nan
        
def fixing_intervals(s0t,s1t,length=1,onlyup=True):
    s0d = s0t[:,1]-s0t[:,0]
    sel0 = s0d <= length
    s1d = s1t[:,1]-s1t[:,0]
    sel1 = s1d <= length
    #s0t[sel1,1] = s1t[sel1,1]
    s1t = s1t[bitwise_not(sel1),:]
    if onlyup:
        return s0t,s1t
    else:
        print("Not implemented the two sideways")
        return s0t,s1t

def fixing_intervals2(x,length=1,onlyup=True):
    xa = array(x)
    s0t = xa[arange(0,xa.shape[0],2),:]
    s1t = xa[arange(1,xa.shape[0],2),:]
    s0d = s0t[:,1]-s0t[:,0]
    sel0 = s0d <= length
    s1d = s1t[:,1]-s1t[:,0]
    sel1 = arange(len(s1d))[s1d <= length]
    if len(sel1) > len(s0t): 
        sel1 = sel1[0:len(s0t)]
        s0t[sel1,1] = s1t[sel1,1]
    elif len(sel1)<len(s0t):
        s0t[sel1,1] = s1t[sel1,1]
    sel1 = s1d<=length    
        

    s1t = s1t[bitwise_not(sel1),:]
    if onlyup:
        return s0t,s1t
    else:
        print("Not implemented the two sideways")
        return s0t,s1t
        
def expanalysis(s0A,s1A,xmin = 1, nbins=51):
    xbins = linspace(0,200,nbins)
    dxbins = (xbins[1]-xbins[0])
    xbins = xbins-(xbins[1]-xbins[0])/2.0
    xbinb = (xbins[1:]+xbins[:-1])/2.0
    h0 = histogram(s0A[:,1]-s0A[:,0],bins=xbins)[0]
    h1 = histogram(s1A[:,1]-s1A[:,0],bins=xbins)[0]
    plot(xbinb,h0/sum(h0)/dxbins,'.-',label="Off state")
    plot(xbinb,h1/sum(h1)/dxbins,'.-',label="On State")
    yscale("log")

    ht = h0/sum(h0)/dxbins; 
    sf = arange(xmin,len(ht))[ht[xmin:]<1e-3][0]
    hy = log(ht[xmin:sf]); hx = xbinb[xmin:sf]
    ps = polyfit(hx,hy,deg=1)
    xseq = linspace(0,50,113)
    lambd = -1./ps[0]*5
    plot(xseq,exp(xseq*ps[0]+ps[1]),'b--',label="$\lambda =%.2f $ ms" % lambd)
    ps0 = ps
    ht = h1/sum(h1)/dxbins; 
    sf = arange(xmin,len(ht))[ht[xmin:]<1e-3][0]
    hy = log(ht[xmin:sf]); hx = xbinb[xmin:sf]
    ps = polyfit(hx,hy,deg=1)
    xseq = linspace(0,50,113)
    lambd = -1./ps[0]*5
    plot(xseq,exp(xseq*ps[0]+ps[1]),'r--',label="$\lambda = %.2f $ ms" % lambd)
    ps1 = ps
    ylim(1e-4,1)
    return ps0,ps1


def histogram_wmemory(x,bins):    
    if isscalar(bins):
        nb = bins-1
        if bins>1:
            minx = min(x)*0.999
            maxx = max(x)*1.001
            dx = (maxx-minx)/(nb)
            xh = arange(minx,maxx,dx)
    else:
        nb = len(bins)-1
        xh = bins
        minx = min(bins)
        maxx = max(bins)
        dx = bins[1]-bins[0]
    h = zeros(nb)
    cache =  {i:[] for i in range(nb)}
    warningOuth = False
    for i,xt in enumerate(x):
        ih = int((xt-minx)/dx)
        if ih<nb:
            h[ih] += 1
            cache[ih].append(i)
        else: 
            warningOuth = True
    return h,xh,cache


        



# Not working correctly. Padding not well implemented
def movingfitwpad(Sr,bins=20,nt=2000,fname="Distribution",par0=array([0.8,-.5,.5,2,2]),normalization=True,tol=1e0,padding="valid"):
    h0 = histogram(Sr,bins=bins)
    hy = h0[0]/sum(h0[0])
    hx = h0[1]
    
    nt2 = int(nt/2)
    if padding == "valid":
        S = array(Sr)
    elif padding == "same":
        S = concatenate((arange(nt/2)*0,Sr,arange(nt/2))).flatten()
        for i in range(nt2):
            S[i] = mean(S[(i+nt2):(i+3*nt2)])
            S[-i-1] = mean(S[(-i-1-nt2):(-i-1-3*nt2)])
    
    if normalization:
        mS = mean(S)
        sS = sqrt(var(S))
    else:
        mS = 0.0
        sS = 1.0
            
    hxt = ((hx[1:]+hx[:-1])/2.-mS)/sS
    hy = hy/(hxt[1]-hxt[0])
    
    def minf(par,hxt,hy): return( sum((dblgaussd(hxt,par)-hy)**2))
    minx  = minimize(minf,par0,args=(hxt,hy))
    
    par0 = minx.x
    part = par0
    para = part
    pars = []
    #max(floor(nt/4),1)
    for i in range(0,len(S)-nt,1):
        St = S[i:(i+nt)]
        h0 = histogram(St,bins=bins)
        
        hx = h0[1]
        hxt = ((hx[1:]+hx[:-1])/2.-mS)/sS
        hy = h0[0]/sum(h0[0])/(hxt[1]-hxt[0])
        minxt = minimize(minf,part,args=(hxt,hy))
        part = minxt.x
        if sqrt(sum((part-para)**2)/sum(para**2))>tol:
            minxt = minimize(minf,par0,args=(hxt,hy))
            part = minxt.x
            if sqrt(sum((part-para)**2)/sum(para**2))>tol:
                part = par0
            
            
        pars.append(part)
    
    # ~ if padding == "same":
        # ~ pars
    return array(pars)

def timeseries(t,S,figname="",figsize=(8,6),nbins=31):
    nullfmt = NullFormatter()         # no labels
    
    # definitions for the axes
    left, width = 0.08, 0.75
    bottom, height = 0.1, 0.8
    bottom_h = left_h = left + width + 0.02

    rect_scatter = [left, bottom, width, height]
    rect_histy = [left_h, bottom, 0.13, height]

    # start with a rectangular Figure
    figure(figname, figsize=figsize)

    axScatter = axes(rect_scatter)
    axHisty = axes(rect_histy)

    # no labels
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the scatter plot:
    axScatter.plot(t, S,'.-',alpha=0.3)

    # now determine nice limits by hand:
    binwidth = 0.25
    axScatter.set_ylim((min(S),max(S)))

    bins = linspace(min(S), max(S), nbins)
    h = histogram(S,bins=bins)
    dx = h[1][1]-h[1][0]
    axHisty.plot(h[0]/sum(h[0])/dx,(h[1][1:]+h[1][:-1])/2.0,alpha=0.4)
    axHisty.fill_betweenx((h[1][1:]+h[1][:-1])/2.0,h[0]/sum(h[0])/dx,alpha=0.4)

    axHisty.set_ylim(axScatter.get_ylim())
    axHisty.set_xlim((0,max(h[0]/sum(h[0])/dx)))
    axHisty.patch.set_visible(False)
    axHisty.axis('off')
    return (axScatter,axHisty)
