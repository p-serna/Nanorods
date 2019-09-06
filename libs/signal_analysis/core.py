from scipy.optimize import minimize,least_squares
from numpy import *
from matplotlib.pylab import *

  

def gaussd(x,par): return(exp(-(x-par[0])**2/2/par[1]**2)/sqrt(2*pi*par[1]**2) )
def dblgaussd(x,par): return(par[0]*gaussd(x,par[1:3])+(1.0-par[0])*gaussd(x,par[3:]))

def dblgausfit(x,y,par0=array([0.5,-1,.5,1,.5])):
    def minf(par): return( sum((dblgaussd(x,par)-y)**2))
    minx  = minimize(minf,par0)
    
    return minx
def gaussian_fit(Sn,bins=linspace(-4,4,101),verbose=False):
    xb = bins
    xh = (xb[1:]+xb[:-1])/2.0
    dx = xb[1]-xb[0]
    ht = histogram(Sn,xb)[0]
    yh = ht/sum(ht)/dx
    if verbose: plot(xh,yh,label=i)
    minx = dblgausfit(xh,yh)
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
