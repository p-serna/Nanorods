from numpy import *
from matplotlib.pylab import *
import scipy.optimize as opt
import pickle
from scipy.optimize import minimize,least_squares
from scipy.stats import linregress

def f0(x,par): return(par[0]+par[1]*x)
def f0e(x,par): return(par[0]+par[1]*x**par[2])
def funcfit(fun,x,y,par0,ey=1.0):
    def minf(par): return( sum((fun(x,par)-y)**2/ey**2))
    minx  = minimize(minf,par0)
    return minx

def minf0(par): return( sum((f0(xt,par[:2])-yt)**2/ey**2))
def minf0e(par): return( sum((f0e(xt,par[:2])-yt)**2/ey**2))
def minf1(par):
    ey2 = ey/yt  
    return( sum((f0(log(xt),par[:2])-log(yt))**2/ey2**2))
def minf2(par):
    if sum(yt>par[2])>0:
        return(inf)
    ey2 = ey/abs(par[2]-yt) 
    return( sum((f0(xt,par[:2])-log(par[2]-yt))**2/ey2**2))
def minf3(par):
    if sum(yt>par[2])>0:
        return(inf)
    ey2 = ey/abs(par[2]-yt)/abs(log(par[2]-yt))  
    return( sum((f0(xt,par[:2])-log(par[2]-yt))**2/ey2**2))



if len(sys.argv)>1:
    wdir = sys.argv[1]
    if wdir[-1] != '/':
        wdir = wdir+'/'
    if len(sys.argv)>2:
        interactive = sys.argv[2]=='1'
    else:
        interactive = False
else:
    wdir = "./"
    interactive = False
    print("No argument, we do it for current folder")




try:
    with open(wdir+'data_msd.pickle', 'rb') as handle:
        data = pickle.load(handle)
except:
    print("No data_msd pickle in the folder!")
    data = []
    return

# This is to see if ax^c+b is better than ax+b

fits = []
for d in data:    
    tmax = 5000
    t = arange(tmax)
    ni,yt,ey2,cnt = d

    #plot(t,yt,'.-')
    sel = (t>5)*(t<4000)
    yt = yt[sel]
    ey2 = ey2[sel]
    xt = t[sel]
    par0 = array([yt[0],(yt[-1]-yt[0])/(xt[-1]-xt[0])])

    fitrw = funcfit(f0,xt,yt,par0,ey2)
    parrw = fitrw.x
    par0 = array([parrw[0]*(1.0+0.01*randn()),parrw[1]*(1.0+0.01*randn()),0.98])
    fitaw = funcfit(f0e,xt,yt,par0,ey2)
    fits.append([fitrw.fun,fitrw.x,fitaw.fun,fitaw.x])

fs = []
ps = []
for f in fits:
    f1,p1,f2,p2 = f
    fs.append([f1,f2])
    ps.append([p1[0],p1[1],p2[0],p2[1],p2[2]])

fs = array(fs)
ps = array(ps)

# ~ sel = isfinite(fs[:,1])
# ~ figure(); hist(ps[sel,4],21)
# ~ figure(); plot(ps[sel,4],'.')
# ~ figure(); plot(ps[sel,2])

def running_mean(x, N):
    cumsumt = cumsum(insert(x, 0, 0)) 
    return( (cumsumt[N:] - cumsumt[:-N]) / float(N))
    
def running_slope(x,y,N):
    cumsumx = cumsum(insert(x, 0, 0))
    xm =  (cumsumx[N:] - cumsumx[:-N]) / float(N)
    cumsumx2 = cumsum(insert(x**2, 0, 0))
    xm2 =  (cumsumx2[N:] - cumsumx2[:-N]) / float(N)
    varx = xm2-xm**2
    cumsumy = cumsum(insert(y, 0, 0))
    ym =  (cumsumy[N:] - cumsumy[:-N]) / float(N)
    cumsumxy = cumsum(insert(x*y, 0, 0))
    xym =  (cumsumxy[N:] - cumsumxy[:-N]) / float(N)
    covxy = xym-xm*ym
    return(covxy/varx)

figure()
cmapt = get_cmap('tab20')
# ~ k = 0
# ~ ons = []
# ~ for i in arange(len(data)):
    # ~ d = data[i]
    # ~ ni,yt,ey2,cnt = d
    
    # ~ ytt = yt*.325**2
    
    # ~ cv = var(ytt)/mean(ytt)**2
    # ~ if cv<0.5:
        # ~ ons.append(1)    
        # ~ plot(t*5,(sqrt(ytt)-sqrt(ytt[1]))/sqrt(max(ytt)),'.-',c=cmapt(k%20),alpha=0.9,markersize=12,label=cv)
        # ~ k+=1
    # ~ else:
        # ~ ons.append(0)



k = 0
for i in range(len(data)):
    d = data[i]
    ni,yt,ey2,cnt = d
   
    xtt = t*5
    ytt = yt*.325**2
    
    dto =  log(5000)/200
    tml = log(t[1:])
    ntl = int(tml[-1]/dto)+1
    t0l = 0
    dd = zeros((ntl,2))
    dd[0,:2] = (t[1],ytt[1])
    for il in arange(1,ntl):
        til = t0l+dto
        sel = (tml>=t0l)*(tml<til)
        xe = t[1:][sel]*5.0e-3
        ye = ytt[1:][sel]
        dd[il,:2] = [mean(xe),mean(ye)]
        t0l = til

    dde = dd[isfinite(dd[:,0]),:]
    dde = dde[1:,:]
    figure(6)
    plot(dde[:,0],dde[:,1]-dde[0,1],'.-',alpha=0.3,label=i)
    xscale("log")
    yscale("log")
    figure(7)
    plot(dde[:,0],(dde[:,1]-dde[0,1])/(dde[-1,1]-dde[0,1]),'.-',alpha=0.3,label=i)
    if i == 0:
        xde = column_stack((dde,dde[:,:1]*0+i))
        xden = column_stack((dde[:,:1],dde[:,1:2]/dde[-1,1],dde[:,:1]*0+i))
    else:
        xde = row_stack((xde,column_stack((dde,dde[:,:1]*0+i))))
        xden = row_stack((xden,column_stack((dde[:,:1],dde[:,1:2]/dde[-1,1],dde[:,:1]*0+i))))
        
    # ~ ds = zeros((ntl,3))
    # ~ nwin = 20
    # ~ for il in range(dde.shape[0]-2):
        # ~ sel = arange(il-nwin+3,il+3)
        # ~ sel = sel[sel>=0]
        # ~ xe = dde[sel,0]
        # ~ ye = dde[sel,1]    
        # ~ lm = linregress(xe,ye)
        # ~ ds[il,0] = mean(xe)    
        # ~ ds[il,1:3] = (lm.slope,lm.stderr)    
        
    ntle = dde.shape[0]
    nwin = 50
    ds = zeros((ntle+nwin-3,3))
    for il in range(ntle+nwin-3):
        sel = arange(il-nwin+3,il+3)
        sel = sel[(sel>=0)*(sel<ntle)]
        xe = dde[sel,0]
        ye = dde[sel,1]    
        lm = linregress(xe,ye)
        ds[il,0] = mean(xe)    
        ds[il,1:3] = (lm.slope,lm.stderr)    
    print(i)
    figure(5)
    #fill_between(ds[:,0],ds[:,1]-ds[:,2],ds[:,1]+ds[:,2],alpha=.5)
    plot(ds[ds[:,1]>1e-5,0],ds[ds[:,1]>1e-5,1],'.-',alpha=0.3,label=i)
    if i == 0:
        xds = column_stack((ds,ds[:,:1]*0+i))
    else:
        xds = row_stack((xds,column_stack((ds,ds[:,:1]*0+i))))
        
xds = array(xds)
xde = array(xde)
xden = array(xden)

yscale("log")
xscale("log")
xlabel("t(s)")
ylabel("$D_{eff}$ ($\mu$ m$^2$/s)")

dsF = column_stack((ds[:,0],ds[:,0]*0,ds[:,0]*0,ds[:,0]*0))


for i in range(len(data)):
    d = data[i]
    ni,yt,ey2,cnt = d


    
    xtt = t*5
    ytt = yt*.325**2
    
    dto =  log(5000)/200
    tml = log(t[1:])
    ntl = int(tml[-1]/dto)+1
    t0l = 0
    dd = zeros((ntl,2))
    dd[0,:2] = (t[1],ytt[1])
    for il in arange(1,ntl):
        til = t0l+dto
        sel = (tml>=t0l)*(tml<til)
        xe = t[1:][sel]*5.0e-3
        ye = ytt[1:][sel]
        dd[il,:2] = [mean(xe),mean(ye)]
        t0l = til

    dde = dd[isfinite(dd[:,0]),:]
    dde = dde[1:,:]
    # ~ ds = zeros((ntl,3))
    # ~ nwin = 20
    # ~ for il in range(dde.shape[0]-2):
        # ~ sel = arange(il-nwin+3,il+3)
        # ~ sel = sel[sel>=0]
        # ~ xe = dde[sel,0]
        # ~ ye = dde[sel,1]    
        # ~ lm = linregress(xe,ye)
        # ~ ds[il,0] = mean(xe)    
        # ~ ds[il,1:3] = (lm.slope,lm.stderr)    
        
    ntle = dde.shape[0]
    nwin = 50

    ds = zeros((ntle+nwin-3,3))
    for il in range(ntle+nwin-3):
        sel = arange(il-nwin+3,il+3)
        sel = sel[(sel>=0)*(sel<ntle)]
        xe = dde[sel,0]
        ye = dde[sel,1]    
        lm = linregress(xe,ye)
        ds[il,0] = mean(xe)    
        ds[il,1:3] = (lm.slope,lm.stderr)    

    # ~ if ons[i]==1:
        # ~ print(i)
    figure(5)
    #fill_between(ds[:,0],ds[:,1]-ds[:,2],ds[:,1]+ds[:,2],alpha=.5)
    
    plot(ds[ds[:,1]>1e-5,0],ds[ds[:,1]>1e-5,1],'.-',alpha=0.3)
    dsF[ds[:,1]>1e-5,1] += log(ds[ds[:,1]>1e-5,1])
    dsF[ds[:,1]>1e-5,2] += (log(ds[ds[:,1]>1e-5,1]))**2
    dsF[ds[:,1]>1e-5,3] += 1

yscale("log")
xscale("log")
xlabel("t(s)")
ylabel("$D_{eff}$ ($\mu$ m$^2$/s)")

dsF[:,1] = dsF[:,1]/dsF[:,3]
dsF[:,2] = (sqrt(dsF[:,2]/dsF[:,3])-dsF[:,1]**2)/sqrt(dsF[:,3])
dsF[:,2] = dsF[:,2]*exp(dsF[:,1])/dsF[:,1]
dsF[:,1] = exp(dsF[:,1])

plot(dsF[:-1,0],dsF[:-1,1],'k-',linewidth = 1.5,alpha=0.9)
fill_between(dsF[:,0],dsF[:,1]-dsF[:,2],dsF[:,1]+dsF[:,2],alpha=.9,color='k')

for i in range(10): close()


dsFn = column_stack((ds[:,0],ds[:,0]*0,ds[:,0]*0,ds[:,0]*0))

k = 0
MSDa = []
Defa = []
for i in range(len(data)):
    d = data[i]
    ni,yt,ey2,cnt = d


    
    xtt = t*5
    ytt = yt*.325**2
    
    dto =  log(5000)/200
    tml = log(t[1:])
    ntl = int(tml[-1]/dto)+1
    t0l = 0
    dd = zeros((ntl,2))
    dd[0,:2] = (t[1],ytt[1])
    for il in arange(1,ntl):
        til = t0l+dto
        sel = (tml>=t0l)*(tml<til)
        xe = t[1:][sel]*5.0e-3
        ye = ytt[1:][sel]
        dd[il,:2] = [mean(xe),mean(ye)]
        t0l = til

    dde = dd[isfinite(dd[:,0]),:]
    dde = dde[1:,:]
    # ~ ds = zeros((ntl,3))
    # ~ nwin = 20
    # ~ for il in range(dde.shape[0]-2):
        # ~ sel = arange(il-nwin+3,il+3)
        # ~ sel = sel[sel>=0]
        # ~ xe = dde[sel,0]
        # ~ ye = dde[sel,1]    
        # ~ lm = linregress(xe,ye)
        # ~ ds[il,0] = mean(xe)    
        # ~ ds[il,1:3] = (lm.slope,lm.stderr)    
    MSDa.append(dde)
    
    ntle = dde.shape[0]
    ds = zeros((ntle+nwin-3,3))
    nwin = 50
    for il in range(ntle+nwin-3):
        sel = arange(il-nwin+3,il+3)
        sel = sel[(sel>=0)*(sel<ntle)]
        xe = dde[sel,0]
        ye = dde[sel,1]    
        lm = linregress(xe,ye)
        ds[il,0] = mean(xe)    
        ds[il,1:3] = (lm.slope,lm.stderr)    
    Defa.append(ds)
    
    if ons[i]==1:
        print(i)
        figure(5)
        #fill_between(ds[:,0],ds[:,1]-ds[:,2],ds[:,1]+ds[:,2],alpha=.5)
        
        plot(ds[ds[:,1]>1e-5,0],ds[ds[:,1]>1e-5,1]/max(ds[:,1]),'.-',alpha=0.9,markersize=12)
        dsFn[ds[:,1]>1e-5,1] += log(ds[ds[:,1]>1e-5,1]/max(ds[:,1]))
        dsFn[ds[:,1]>1e-5,2] += (log(ds[ds[:,1]>1e-5,1]/max(ds[:,1])))**2
        dsFn[ds[:,1]>1e-5,3] += 1

yscale("log")
xscale("log")
xlabel("t(s)")
ylabel("$D_{eff}$ (normalized)")
savefig("DiffusionCoeffTnorm.png",transparent=True)

dsFn[:,1] = dsFn[:,1]/dsFn[:,3]
dsFn[:,2] = (sqrt(dsFn[:,2]/dsFn[:,3])-dsFn[:,1]**2)/sqrt(dsFn[:,3])
dsFn[:,2] = dsFn[:,2]*exp(dsFn[:,1])/dsFn[:,1]
dsFn[:,1] = exp(dsFn[:,1])

plot(dsFn[:-1,0],dsFn[:-1,1],'k-',linewidth = 2,alpha=0.9)
ts = 10**linspace(-2,1.3); plot(ts,(ts/1e-2)**-0.7,'k--' )
savefig("DiffusionCoeffTnorm2.png",transparent=True)

close()


figure()
cmapt = get_cmap('tab20')
k = 0
for i in arange(len(data))[array(ons)==1]:
    d = MSDa[i]
    
    ytt = sqrt(d[:,1])
    plot(d[:,0],1e3*(ytt-ytt[1]),'.-',c=cmapt(k%20),alpha=0.8,label=i)
    k+=1
  
xlabel("t(s)")
ylabel("$\sqrt{<MSD>}-\sqrt{A_0}$(nm)")
xscale("log")
yscale("log")
ylim(2e-1,1.1e3)
#plot(t[t<25],t[t<25]*0+400,'k--',alpha=0.6)
savefig("sMSDmA0smoothed.png",transparent=True,dpi=800)
close()

selection = array([43,13,29,22,19])
figure()
cmapt = get_cmap('tab20')
k = 0
for i in arange(len(data))[array(ons)==1]:
    if i in selection:
        d = MSDa[i]
        
        ytt = d[:,1]
        plot(d[:,0],(ytt-ytt[1]),'.-',c=cmapt(k%20),alpha=0.8,label=i)
        xlabel("t(s)")
        ylabel("$\langle MSD\\rangle$($\mu$m$^2$)")
        xlim(-0.1,5)
        savefig("MSDmA0_"+str(i).zfill(2)+".png",transparent=True,dpi=800)
        

    k+=1


  


figure()
cmapt = get_cmap('tab20')
k = 0
for i in arange(len(data))[array(ons)==1]:
    d = MSDa[i]
    
    ytt = sqrt(d[:,1])
    plot(d[:,0],(ytt-ytt[1])/ytt[-1],'.-',c=cmapt(k%20),alpha=0.8,label=i)
    k+=1
  
xlabel("t(s)")
ylabel("$\sqrt{<MSD>}-\sqrt{A_0}$(normalised)")
xscale("log")
yscale("log")
ylim(1e-3,1.1)
ts = 10**linspace(-2,1.2)
plot(ts,(ts/25.0)**0.5,'k--',alpha=0.8)
text(0.4,0.07,"Random\n Walk")
text(2,0.07,"Ballistic?")

savefig("sMSDmA0smoothedN.png",transparent=True,dpi=800)
close()



figure(6)
#fill_between(ds[:,0],ds[:,1]-ds[:,2],ds[:,1]+ds[:,2],alpha=.5)
plot(dd[1:,0],dd[1:,1],'.-',alpha=0.3)
ylabel("MSD = <R^2> ($\mu$ m$^2$)")
xlabel("t(s)")
savefig("MSDpart3.png")



figure(7)
#fill_between(ds[:,0],ds[:,1]-ds[:,2],ds[:,1]+ds[:,2],alpha=.5)
plot(dd[1:,0],dd[1:,1],'.-',alpha=0.3)
ylabel("MSD = <R^2> ($\mu$ m$^2$)")
xlabel("t(s)")
xscale("log")
xscale("log")
savefig("MSDpart3dlog.png")





figure(8)
#fill_between(ds[:,0],ds[:,1]-ds[:,2],ds[:,1]+ds[:,2],alpha=.5)
plot(dd[2:,0],sqrt(dd[2:,1])-sqrt(dd[1,1]),'.-',alpha=0.3)
title("loc. precision = "+str(int(sqrt(dd[1,1])*1000))+" nm")
ylabel("sqrt(<R^2>)-loc. prec ($\mu$ m$)")
xlabel("t(s)")
xscale("log")
yscale("log")
xs = 10**linspace(-2,1.2,100)
cf =  0.5; plot(xs,xs**cf*0.0563016)
savefig("sMSDpart3dlog.png")
