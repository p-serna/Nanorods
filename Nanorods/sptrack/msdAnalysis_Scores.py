import sys
import os

wdirs = ["/mnt/data/Anastasia/18_11_29_pd23_11_div6_25Hzsqwave/"]

dfiles = []
for dirt in wdirs:
    basedir = dirt
    files = os.listdir(basedir)
    if dirt[-1] != '/':
        dirt = dirt+'/'
    for f in files:
        if f[-4:]=='.tif': 
            try:
                i = int(f[-5])
                dfiles.append(dirt+f)
            except:
                pass
cfile = dfiles[0]
cfile
cfile = dfiles[1]
wdir = ''
cf2 = cfile.split(".")[0].split("/")
for fs in cf2[:-1]:
    wdir = wdir+fs+'/'
wdir = wdir+cf2[-1]+'output/sptrack/'



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

tmax = data[0][1].shape[0]

t = arange(tmax)


figure()
cmapt = get_cmap('tab20')

counter = 0
k = 0
for i in range(len(data)):
    d = data[i]
    ni,yt,ey2,cnt = d
   
    xtt = t*10.0e-3
    ytt = yt*.325**2
    
    dto =  log(tmax)/200
    tml = log(t[1:])
    ntl = int(tml[-1]/dto)+1
    t0l = 0
    dd = zeros((ntl,2))-1
    dd[0,:2] = (t[1],ytt[1])
    for il in arange(1,ntl):
        til = t0l+dto
        sel = (tml>=t0l)*(tml<til)
        if sel.sum()>0:
            xe = xtt[1:][sel]
            ye = ytt[1:][sel]
            dd[il,:2] = [mean(xe),mean(ye)]
        t0l = til

    dde = dd[dd[:,0]>0,:]
    dde = dde[1:,:]
    
    # ~ figure(6)
    # ~ plot(dde[:,0],dde[:,1]-dde[0,1],'.-',alpha=0.3,label=i)
    # ~ xscale("log")
    # ~ yscale("log")
    
    
    # ~ figure(7)
    # ~ plot(dde[:,0],(dde[:,1]-dde[0,1])/(dde[-1,1]-dde[0,1]),'.-',alpha=0.3,label=i)
    
    if i == 0:
        xde = column_stack((dde,dde[:,:1]*0+i))
        xden = column_stack((dde[:,:1],(dde[:,1:2]-dde[0,1])/(dde[-1,1]-dde[0,1]),dde[:,:1]*0+i))
    else:
        xde = row_stack((xde,column_stack((dde,dde[:,:1]*0+i))))
        xden = row_stack((xden,column_stack((dde[:,:1],(dde[:,1:2]-dde[0,1])/(dde[-1,1]-dde[0,1]),dde[:,:1]*0+i))))
        
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
    
    xt = ds
    xt = xt[xt[:,1]>0,:]
    De = array([mean(xt[xt[:,0]<1e-1,1]),exp(mean(log(xt[xt[:,0]<1e-1,1])))])
    if log10(De[0])< -2.5:       
        # ~ figure(5)
        # ~ #fill_between(ds[:,0],ds[:,1]-ds[:,2],ds[:,1]+ds[:,2],alpha=.5)
        # ~ plot(ds[ds[:,1]>1e-5,0],ds[ds[:,1]>1e-5,1],'.-',alpha=0.3,label=i)
        if counter< 100:
            figure(9)
            plot(dde[:,0],dde[:,1]-dde[0,1],'.-',alpha=0.3,label=i)
            counter +=1
        # ~ if counter>15:
            # ~ break
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

ts = array([100,200,400,800,1600,3200,6400,12800])/1000.0
xt = xde[abs(xde[:,-1]-i)<1e-3,0]
idxts = zeros(ts.shape[0],dtype=int)
idx = arange(xt.shape[0])
for i in range(len(ts)):
    idxts[i] = idx[argmin(abs(xt-ts[i]))]

alphas = zeros((len(data),7))
alphalm = zeros((len(data),2))

Des = zeros((len(data),3))
Dep = zeros((len(data),6))

for i in range(len(data)):
    xt = xds[abs(xds[:,-1]-i)<1e-3,:]
    xt = xt[xt[:,1]>0,:]
    De = array([mean(xt[xt[:,0]<2e-1,1]),mean(xt[xt[:,0]<2e-1,2]),max(xt[xt[:,0]<2e-1,2])])
    sel = (xt[:,0]>5e-1)*(xt[:,0]<1.5e0)
    De2 = array([mean(xt[sel,1]),mean(xt[sel,2]),max(xt[sel,2])])
    sel = (xt[:,0]>5e0)*(xt[:,0]<1.5e1)
    De3 = array([mean(xt[sel,1]),mean(xt[sel,2]),max(xt[sel,2])])
    
    xt = xde[abs(xde[:,-1]-i)<1e-3,:]    
    msds = xt[idxts,:]
    msds = msds[1:,:]/msds[:-1,:]
    alpha = log(msds[:,1])/log(msds[:,0])
    lm = linregress(log(msds[:,0]).cumsum(),alpha)
    alphas[i,:] = alpha
    alphalm[i,:] = (lm.slope,lm.intercept)
    Des[i,:] = De
    Dep[i,:] = concatenate((De2,De3))
    
figure()
plot(Des[:,0],Des[:,1]/Des[:,0],'.',alpha=0.3)
xscale("log")
yscale("log")
plot(Des[:,0],Des[:,2]/Des[:,0],'.',alpha=0.3)
figure(); plot(Des[:,0],Des[:,1]/Des[:,0],'.',alpha=0.3)
xscale("log")

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

amps = []
for ni,name in enumerate(dfiles):
    popts = load(name)
    posx = popts[:,5]
    posy = popts[:,6]
    amp = popts[:,0]
    th = 0
    th2 = 13.0
    sel = (popts[:,-1]==0)*(posx>0.1)*(posy>0.1)*(posx<4.9)*(posy<4.9)*(amp>exp(th))*(amp<exp(th2))
    Tf = len(posx)
    if sel.sum() > 1000:
        amps.append([mean(amp[sel]),mean(popts[sel,4])])
ampst = array(amps)

figure(); plot(Des[:,0],ampst[:,0],'.',alpha=0.3); xscale('log'); yscale('log')
xlabel("D_0")
ylabel("mean Amplitude")
figure(); plot(Des[:,1]/Des[:,0],ampst[:,0],'.',alpha=0.3); xscale('log'); yscale('log')
ylabel("mean Amplitude")
xlabel("$\Delta D_0/D_0$")

figure(); plot(Des[:,0],ampst[:,0]/ampst[:,1],'.',alpha=0.3); xscale('log'); yscale('log')
xlabel("D_0")
ylabel("mean Amplitude/bg")

figure(); plot(Des[:,1]/Des[:,0],ampst[:,0]/ampst[:,1],'.',alpha=0.3); xscale('log'); yscale('log')
ylabel("mean Amplitude / bg")
xlabel("$\Delta D_0/D_0$")

figure();
hist(log10(Des[isfinite(Des[:,1]),1]),51)
xlabel("log10(D_0)")



figure();
hist(alphas[isfinite(alphas[:,0]),0],51)
xlabel("beta")


figure();
hist(alphas[isfinite(alphas[:,0]),0],51)
hist(alphalm[isfinite(alphalm[:,1]),1],51)

xlabel("beta")


figure();
hist(alphalm[isfinite(alphalm[:,0]),0],51)

xlabel("beta")

figure();am = alphas.mean(axis=1); hist(am[isfinite(am)],31)

masked = loadtxt(wdir+"../ROIs_Inmask.dat")


savetxt(outputdir+"ROIsfeatures.dat",column_stack((Des,am,alphalm,alphas)),header='# <Deff(t=0)> <Deff(t=0)>_log <alpha> lm_alpha: slope intercept, alpha(0:8)')












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
