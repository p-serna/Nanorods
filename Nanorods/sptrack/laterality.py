from numpy import *
from matplotlib.pylab import *
import os
#from core_libs import *
import scipy.optimize as opt
from scipy.optimize import minimize
import pickle

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


dirt = wdir
basedir = dirt
files = os.listdir(basedir)

dfiles = []
for f in files:
    if f[-4:]=='.npy'and f[:4]=='posh': dfiles.append(basedir+f)

dfiles.sort()



# ~ if interactive:
    # ~ ion()
    # ~ fig = figure(1)
    # ~ show()
# ~ ths = []
# ~ for cfile in dfiles:
    # ~ popts = load(cfile)
    # ~ posx = popts[:,5]
    # ~ posy = popts[:,6]
    # ~ amp = popts[:,0]
    # ~ if interactive:
        # ~ h = hist(log(amp),arange(log(min(amp))-0.02,log(max(amp))+0.02,0.02))
    # ~ else:
        # ~ h = histogram(log(amp),arange(log(min(amp))-0.02,log(max(amp))+0.02,0.02))
    
    # ~ yh = h[0]/sum(h[0]); xh = h[1]; xh = .5*(xh[1:]+xh[:-1]); wy = 1.0/(sqrt(h[0])/sum(h[0]))**2
    # ~ sel = yh>0
    # ~ if sum(sel)>0:
        # ~ yh = yh[sel]
        # ~ xh = xh[sel]
        # ~ wy = wy[sel]
    # ~ me = sum(xh*yh)
    # ~ se = sqrt(sum(xh**2*yh)-me**2)
    # ~ yh = yh/0.05
    # ~ wy = wy*0.05**2
    # ~ par0 = array([0.5,me-se,se/2.0,me+se,se/2.0])
    # ~ gfit = dblgausfit(xh,yh,wy,par0)
    # ~ par = gfit.x
    # ~ xs = arange(xh[0],xh[-1],2e-3)
    # ~ if interactive:
        # ~ xs = arange(xh[0],xh[-1],2e-3)
        # ~ plot(xs,dblgaussd(xs,par)*sum(h[0])*0.05,'k-')
        # ~ plot(xs,gaussd(xs,par[1:3])*par[0]*sum(h[0])*0.05,'b--')
        # ~ plot(xs,gaussd(xs,par[3:])*(1-par[0])*sum(h[0])*0.05,'r--')
    
    # ~ # If they are two gaussians:
    # ~ if min(par[1],par[3])>10.5:
        # ~ th = 10.5
    # ~ else:
        # ~ th = par[1]+par[2]
    # ~ if sum(yh[xh>th])<0.5:
        # ~ th = me
        
    
    # ~ if interactive:
        # ~ vlines(th,0,max(h[0]),'k') 
        # ~ fig.show()
        # ~ thok = input("Threshold")
        # ~ if thok != "":
            # ~ th = float(thok)
        # ~ fig.clear()
    
        
    # ~ #sel = (popts[-1,:]==0)*(popts[1,:]>0.32)*(popts[2,:]>0.32)*(
    # ~ #    popts[5,:]>0.1)*(popts[6,:]>0.1)*(popts[5,:]<4.9)*(popts[6,:]<4.9)*(popts[0,:]>exp(th))
    # ~ ths.append(th)
# ~ ths = array(ths)
# ~ save(wdir+"thresholds.npy",ths)
    
data = []
for ni,name in enumerate(dfiles):
    popts = load(name)
    posx = popts[:,5]
    posy = popts[:,6]
    amp = popts[:,0]
    th = 0
    sel = (popts[:,-1]==0)*(posx>0.1)*(posy>0.1)*(posx<4.9)*(posy<4.9)*(amp>exp(th))
    Tf = len(posx)
    
    ts = arange(Tf)[sel]
    Tfn = ts[-1]
    #posx = posx[sel]
    #posy = posy[sel]
    lp = sum(sel)

    tmax = 5000
    msd = zeros(tmax)
    msd2 = zeros(tmax)
    cnt = zeros(tmax)
    dx = posx[1:]-posx[:-1]
    dy = posy[1:]-posy[:-1]
    dxn = dx/sqrt(dx**2+dy**2)
    dyn = dy/sqrt(dx**2+dy**2)
    #plot(dxn,dyn,'.-')
    figure()
    hist(arctan2(dyn,dxn),31)
    figure()
    arrow(dxn*0,dyn*0,dx,dy,alpha=0.2)
    #t = arange(tmax)

    msd01 = array(msd/cnt)
    msd01b = array(msd2/cnt)
    cnt01 = array(cnt)
    data.append([ni,msd01,msd01b,cnt01])
    
    # ~ figure(2)

    # ~ plot(t[2:],sqrt((msd01[2:]-msd01[1])/msd01[-1]),alpha=0.5)
    # ~ figure(3)

    # ~ plot(t[2:],(sqrt(msd01[2:])-sqrt(msd01[1]))/sqrt(msd01[-1]),alpha=0.5)
    print("File ",ni," of ", len(dfiles))
# ~ xscale("log")
# ~ yscale("log")


with open(wdir+'data_msd.pickle', 'wb') as handle:
    pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

# ~ with open('data.pickle', 'rb') as handle:
    # ~ b = pickle.load(handle)

# ~ idcs = []
# ~ for i,n in enumerate(namesm):
    # ~ if n[-7:-3] == '_red':
        # ~ for j in range(i,len(namesm)):
            # ~ nb = namesm[j]
            # ~ if nb[-7:-3] == 'blue' and nb[-2:]==n[-2:]:
                # ~ idcs.append([i,j])
                # ~ break
# ~ idcs = array(idcs)
# ~ save("sptrack/indices.npy",idcs)
     
    
