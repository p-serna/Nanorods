from numpy import *
from matplotlib.pylab import *
import os
#from core_libs import *
import scipy.optimize as opt
from scipy.optimize import minimize
import pickle

def gaussd(x,par): return(exp(-(x-par[0])**2/2/par[1]**2)/sqrt(2*pi*par[1]**2) )
def dblgaussd(x,par): return(par[0]*gaussd(x,par[1:3])+(1.0-par[0])*gaussd(x,par[3:]))
def minf(par,x,y,wy): return( sum(wy*(dblgaussd(x,par)-y)**2)/sum(wy))

def dblgausfit(x,y,wy=array([1.0]),par0=array([0.5,-1,.5,1,.5])):
    minx  = minimize(minf,par0,args=(x,y,wy))
    return minx

def msdestimate0(name,th = 0, th2 = 15.0):
    popts = load(name)
    posx = popts[:,5]
    posy = popts[:,6]
    amp = popts[:,0]
    sel = (popts[:,-1]==0)*(posx>0.1)*(posy>0.1)*(posx<4.9)*(posy<4.9)*(amp>exp(th))*(amp<exp(th2))
    Tf = len(posx)
    
    if sel.sum() > 1000:
        ts = arange(Tf)[sel]
        Tfn = ts[-1]
        #posx = posx[sel]
        #posy = posy[sel]
        lp = sum(sel)

        tmax = popts.shape[0]//2
        msd = zeros(tmax)
        msd2 = zeros(tmax)
        cnt = zeros(tmax)
        for i,t in enumerate(ts[:-5]):
            #print(i,t)
            i2 = min(i+tmax,lp)
            sel = ts[i:i2]
            t2 = ts[i:i2]-t
            sel = sel[t2<tmax]
            t2 = t2[t2< tmax]
            
            xs = posx[sel]
            ys = posy[sel]
            msd[t2] += (xs-xs[0])**2+(ys-ys[0])**2
            msd2[t2] += ((xs-xs[0])**2+(ys-ys[0])**2)**2
            cnt[t2] += 1

        #t = arange(tmax)

        msd01 = array(msd/cnt)
        msd01b = array(msd2/cnt)
        cnt01 = array(cnt)
        return(msd01,msd01b,cnt01)

def msdestimate(name,th = 0, th2 = 15.0,twoG = False):
    popts = load(name)
    posx = popts[:,5]
    posy = popts[:,6]
    amp = popts[:,0]
    sel = (popts[:,-1]==0)*(posx>0.1)*(posy>0.1)*(posx<4.9)*(posy<4.9)*(amp>exp(th))*(amp<exp(th2))
    Tf = len(posx)
    
    if sel.sum() > 1000:
        if twoG:
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
            if dgt[5]<dgt[3]:
                temp= ()
            
            sel2 = sel*(amb>min(dgt[5]+2*dgt[6],dgt[3]+2*dgt[4]))
            if sel2.sum()<500:
                sel2 = sel*(amb>max(dgt[5]-2*dgt[6],dgt[3]-2*dgt[4]))
            if sel2.sum()>500:
                sel = sel2
        
        #ts is the index of posx, posy
        ts = arange(Tf)[sel]
        Tfn = ts[-1]
        lp = sum(sel)
        

        tmax = popts.shape[0]//2
        msd = zeros(tmax)
        msd2 = zeros(tmax)
        cnt = zeros(tmax)
        for i,t in enumerate(ts[:-5]):
            #print(i,t)
            # For each i, and and t (= i or i+ #dropped frames)
            # we get i2 = i + maximum time used to estimate MSD 
            # clipped to the array length
            i2 = min(i+tmax,lp)
            # Now we get the indices for them
            sel = ts[i:i2]
            # Difference in time for those indices
            t2 = ts[i:i2]-t
            # We keep only indices for those s
            sel = sel[t2<tmax]
            t2 = t2[t2< tmax]
            
            xs = posx[sel]
            ys = posy[sel]
            msd[t2] += (xs-xs[0])**2+(ys-ys[0])**2
            msd2[t2] += ((xs-xs[0])**2+(ys-ys[0])**2)**2
            cnt[t2] += 1

        #t = arange(tmax)

        msd01 = array(msd/cnt)
        msd01b = array(msd2/cnt)
        cnt01 = array(cnt)
        return(msd01,msd01b,cnt01)


def msdanalysis(dfiles):
    data = []
    for ni,name in enumerate(dfiles):
        msd01,msd01b,cnt01 = msdestimate(name,th = 0, th2 = 15.0,twoG = True)
        data.append([ni,msd01,msd01b,cnt01])
        #print("File ",ni," of ", len(dfiles))
    return(data)
