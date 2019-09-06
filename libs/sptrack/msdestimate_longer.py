from numpy import *
from matplotlib.pylab import *
import os
#from core_libs import *
import scipy.optimize as opt
from scipy.optimize import minimize
import pickle
import multiprocessing as mp

def gaussd(x,par): return(exp(-(x-par[0])**2/2/par[1]**2)/sqrt(2*pi*par[1]**2) )
def dblgaussd(x,par): return(par[0]*gaussd(x,par[1:3])+(1.0-par[0])*gaussd(x,par[3:]))

def dblgausfit(x,y,wy=1.0,par0=array([0.5,-1,.5,1,.5])):
    def minf(par): return( sum(wy*(dblgaussd(x,par)-y)**2)/sum(wy))
    minx  = minimize(minf,par0)
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

def msdestimate(name,th = 0, th2 = 15.0,twoG = False,minlen = 500):
    popts = load(name)
    posx = popts[:,5]
    posy = popts[:,6]
    amp = popts[:,0]
    sel = (popts[:,-1]==0)*(posx>0.1)*(posy>0.1)*(posx<4.9)*(posy<4.9)*(amp>exp(th))*(amp<exp(th2))
    Tf = len(posx)
    
    if sel.sum() > minlen:
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
            xs = xs-xs[0]
            ys = ys-ys[0]
            xs = xs*xs
            ys = ys*ys
            msd[t2] += xs+ys
            msd2[t2] += (xs+ys)*(xs+ys)
            cnt[t2] += 1

        #t = arange(tmax)

        msd01 = msd/cnt
        msd01b = msd2/cnt
        cnt01 = cnt
    else:
        #tmax = 3000
        msd01 = nan
        msd01b = nan
        cnt01 = nan
    return(msd01,msd01b,cnt01)
    
def msdanalysis(dfiles):
    data = []
    for ni,name in enumerate(dfiles):
        msd01,msd01b,cnt01 = msdestimate(name,th = 0, th2 = 15.0,twoG = True)
        data.append([ni,msd01,msd01b,cnt01])
        print("File ",ni," of ", len(dfiles))
    return(data)

def msdtemp(name):
    return(list(msdestimate(name,th = 0, th2 = 15.0,twoG = True)))    
    
def msdanalysisa(dfiles,nthreads = 32):
    pool = mp.Pool(nthreads)

    results = pool.map(msdtemp, dfiles)
    pool.close()
    pool.join()
    data = [ [i,r[0],r[1],r[2]] for i,r in enumerate(results)] 
    return(data)

def main():
    if len(sys.argv)>1:
        wdir = sys.argv[1]
        if wdir[-1] != '/':
            wdir = wdir+'/'
        if len(sys.argv)>2:
            driftcorrected = sys.argv[2]=='1'
        else:
            driftcorrected = True
    else:
        wdir = "./"
        driftcorrected = True
        print("No argument, we do it for current folder")


    #driftcorrected = True

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

    if len(dfiles)== 0:
        print("0 files found, do you want drift corrected? Did you run it?")
    else:
        print(len(dfiles)," files found, let's go!")
        
    dfiles.sort()

    print("In ",wdir)
    data = msdanalysisa(dfiles,nthreads = 32)
    print("Finished")
    # ~ data = []
    # ~ for ni,name in enumerate(dfiles):
        # ~ msd01,msd01b,cnt01 = msdestimate(name,th = 0, th2 = 15.0,twoG = True)
        # ~ data.append([ni,msd01,msd01b,cnt01])
        # ~ print("File ",ni," of ", len(dfiles))


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
        


        # ~ popts = load(name)
        # ~ posx = popts[:,5]
        # ~ posy = popts[:,6]
        # ~ amp = popts[:,0]
        # ~ th = 0
        # ~ th2 = 13.0
        # ~ sel = (popts[:,-1]==0)*(posx>0.1)*(posy>0.1)*(posx<4.9)*(posy<4.9)*(amp>exp(th))*(amp<exp(th2))
        # ~ Tf = len(posx)
        
        # ~ if sel.sum() > 1000:
            # ~ ts = arange(Tf)[sel]
            # ~ Tfn = ts[-1]
            # ~ #posx = posx[sel]
            # ~ #posy = posy[sel]
            # ~ lp = sum(sel)

            # ~ tmax = popts.shape[0]//2
            # ~ msd = zeros(tmax)
            # ~ msd2 = zeros(tmax)
            # ~ cnt = zeros(tmax)
            # ~ for i,t in enumerate(ts[:-5]):
                # ~ #print(i,t)
                # ~ i2 = min(i+tmax,lp)
                # ~ sel = ts[i:i2]
                # ~ t2 = ts[i:i2]-t
                # ~ sel = sel[t2<tmax]
                # ~ t2 = t2[t2< tmax]
                
                # ~ xs = posx[sel]
                # ~ ys = posy[sel]
                # ~ msd[t2] += (xs-xs[0])**2+(ys-ys[0])**2
                # ~ msd2[t2] += ((xs-xs[0])**2+(ys-ys[0])**2)**2
                # ~ cnt[t2] += 1

            # ~ #t = arange(tmax)

            # ~ msd01 = array(msd/cnt)
            # ~ msd01b = array(msd2/cnt)
            # ~ cnt01 = array(cnt)
            # ~ data.append([ni,msd01,msd01b,cnt01])
            
            # # ~ figure(2)

            # # ~ plot(t[2:],sqrt((msd01[2:]-msd01[1])/msd01[-1]),alpha=0.5)
            # # ~ figure(3)

            ## ~ plot(t[2:],(sqrt(msd01[2:])-sqrt(msd01[1]))/sqrt(msd01[-1]),alpha=0.5)
        

    if driftcorrected:
        with open(wdir+'data_msdB.pickle', 'wb') as handle:
            pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with open(wdir+'data_msd_U.pickle', 'wb') as handle:
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
         
if __name__ == "__main__":
    main()    
