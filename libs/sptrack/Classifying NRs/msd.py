# ~ s1 = column_stack((dgps[:,6],dgps[:,4])).max(axis=1)
# ~ dG = (dgps[:,3]-dgps[:,5])/s1

# ~ dG[dG>20] = 0
from numpy import *
from scipy.stats import linregress


# ~ cmapt = get_cmap('tab20')
def linregress1(x,y):
    xm = x.mean()
    ym = y.mean()
    x2m = (x**2).sum()
    y2m = (y**2).sum()
    xym = (x*y).sum()
    n = x.shape[0]
    ssx = x2m-n*xm**2
    #ssy = y2m-n*ym**2
    ssxy = xym-n*xm*ym
    b = ssxy/ssx
    a = ym-b*xm
    return((a,b))
    
def linregress2(x,y):
    xm = x.mean()
    ym = y.mean()
    x2m = (x**2).sum()
    y2m = (y**2).sum()
    xym = (x*y).sum()
    n = x.shape[0]
    ssx = x2m-n*xm**2
    ssy = y2m-n*ym**2
    ssxy = xym-n*xm*ym
    b = ssxy/ssx
    seb = sqrt((ssx*ssy-ssxy**2)/((n-2)))/ssx
    return((b,seb))
    

# ~ cmapt = get_cmap('tab20')

def extractmsd(dataL,verbose=True):
    ''' Add a proper description here: it outputs averaged msd, and slopes, etc
    
    '''
    iNR = 0

    for idat,data in enumerate(dataL):
        tmax = data[0][1].shape[0]

        t = arange(tmax)
        counter = zeros(8)
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
                    dd[il,:2] = [xe.mean(),ye.mean()]
                t0l = til

            dde = dd[dd[:,0]>0,:]
            dde = dde[1:,:]

            if idat == 0 and i ==0:
                xde = column_stack((dde,dde[:,:1]*0+i,dde[:,:1]*0+idat))
                xden = column_stack((dde[:,:1],(dde[:,1:2]-dde[0,1])/(dde[-1,1]-dde[0,1]),dde[:,:1]*0+i,dde[:,:1]*0+idat))
            else:
                xde = row_stack((xde,column_stack((dde,dde[:,:1]*0+i,dde[:,:1]*0+idat))))
                xden = row_stack((xden,column_stack((dde[:,:1],(dde[:,1:2]-dde[0,1])/(dde[-1,1]-dde[0,1]),dde[:,:1]*0+i,dde[:,:1]*0+idat))))

            ntle = dde.shape[0]
            nwin = 50
            ds = zeros((ntle+nwin-3,3))
            for il in range(ntle+nwin-3):
                sel = arange(il-nwin+3,il+3)
                sel = sel[(sel>=0)*(sel<ntle)]
                xe = dde[sel,0]
                ye = dde[sel,1]    
                lm = linregress2(xe,ye)
                ds[il,0] = xe.mean()    
                ds[il,1:3] = (lm[0],lm[1])    


            if idat == 0 and i ==0:
                xds = column_stack((ds,ds[:,:1]*0+i,ds[:,:1]*0+idat))
            else:
                xds = row_stack((xds,column_stack((ds,ds[:,:1]*0+i,ds[:,:1]*0+idat))))
            
            iNR +=1
        if verbose:
            print(idat)

    xds = array(xds)
    xde = array(xde)
    xden = array(xden)
    return((xds,xde,xden))

def extractDe(dataL,xde,xds,verbose=True, ts = array([100,200,400,
  800,1600,3200,6400,12800])/1000.0):
    for idat,data in enumerate(dataL):
        xt = xde[(abs(xde[:,-2]-0)<1e-3)*(abs(xde[:,-1]-idat)<1e-3),0]
        idxts = zeros(ts.shape[0],dtype=int)
        idx = arange(xt.shape[0])
        for i in range(len(ts)):
            idxts[i] = idx[argmin(abs(xt-ts[i]))]

        alphas = zeros((len(data),7))
        alphalm = zeros((len(data),2))

        Det0 = zeros((len(data),3))
        Det1 = zeros((len(data),6))

        for i in range(len(data)):
            xt = xds[(abs(xds[:,-2]-i)<1e-3)*(abs(xds[:,-1]-idat)<1e-3),:]
            xt = xt[xt[:,1]>0,:]
            try:
                De = array([mean(xt[xt[:,0]<2e-1,1]),mean(xt[xt[:,0]<2e-1,2]),max(xt[xt[:,0]<2e-1,2])])
                sel = (xt[:,0]>2e-1)*(xt[:,0]<1.5e0)
                De2 = array([mean(xt[sel,1]),mean(xt[sel,2]),max(xt[sel,2])])
                sel = (xt[:,0]>1.5e0)*(xt[:,0]<1.5e1)
                De3 = array([mean(xt[sel,1]),mean(xt[sel,2]),max(xt[sel,2])])

                xt = xde[(abs(xde[:,-2]-i)<1e-3)*(abs(xde[:,-1]-idat)<1e-3),:]    
                msds = xt[idxts,:]
                msds = msds[1:,:]/msds[:-1,:]
                alpha = log(msds[:,1])/log(msds[:,0])
                lm = linregress1(log(msds[:,0]).cumsum(),alpha)
                alphas[i,:] = alpha
                alphalm[i,:] = (lm[1],lm[0])
                Det0[i,:] = De
                Det1[i,:] = concatenate((De2,De3))
            except:
                pass
        if idat ==0:
            Des = column_stack((Det0,Det0[:,0]*0+idat))
            Dep = column_stack((Det1,Det1[:,0]*0+idat))
        else:
            Des = row_stack((Des,column_stack((Det0,Det0[:,0]*0+idat))))
            Dep = row_stack((Dep,column_stack((Det1,Det1[:,0]*0+idat))))
        if verbose:
            print(idat)
    return((Des,Dep))


def main():
    for idat,data in enumerate(dataL):
        tmax = data[0][1].shape[0]

        t = arange(tmax)


        for ic in range(7):
            figure(ic)
            title(conditionlabel(ic))


        counter = zeros(8)
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

            if idat == 0 and i ==0:
                xde = column_stack((dde,dde[:,:1]*0+i,dde[:,:1]*0+idat))
                xden = column_stack((dde[:,:1],(dde[:,1:2]-dde[0,1])/(dde[-1,1]-dde[0,1]),dde[:,:1]*0+i,dde[:,:1]*0+idat))
            else:
                xde = row_stack((xde,column_stack((dde,dde[:,:1]*0+i,dde[:,:1]*0+idat))))
                xden = row_stack((xden,column_stack((dde[:,:1],(dde[:,1:2]-dde[0,1])/(dde[-1,1]-dde[0,1]),dde[:,:1]*0+i,dde[:,:1]*0+idat))))

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

            #print(i)

            xt = ds
            xt = xt[xt[:,1]>0,:]
            De = array([mean(xt[xt[:,0]<1e-1,1]),exp(mean(log(xt[xt[:,0]<1e-1,1])))])

            ic = conditionnumber(ampst[i,0]/ampst[i,1],De[0],dG[iNR])
            if counter[ic]< 20:
                figure(ic)
                plot(dde[:,0],dde[:,1]-dde[0,1],'.-',alpha=0.3,label=i)
                counter[ic] +=1

            if idat == 0 and i ==0:
                xds = column_stack((ds,ds[:,:1]*0+i,ds[:,:1]*0+idat))
            else:
                xds = row_stack((xds,column_stack((ds,ds[:,:1]*0+i,ds[:,:1]*0+idat))))
            
            iNR +=1
        print(idat)
