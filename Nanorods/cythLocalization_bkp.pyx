#  python setup.py build_ext --inplace


import numpy as np
cimport numpy as np
ctypedef np.float64_t dtype_t
from scipy.optimize import curve_fit, minimize
cimport cython

@cython.boundscheck(False)
@cython.wraparound(False)
def twoD_Gaussian(xdata, lamplitude = 1.0, sx2 = 1.0, sy2= 1.0, theta = 0.0, offset= 0.0, x0 = 0.0, y0 = 0.0, ):
    x, y = xdata
    x0 = x0
    y0 = y0    
    a = (np.cos(theta)**2)/(2*sx2) + (np.sin(theta)**2)/(2*sy2)
    b = -(np.sin(2*theta))/(4*sx2) + (np.sin(2*theta))/(4*sy2)
    c = (np.sin(theta)**2)/(2*sx2) + (np.cos(theta)**2)/(2*sy2)
    g = offset + np.exp(lamplitude - (a*((x-x0)**2) + 2*b*(x-x0)*(y-y0) 
                            + c*((y-y0)**2)))
    return g.ravel()

@cython.boundscheck(False)
@cython.wraparound(False)
def twoD_GaussianA(xdata, lamplitude = 1.0, sx2 = 1.0, sy2= 1.0, theta = 0.0, offset= 0.0, x0 = 0.0, y0 = 0.0, ):
    x, y = xdata
    cth = np.cos(theta)
    sth = np.sin(theta)
    s2th = 2.0*sth*cth
    a = (cth**2)/(2*sx2) + (sth**2)/(2*sy2)
    b = -(s2th)/(4*sx2) + (s2th)/(4*sy2)
    c = (sth**2)/(2*sx2) + (cth**2)/(2*sy2)
    g = np.exp(offset) + np.exp(lamplitude - (a*((x-x0)**2) + 2*b*(x-x0)*(y-y0) 
                            + c*((y-y0)**2)))
    return g

      
@cython.boundscheck(False)
@cython.wraparound(False)
def twoD_Gaussianfxthpi4(xdata, lamplitude = 1.0, sx2 = 1.0, sy2= 1.0,  offset= 0.0, x0 = 0.0, y0 = 0.0, ):
    x, y = xdata
    cth = np.sqrt(0.5)
    sth = np.sqrt(0.5)
    s2th = 2.0*sth*cth
    a = (cth**2)/(2*sx2) + (sth**2)/(2*sy2)
    b = -(s2th)/(4*sx2) + (s2th)/(4*sy2)
    c = (sth**2)/(2*sx2) + (cth**2)/(2*sy2)
    g = np.exp(offset) + np.exp(lamplitude - (a*((x-x0)**2) + 2*b*(x-x0)*(y-y0) 
                            + c*((y-y0)**2)))
    return g
    
@cython.boundscheck(False)
@cython.wraparound(False)
def twoD_Gaussianfxth(xdata, lamplitude = 1.0, sx2 = 1.0, sy2= 1.0,  offset= 0.0, x0 = 0.0, y0 = 0.0, ):
    x, y = xdata
    cth = 1.0
    sth = 0.0
    s2th = 2.0*sth*cth
    a = (cth**2)/(2*sx2) + (sth**2)/(2*sy2)
    b = -(s2th)/(4*sx2) + (s2th)/(4*sy2)
    c = (sth**2)/(2*sx2) + (cth**2)/(2*sy2)
    g = np.exp(offset) + np.exp(lamplitude - (a*((x-x0)**2) + 2*b*(x-x0)*(y-y0) 
                            + c*((y-y0)**2)))
    return g
    
@cython.boundscheck(False)
@cython.wraparound(False)
def twoD_Gaussianfxsh(xdata, lamplitude = 1.0,  offset= 0.0, x0 = 0.0, y0 = 0.0, sx2 = 1.0, sy2= 1.0):
    x, y = xdata
    cth = 1.0
    sth = 0.0
    s2th = 2.0*sth*cth
    a = (cth**2)/(2*sx2) + (sth**2)/(2*sy2)
    b = -(s2th)/(4*sx2) + (s2th)/(4*sy2)
    c = (sth**2)/(2*sx2) + (cth**2)/(2*sy2)
    g = np.exp(offset) + np.exp(lamplitude - (a*((x-x0)**2) + 2*b*(x-x0)*(y-y0) 
                            + c*((y-y0)**2)))
    return g
   
@cython.boundscheck(False)
@cython.wraparound(False)
def twoD_GaussianV(xdata, lamplitude = np.ones((12000,25)), sx2 = 1.0, sy2= 1.0, theta = 0.0, offset= 0.0, x0 = np.ones((12000,25)), y0 = np.ones((12000,25)) ):
    x, y = xdata
    cth = np.cos(theta)
    sth = np.sin(theta)
    s2th = 2.0*sth*cth
    a = (cth**2)/(2*sx2) + (sth**2)/(2*sy2)
    b = -(s2th)/(4*sx2) + (s2th)/(4*sy2)
    c = (sth**2)/(2*sx2) + (cth**2)/(2*sy2)
    x0p = np.tile(x0,(25,1)).flatten()
    y0p = np.tile(y0,(25,1)).flatten()
    g = np.exp(np.tile(offset,(25,1)).transpose().flatten()) + np.exp(np.tile(lamplitude,(25,1)).transpose().flatten() - (a*((x-x0p)**2) + 2*b*(x-x0p)*(y-y0p) 
                            + c*((y-y0p)**2)))
    return g.ravel()
    

@cython.boundscheck(False)
@cython.wraparound(False)
def minimizefV(xp,roi=np.ones(25),x=np.array([np.arange(5) for j in range(5)]),y = np.transpose(np.array([np.arange(5) for j in range(5)]))):
    xp = np.array(xp)
    yt = twoD_GaussianV((x,y),xp[0:12000],xp[48000],xp[48001],xp[48002],xp[12000:24000],xp[24000:36000],xp[36000:48000])
    yt =  (yt-roi)/np.sqrt(roi)
    return(np.dot(yt,yt))

#~ popt = minimize(minimizef, args = (testm2p.flatten(),x, y), x0=(mmax,mmin,arx,ary,1.5,1.5,0.1))

@cython.boundscheck(False)
@cython.wraparound(False)
def fitGausV(testm2):
    
    ar = np.array(testm2[:,25],dtype=int)

    mmax = testm2[np.arange(12000),ar]
    mmin = testm2[np.arange(12000),np.array(testm2[:,26],dtype=int)]

    #The amplitude!
    mmax = np.log(mmax-mmin)
    mmin = np.log(mmin)
    arx =ar%5
    ary = ar//5
    #      Amplitude, sx, sy, angle, bg, x, y
    #seed2 = (mmax-mmin,1.5,1.5,0.1,mmin,arx,ary)                 
    x =np.array([np.arange(5) for j in range(5)])
    y =  np.transpose(x).flatten()
    y =  np.transpose(x).flatten()
    x = np.tile(x.flatten(),(12000))
    y = np.tile(y,(12000))
        
    try:
        #print("Hello")
        popt =  minimize(minimizefV, args = (testm2[:,:25].flatten(),x, y), x0=np.array(tuple(np.concatenate((mmax,mmin,arx,ary,[1.5],[1.5],[0.1])))))
        #print("Hello2")
        popt = np.concatenate((popt.x,[0]))
    except Exception as e: 
        #print(e)
        popt = np.concatenate((np.array([np.log(mmax-mmin),1.5,1.5,0.1,np.log(mmin),arx,ary]),[1]))    
    return(popt) 
    
    
 
@cython.boundscheck(False)
@cython.wraparound(False)
def TwodGJacF(xdata, lamplitude = 1.0, sx2 = 1.0, sy2= 1.0, theta = 0.0, offset= 0.0, x0 = 0.0, y0 = 0.0, ):
    x, y = xdata
    x = 1.0*x-x0
    y = 1.0*y-y0
    cth = np.cos(theta)
    sth = np.sin(theta)
    cth2 = cth*cth
    sth2 = sth*sth
    #s2th = np.sin(2*theta)
    s2th = 2*sth*cth
    c2th = cth2-sth2
    a0 = cth2
    b0 = sth2
    c0 = (sx2-sy2)*s2th
    jac = np.zeros((25,7))
    a = cth2/2.0/sx2 + sth2/2.0/sy2
    b = s2th*(-1.0/(sx2) + 1.0/(sy2))/2.0
    c = sth2/(2.0*sx2) + cth2/(2.0*sy2)
    #First one is the amplitude
    jac[:,0] = np.exp(lamplitude - (a*(x**2) + b*x*y+ c*(y**2)))
    jac[:,1] = jac[:,0]*(x*cth-y*sth)**2/2.0/sx2**2
    jac[:,2] = jac[:,0]*(y*cth+x*sth)**2/2.0/sy2**2
    jac[:,3] = jac[:,0]*(x*y*c2th+(x-y)*(x+y)*cth*sth)*(sy2-sx2)/sx2/sy2
    jac[:,4] = np.exp(offset)
    jac[:,5] = (a0*x/sx2+b0*x/sy2+c0*y/sx2/sy2/2.0)*jac[:,0]
    jac[:,6] = (a0*y/sy2+b0*y/sx2+c0*x/sx2/sy2/2.0)*jac[:,0]

    return(jac)

@cython.boundscheck(False)
@cython.wraparound(False)
def minimizef(xp,roi=np.ones(25),x=np.array([np.arange(5) for j in range(5)]),y = np.transpose(np.array([np.arange(5) for j in range(5)]))):
    yt = twoD_GaussianA((x,y),xp[0],xp[1],xp[2],xp[3],xp[4],xp[5],xp[6])
    yt =  (yt-roi)/np.sqrt(roi)
    return(np.dot(yt,yt))
    
@cython.boundscheck(False)
@cython.wraparound(False)
def jacfun(xp,roi=np.ones(25),x=np.array([np.arange(5) for j in range(5)]),y = np.transpose(np.array([np.arange(5) for j in range(5)]))):
    yt = twoD_GaussianA((x,y),xp[0],xp[1],xp[2],xp[3],xp[4],xp[5],xp[6])
    yt2 =  2*(yt-roi)/roi
    jac = TwodGJacF((x,y),xp[0],xp[1],xp[2],xp[3],xp[4],xp[5],xp[6])
    #jac = np.concatenate((np.dot(yt2,jac),[0]))
    jac = np.dot(yt2,jac)
    return(jac)
    
    
@cython.boundscheck(False)
@cython.wraparound(False)
def minimizegJac(roi):
    ar = int(roi[25])
    mmax = roi[ar]
    mmin = roi[int(roi[26])]

    arx = float(ar%5)
    ary = float(ar//5)
    #      Amplitude, sx, sy, angle, bg, x, y
    #seed2 = (mmax-mmin,1.5,1.5,0.1,mmin,arx,ary)                 
    x =np.array([np.arange(5) for j in range(5)])
    y =  np.transpose(x).flatten()
    x = x.flatten()
    try:
        minx = minimize(minimizef,x0=(np.log(mmax-mmin),1.5,1.5,0.1,np.log(mmin),arx,ary),args=(roi[:25],x,y),jac=jacfun)
        #print("Hello2")
        popt = np.concatenate((minx.x,[0]))
    except Exception as e: 
        #print(e)
        popt = np.concatenate((np.array([np.log(mmax-mmin),1.5,1.5,0.1,np.log(mmin),arx,ary]),[1]))    
    return(popt) 

@cython.boundscheck(False)
@cython.wraparound(False)
def minimizeg(roi):
    ar = int(roi[25])
    mmax = roi[ar]
    mmin = roi[int(roi[26])]

    arx = float(ar%5)
    ary = float(ar//5)
    #      Amplitude, sx, sy, angle, bg, x, y
    #seed2 = (mmax-mmin,1.5,1.5,0.1,mmin,arx,ary)                 
    x =np.array([np.arange(5) for j in range(5)])
    y =  np.transpose(x).flatten()
    x = x.flatten()
    try:
        minx = minimize(minimizef,x0=(np.log(mmax-mmin),1.5,1.5,0.1,np.log(mmin),arx,ary),args=(roi[:25],x,y))
        #print("Hello2")
        popt = np.concatenate((minx.x,[0]))
    except Exception as e: 
        #print(e)
        popt = np.concatenate((np.array([np.log(mmax-mmin),1.5,1.5,0.1,np.log(mmin),arx,ary]),[1]))    
    return(popt) 
    
@cython.boundscheck(False)
@cython.wraparound(False)
def fitGaus(roi):
    ar = int(roi[25])
    mmax = roi[ar]
    mmin = roi[int(roi[26])]

    arx = float(ar%5)
    ary = float(ar//5)
    #      Amplitude, sx, sy, angle, bg, x, y
    #seed2 = (mmax-mmin,1.5,1.5,0.1,mmin,arx,ary)                 
    x =np.array([np.arange(5) for j in range(5)])
    y =  np.transpose(x).flatten()
    x = x.flatten()
    try:
        #print("Hello")
        popt, pcov = curve_fit(twoD_GaussianA, (x, y), roi[:25], p0=(np.log(mmax-mmin),1.5,1.5,0.1,np.log(mmin),arx,ary))
        #print("Hello2")
        popt = np.concatenate((popt,[0]))
    except Exception as e: 
        #print(e)
        popt = np.concatenate((np.array([np.log(mmax-mmin),1.5,1.5,0.1,np.log(mmin),arx,ary]),[1]))    
    return(popt) 
    
@cython.boundscheck(False)
@cython.wraparound(False)
def fitGausfxth(roi):
    ar = int(roi[25])
    mmax = roi[ar]
    mmin = roi[int(roi[26])]

    arx = float(ar%5)
    ary = float(ar//5)
    #      Amplitude, sx, sy, bg, x, y
    #seed2 = (mmax-mmin,1.5,1.5,mmin,arx,ary)                 
    x =np.array([np.arange(5) for j in range(5)])
    y =  np.transpose(x).flatten()
    x = x.flatten()
    try:
        #print("Hello")
        popt, pcov = curve_fit(twoD_Gaussianfxth, (x, y), roi[:25], p0=(np.log(mmax-mmin),1.5,1.5,np.log(mmin),arx,ary))
        #print("Hello2")
        popt = np.concatenate((popt,[0]))
    except Exception as e: 
        #print(e)
        popt = np.concatenate((np.array([np.log(mmax-mmin),1.5,1.5,np.log(mmin),arx,ary]),[1]))    
    return(popt) 
    
@cython.boundscheck(False)
@cython.wraparound(False)
def fitGausfxthpi4(roi):
    ar = int(roi[25])
    mmax = roi[ar]
    mmin = roi[int(roi[26])]

    arx = float(ar%5)
    ary = float(ar//5)
    #      Amplitude, sx, sy, bg, x, y
    #seed2 = (mmax-mmin,1.5,1.5,mmin,arx,ary)                 
    x =np.array([np.arange(5) for j in range(5)])
    y =  np.transpose(x).flatten()
    x = x.flatten()
    try:
        #print("Hello")
        popt, pcov = curve_fit(twoD_Gaussianfxthpi4, (x, y), roi[:25], p0=(np.log(mmax-mmin),1.5,1.5,np.log(mmin),arx,ary))
        #print("Hello2")
        popt = np.concatenate((popt,[0]))
    except Exception as e: 
        #print(e)
        popt = np.concatenate((np.array([np.log(mmax-mmin),1.5,1.5,np.log(mmin),arx,ary]),[1]))    
    return(popt) 

   
@cython.boundscheck(False)
@cython.wraparound(False)
def fitGausfxsh(roi):
    s1, s2 = (roi[27],roi[28])
    ar = int(roi[25])
    mmax = roi[ar]
    mmin = roi[int(roi[26])]

    arx = float(ar%5)
    ary = float(ar//5)
    #      Amplitude, bg, x, y, sx, sy
    #seed2 = (mmax-mmin,mmin,arx,ary,s1,s2)                 
    x =np.array([np.arange(5) for j in range(5)])
    y =  np.transpose(x).flatten()
    x = x.flatten()
    try:
        #print("Hello")
        def gausfixsh(x,a,o,x0,y0): return(twoD_Gaussianfxsh(x, a,  o, x0, y0, s1, s2))
        
        popt, pcov = curve_fit(gausfixsh, (x, y), roi[:25], p0=(np.log(mmax-mmin),np.log(mmin),arx,ary))
        #print("Hello2")
        popt = np.concatenate((popt,[s1,s2,0]))
    except Exception as e: 
        #print(e)
        popt = np.concatenate((np.array([np.log(mmax-mmin),np.log(mmin),arx,ary,s1,s2]),[1]))    
    return(popt) 

@cython.boundscheck(False)
@cython.wraparound(False)
def fitGausJac(roi):
    ar = int(roi[25])
    mmax = roi[ar]
    mmin = roi[int(roi[26])]

    arx = float(ar%5)
    ary = float(ar//5)
    #      Amplitude, sx, sy, angle, bg, x, y
    #seed2 = (mmax-mmin,1.5,1.5,0.1,mmin,arx,ary)                 
    x =np.array([np.arange(5) for j in range(5)])
    y =  np.transpose(x).flatten()
    x = x.flatten()
    try:
        #print("Hello")
        popt, pcov = curve_fit(twoD_GaussianA, (x, y), roi[:25], p0=(np.log(mmax-mmin),1.5,1.5,0.1,np.log(mmin),arx,ary),jac=TwodGJacF)
        #print("Hello2")
        popt = np.concatenate((popt,[0]))
    except Exception as e: 
        #print(e)
        popt = np.concatenate((np.array([np.log(mmax-mmin),1.5,1.5,0.1,np.log(mmin),arx,ary]),[1]))    
    return(popt) 



import multiprocessing as mp

@cython.boundscheck(False)
@cython.wraparound(False)
def fitvideo(video):
    return(list(map(fitGaus,video)))
    

def fitvideoP(video,nthreads=32):
    pool = mp.Pool(nthreads)
    result = pool.map(fitGaus, video)
    pool.close()
    pool.join()
    result = np.array(result)
    result[:,0] = np.exp(result[:,0])
    result[:,1:3] = np.sqrt(result[:,1:3])
    return(result)

def fitvideoPfxth(video,nthreads=32):
    pool = mp.Pool(nthreads)
    result = pool.map(fitGausfxth, video)
    pool.close()
    pool.join()
    result = np.array(result)
    result[:,0] = np.exp(result[:,0])
    result[:,1:3] = np.sqrt(result[:,1:3])
    return(result)
    
def fitvideoPfxthpi4(video,nthreads=32):
    pool = mp.Pool(nthreads)
    result = pool.map(fitGausfxthpi4, video)
    pool.close()
    pool.join()
    result = np.array(result)
    result[:,0] = np.exp(result[:,0])
    result[:,1:3] = np.sqrt(result[:,1:3])
    return(result)
    
def fitvideoPfxsh(video,nthreads=32):
    pool = mp.Pool(nthreads)
    result = pool.map(fitGausfxsh, video)
    pool.close()
    pool.join()
    result = np.array(result)
    result[:,0:2] = np.exp(result[:,0:2])
    result[:,4:6] = np.sqrt(result[:,4:6])
    return(result)
