from numpy import *
from matplotlib.pylab import *
import scipy.optimize as opt

def fitGaus(roi):
    #print(roi[0])
    ar = int(roi[-2])
    mmax = roi[ar]
    mmin = roi[int(roi[-1])]

    arx = float(ar%5)
    ary = float(ar//5)
    #      Amplitude, sx, sy, angle, bg, x, y
    seed2 = (mmax-mmin,1.5,1.5,0.1,mmin,arx,ary)                 
    y = x = reshape(arange(0,25)%5,(5,5))
    try:
        popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), roi[:25], p0=seed2, \
        bounds=(array([0,0.3,0.3,0,0,0,0]),array([10000/5e-3,2,2,2*pi,5000/5e-3,5,5])))
        popt = concatenate((popt,[0]))

    except:
        popt = concatenate((seed,[1]))    
    return(popt) 
    
# popts = list(map(fitGaus,transpose(testm)))

def bruteforce(testm):
    popts = zeros((8,12000))
    arga = argmax(testm,axis=0)
    argi = argmin(testm,axis=0)
    for j,k in enumerate(arange(testm.shape[1])[:12000]):
        temp = testm[:,k]
        try:
            mmax = temp[arga[k]]
            mmin = temp[argi[k]]

            ar = arga[k]
            arx = float(ar%5)
            ary = float(ar//5)
            #      Amplitude, sx, sy, angle, bg, x, y
            seed2 = (mmax-mmin,1.5,1.5,0.1,mmin,arx,ary)                 
            popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), temp, p0=seed2, \
            bounds=(array([0,0.3,0.3,0,0,0,0]),array([10000/5e-3,2,2,2*pi,5000/5e-3,5,5])))
            popts[:7,j] = popt
        except:
            popts[:7,j] = seed2
            popts[-1,j] = 1
    return(popts)
                               
def twoD_Gaussian(xdata, amplitude = 1.0, sx = 1.0, sy= 1.0, theta = 0.0, offset= 0.0, x0 = 0.0, y0 = 0.0, ):
    x, y = xdata
    x0 = x0
    y0 = y0    
    a = (cos(theta)**2)/(2*sx**2) + (sin(theta)**2)/(2*sy**2)
    b = -(sin(2*theta))/(4*sx**2) + (sin(2*theta))/(4*sy**2)
    c = (sin(theta)**2)/(2*sx**2) + (cos(theta)**2)/(2*sy**2)
    g = offset + amplitude*exp( - (a*((x-x0)**2) + 2*b*(x-x0)*(y-y0) 
                            + c*((y-y0)**2)))
    return g.ravel()
    
