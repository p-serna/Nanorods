from numpy import *
from matplotlib.pylab import *
import pytiff
from scipy.optimize import minimize,least_squares
import os
import sys


def f0(x,par): return(exp(par[0])+par[1]*x+exp(par[2]-x/par[3]))



def funcfit(fun,x,y,par0,ey=1.0):
    def minf(par): return( sum((fun(x,par)-y)**2/ey**2))
    minx  = minimize(minf,par0)
    return minx


def readBigtifFile(fname,dimensions = None):
    with pytiff.Tiff(fname) as handle:
        tags = handle.read_tags()
        offset = tags['strip_offsets'][0]
        nframes = int(tags['image_description'].split()[1][7:])
        width, height = (tags['image_width'][0],tags['image_length'][0])
    if dimensions is None:
        dimensions = (nframes,height,width)
    with open(fname,"rb") as file:
        file.seek(offset)
        temp = fromfile(file,dtype=">u2")
    #print(temp.shape)
    
    return(temp.reshape(dimensions))

f = sys.argv[1]

movie = readBigtifFile(f)
nframes, height, width = movie.shape
im1 = movie[0,:,:]/900
im1[im1>1] = 1
xmax = log10(max(movie[array([0,1,-1]),:,:].flatten()))
xmin = log10(min(movie[array([0,1,-1]),:,:].flatten()))
dx = (xmax-xmin)/254
xbreaks = arange(xmin-dx/2.0,xmax+dx/2.0,dx)
xmid = 10**(.5*(xbreaks[1:]+xbreaks[:-1]))
h0 = histogram(log10(movie[0,:,:].flatten()),bins=xbreaks)[0]
h1 = histogram(log10(movie[1,:,:].flatten()),bins=xbreaks)[0]
hf = histogram(log10(movie[-1,:,:].flatten()),bins=xbreaks)[0]
figure("I distribution")
plot(xmid,h0*10**dx/sum(h0),'.-',label="t=0")
plot(xmid,h1*10**dx/sum(h1),'.-',label="t=1")
plot(xmid,hf*10**dx/sum(hf),'.-',label="t=12000")
yscale("log")
xscale("log")
legend()
xlabel("Intensity (counts)")
ylabel("P(I)")

savefig("Idist_"+f[2:-4]+".png")
close("I distribution")

#imshow(im1)

#totalF = sum(movie.reshape(nframes,height*width),axis=-1)
ts = arange(nframes)

width2 = width//2
totalFA = array(sum(sum(movie[:,:,:width2],axis=-1),axis=-1),dtype=float)
totalFB = array(sum(sum(movie[:,:,width2:],axis=-1),axis=-1),dtype=float)


cf0 = (totalFA[-1]-totalFA[4000])/8000.0
cf1 = totalFA[4000]-cf0*4000.0
par0 = array([log(cf1),cf0,log(totalFA[1]-cf1),3000.0])
fitA = funcfit(f0,array(ts[1:],dtype=float),totalFA[1:],par0,ey = totalFA[1:]*0.005)
cf0 = (totalFB[-1]-totalFB[4000])/8000.0
cf1 = totalFB[4000]-cf0*4000.0
par0 = array([log(cf1),cf0,log(totalFA[1]-cf1),3000.0])
fitB = funcfit(f0,array(ts[1:],dtype=float),totalFB[1:],par0,ey = totalFB[1:]*0.005)


figure("Bleaching")
plot(totalFA)
plot(ts,f0(ts,fitA.x),'k--')
plot(totalFB) 
plot(ts,f0(ts,fitB.x),'k--')
xlabel("frames")
ylabel("Total intensity")
savefig("TotI_"+f[2:-4]+".png")
close("Bleaching")

x = arange(width2)
y = arange(height)
XX, YY = meshgrid(x, y)

xs = zeros((nframes,4))
for i in range(nframes):
    im1 = movie[i,:,:]/900.0
    im1[im1>1] = 1
    xcm = sum((XX*im1[:,:width2]).flatten())/sum(im1[:,:width2])
    ycm = sum((YY*im1[:,:width2]).flatten())/sum(im1[:,:width2])

    xcm2 = sum((XX*im1[:,width2:]).flatten())/sum(im1[:,width2:])
    ycm2 = sum((YY*im1[:,width2:]).flatten())/sum(im1[:,width2:])
    xs[i,:] = array([xcm,ycm,xcm2,ycm2])

figure("Avg of X")
plot(xs[:,0])
plot(xs[:,2])
xlabel("frames")
ylabel("<x>")
xminA = min(xs[:,0])
xminB = min(xs[:,2])
xmaxA = max(xs[:,0])
xmaxB = max(xs[:,2])
dxA = xmaxA-xminA
dxB = xmaxB-xminB
if dxA>dxB:
    dx = dxA
    xmin = xminA
    xmax = xmaxA
else:
    dx = dxB
    xmin = xminB
    xmax = xmaxB
    
plot(ts,ts*0+xmin,'k:')
plot(ts,ts*0+xmax,'k:')
arrow(12200,xmin,0.0,xmax-xmin)
dx = (xmax-xmin)
text(12000,xmin+(xmax-xmin)*.2,'''maximum drift = {:.2f} pixel\n
                    ~{:.1f} nm'''.format(dx,325*dx), horizontalalignment='right')
    
savefig("AvgX_"+f[2:-4]+".png")
close("Avg of X")
dxX = dx
figure("Avg of Y")
plot(xs[:,1])
plot(xs[:,3])

xlabel("frames")
ylabel("<y>")
xminA = min(xs[:,1])
xminB = min(xs[:,3])
xmaxA = max(xs[:,1])
xmaxB = max(xs[:,3])
dxA = xmaxA-xminA
dxB = xmaxB-xminB
if dxA>dxB:
    dx = dxA
    xmin = xminA
    xmax = xmaxA
else:
    dx = dxB
    xmin = xminB
    xmax = xmaxB
    
plot(ts,ts*0+xmin,'k:')
plot(ts,ts*0+xmax,'k:')
arrow(12200,xmin,0.0,xmax-xmin)
dx = (xmax-xmin)
text(12000,xmin+(xmax-xmin)*.2,'''maximum drift = {:.2f} pixel\n
                    ~{:.1f} nm'''.format(dx,325*dx), horizontalalignment='right')
    
savefig("AvgY_"+f[2:-4]+".png")
close("Avg of Y")
dxY = dx

# ~ xs2 = zeros((12000,4))
# ~ for i in range(12000):
    # ~ im1 = movie[i,:,:]/890.0
    # ~ im1[im1>1] = 0
    # ~ xcm = sum((XX*im1[:,:512]).flatten())/sum(im1[:,:512])
    # ~ ycm = sum((YY*im1[:,:512]).flatten())/sum(im1[:,:512])

    # ~ xcm2 = sum((XX*im1[:,512:]).flatten())/sum(im1[:,512:])
    # ~ ycm2 = sum((YY*im1[:,512:]).flatten())/sum(im1[:,512:])
    # ~ xs2[i,:] = array([xcm,ycm,xcm2,ycm2])



figure(5)
xst = xs[:,0]-min(xs[:,0])
cf0 = (log(xst[1000])-log(xst[100]))/log(1000/100)
plot(xst,'.-')
plot(arange(10000),xst[1000]*arange(10000)**cf0/1e3**cf0)

xst = xs[:,2]-min(xs[:,2])
cf1 = (log(xst[1000])-log(xst[100]))/log(1000/100)
plot(xst,'.-')
plot(arange(10000),xst[1000]*arange(10000)**cf1/1e3**cf1)
xscale("log")
yscale("log")
xlabel("frames")
ylabel("<x>")
print(cf0,cf1)
savefig("AvXt_"+f[2:-4]+".png")
#cf = 0.5; plot(ts,ts**cf/1e3**cf*0.47912,'C1--')
close(5)
data.append([fitA.x,fitB.x,dxX,dxY,cf0,cf1])

    # ~ from scipy import signal

    # ~ im1 = movie[1,:,512:]/890.0
    # ~ im1[im1>1] = 0
    # ~ im2 = movie[-1,:,512:]/890.0
    # ~ im2[im2>1] = 0

    # ~ crscor = signal.correlate2d(im1,im2)

    # ~ csh = crscor.shape
    # ~ cmax = argmax(crscor.flatten())
    # ~ ix,iy = int(cmax%csh[1]),int(cmax/csh[1])
    # ~ offset = [ix - (csh[1]+1)/2,iy - (csh[0]+1)/2]

    # ~ figure(1)
    # ~ imshow(im1)
    # ~ figure(2)
    # ~ imshow(im2)
    # ~ figure(3)
    # ~ imshow(crscor)
