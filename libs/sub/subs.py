from numpy import *
from matplotlib.pylab import *
#import pytiff
from scipy.optimize import minimize,least_squares
import cv2
from codecs import open as copen

def f0(x,par): return(exp(par[0])+par[1]*x+exp(par[2]-x/par[3]))
 
def funcfit(fun,x,y,par0,ey=1.0):
    def minf(par): return( sum((fun(x,par)-y)**2/ey**2))
    minx  = minimize(minf,par0)
    return minx


# ~ cmask = pytiff.Tiff("cell3_Mask_CMOS.tif")       
# ~ cm = cmask[:512,:512]

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
    # print(temp.shape)
    
    return(temp.reshape(dimensions))

def readtifFile(fname,dimensions = None):
    with pytiff.Tiff(f) as handle:
        tags = handle.read_tags()
        offset = tags['strip_offsets'][0]
        nframes = int(tags['image_description'].split()[1][7:])
        width, height = (tags['image_width'][0],tags['image_length'][0])
        movie = zeros((nframes,height,width))
        for i,pages in enumerate(handle):
            movie[i,:,:] = pages
    
    return(movie)

def readtifImage(fname):    
    return(pytiff.Tiff(fname)[:,:])

def gettimes(fname,nframes=None):
    if nframes is None:
        try:
            with pytiff.Tiff(fname) as handle:
                tags = handle.read_tags()
                nframes = int(tags['image_description'].split()[2][7:])
        except:
            print("I could not get the number of frames, please provide it")
            return

    with copen(fname,"r","windows-1252") as f:
        j = 0
        times = zeros((nframes))
        while True:
            try:
                line = f.readline()
                # print(line)
                linesp = line.replace('\x00','').strip().split()
                if len(linesp)==5:
                    if linesp[2]=='Time_From_Last':
                        k,t = linesp[1],linesp[-1]
                        # print(int(k),float(t))
                        times[int(k)-1] = float(t)
                        j = j+1
            except:
                break
            if j>=nframes:
                break
    return(times)




def extractpeaks2d(img,kernel=7,width=None,height=None):
    if width is None or height is None:
        height, width = img.shape
    if type(kernel) is int:
        kernelm = ones((kernel,kernel),float32)/float(kernel)**2
    elif type(kernel) is ndarray:
        kernelm = kernel
    else:
        print("Error kernel")
        return()
    dilated = cv2.dilate(img,kernel=kernelm)

    pos = arange(width*height)[dilated.flatten() == img.flatten()]
    prspctpeaks7 = img.flatten()[pos]
    sel = (-prspctpeaks7).argsort()
    #                          x                y : imA[y,x]
    posA = column_stack((pos[sel]%width,pos[sel]//width))
    return(posA)

def extractpeaks2dTH(img,kernel=7,width=None,height=None):
    if width is None or height is None:
        height, width = img.shape
    if type(kernel) is int:
        kernelm = ones((kernel,kernel),float32)/49.0
    elif type(kernel) is numpy.ndarray:
        kernelm = kernel
    else:
        print("Error kernel")
        return()
    tophat = cv2.morphologyEx(img, cv2.MORPH_TOPHAT, kernel)
    dilated = cv2.dilate(tophat,kernel=kernelm)
    pos = arange(width*height)[dilated.flatten() == img.flatten()]
    prspctpeaks7 = img.flatten()[pos]
    sel = (-prspctpeaks7).argsort()
    #                          x                y : imA[y,x]
    posA = column_stack((pos[sel]%width,pos[sel]//width))
    return(posA)

def stripborder(pos,height = None,width = None,minx = 0,miny = 0):
    if height is None:
        height = max(pos[:,1])
    if width is None:
        width = max(pos[:,0])
    sel = (pos[:,0] <width-1)*(pos[:,0] >minx)*(pos[:,1] <height-1)*(pos[:,1] >miny)
    return(pos[sel,:])


def selROI(movie,x,y,t = 0,tw = None,xlim = (0,512),ylim= (0,256),size = 1):
    if x<xlim[0]+size:
        x0,xf = (xlim[0],xlim[0]+2*size+1)
    elif x>xlim[1]-size-1:
        x0,xf = (xlim[1]-2*size-1,xlim[1])
    else:
        x0,xf = (x-size,x+size+1)
    if y<ylim[0]+size:
        y0,yf = (ylim[0],ylim[0]+2*size+1)
    elif y>ylim[1]-size-1:
        y0,yf = (ylim[1]-2*size-1,ylim[1])
    else:
        y0,yf = (y-size,y+size+1)
    
    if tw is None:
        return(movie[t:,y0:yf,x0:xf])
    else:
        return(movie[t:(t+tw),y0:yf,x0:xf])       

    
def visualization(imA,pos=None,widthr=5,heightr=5,contrastd=1.0,contrastu=2.5,figname="",color='red',figsize = (7,4.5)):
    fig = figure(figname,figsize = figsize)
    ax = fig.add_subplot(111)
    imt = imA*1.0; me = mean(imt.flatten()); sd = std(imt.flatten());
    imt[imt>me+contrastu*sd] = me+contrastu*sd; imt[imt<me-contrastd*sd] = me-contrastd*sd;
    ax.imshow(imt,cmap='gray')
    if pos is None:
        pass
    else:
        for a_x, a_y in pos:
            ax.add_patch(Rectangle(xy=(a_x-widthr/2, a_y-heightr/2) ,width=widthr, height=heightr, linewidth=1.5, color=color,alpha=0.95, fill=False))
    axis('off')
    subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0)
    margins(0,0)
    gca().xaxis.set_major_locator(NullLocator())
    gca().yaxis.set_major_locator(NullLocator())

    return((fig,ax))

def projectmovie(cfile):
    # ----------------------------------------------------------------
    # We read the file.     

    movie = readBigtifFile(cfile)
    tduration, mheight,mwidth = movie.shape
   
    return(movie.mean(axis=0))

    
def running_mean(x, N):
    cumsumt = cumsum(insert(x, 0, 0)) 
    return (cumsumt[N:] - cumsumt[:-N]) / float(N)

