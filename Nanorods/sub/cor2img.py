from numpy import *
from matplotlib.pylab import *
from scipy.optimize import minimize

# =============================================
#  Functions to calculate minimum distances
#
def mindist(p0,vecpos):
    '''
    - mindist calculate the minimum squared distance of a vector p0 to a set of
    vectors vecpos.
    - it returns minimum distance and the position of the vector that minimize it.
    '''
    dst = ((vecpos-p0)*(vecpos-p0)).sum(axis=1)
    return(dst.min(),dst.argmin())

def vmindist(vecA,vecB,arg = False):
    '''
    cordist calculate vector of minima squared distances of vector from vector set vecA to vectir set vecB.
    It returns ta vector of the (minimal) squared distances divided by the number of vectors in vecA.
    '''
    dA = vecA.shape
    dB = vecB.shape
    tt1 = tile(vecA.transpose(),(dB[0],1,1))
    ttn = tt1.transpose()-vecB.transpose()
    ttn = (ttn**2).sum(axis=1)
    tmin = ttn.min(axis=1)
    if arg:
        targ = ttn.argmin(axis=1)
        return(column_stack((tmin,targ)))
    else:
        return(tmin)
    
    
def cordist(vecA,vecB):
    '''
    cordist calculate the sum of the minima squared distances of vector from vector set vecA to vectir set vecB.
    It returns the sum of the (minimal) squared distances divided by the number of vectors in vecA.
    '''
    dA = vecA.shape
    dB = vecB.shape
    tt1 = tile(vecA.transpose(),(dB[0],1,1))
    ttn = tt1.transpose()-vecB.transpose()
    tmin = (ttn**2).sum(axis=1).min(axis=1)
    #dst = ((vecpos-p0)*(vecpos-p0)).sum(axis=1)
    return(tmin.sum()/tmin.shape[0])


# =============================================
#  Functions to translate, translate and multiply, and translate, multiply and rotate.
#
def mintras(v,vecA,vecB):
    vecAt = vecA+v
    return(cordist(vecAt,vecB))

def minlr(v,vecA,vecB):
    v0 = v[:2]
    v1 = v[2:4]
    vecAt = vecA*v1+v0
    return(cordist(vecAt,vecB))

def minfr(v,vecA,vecB):
    v0 = v[:2]
    v1 = v[2:4]
    v2 = array([[0,v[4]],[v[5],0]])
    vecA2 = dot(v2,vecA.transpose()).transpose()
    vecAt = vecA*v1+vecA2+v0
    return(cordist(vecAt,vecB))

def minfrpr(v,vecA,vecB):
    v0 = v[:2]
    v1 = v[2:4]
    v2 = array([[0,v[4]],[v[5],0]])
    vecA2 = dot(v2,vecA.transpose()).transpose()
    vecAt = vecA*v1+vecA2+v0
    vmin = vmindist(vecAt,vecB,arg=False)

    return(sum(vmin<3**2)/vmin.shape[0])
    
# =============================================
#  Function that minimizes the chi^2 defined by cordist
#
def transfpar(vecA,vecB,transform = None,x0 = None):
    '''
    transfpar minimize the transformation to vecA to obtain minimal distances to vecB
     Inputs:
      vecA - vector with shape (nA,spatial dimensions) (tested with spatial dim =2)
      vecB - vector with shape (nB,spatial dimensions) (tested with spatial dim =2)
      transform - if default or < 1 it only minimize when transformation is a translation
                    ([x,y])
                - if 1, translation and dilation with a matrix [[ax,0],[0,ay]]
                - if 2 or more, matrix [[ax,bxy],[cxy,ay]]
     Output:
      array([x,y,ax,ay,bxy,cxy])
    '''
    if x0 is None:
        x0 = array([randn(),randn()])
    minxt = minimize(mintras,x0=x0,args=(vecA,vecB))
    xf = minxt.x
    if transform>0:
        minxl = minimize(minlr,x0=array([xf[0],xf[1],1.0,1.0]),args=(vecA,vecB))
        xf = minxl.x    
    if transform>1:
        minxf = minimize(minfr,x0=concatenate((xf,array([0.0,0.0]))),args=(vecA,vecB))
        xf = minxf.x
    return(xf)
    
# =============================================
#  Function that gives ROIs that appear in both
#

def coincidROI(vecA,vecB,err = 3):
    '''
    It returns positions that are coincident in vecA first column and vecB second col
    '''
    vmin = vmindist(vecA,vecB,arg = True)
    vAid = arange(vecA.shape[0])
    vAid = vAid[vmin[:,0]<err**2]
    selection = vmin[vmin[:,0]<err**2,1]
    return(column_stack((vAid,selection)) )

# =============================================
#  Main function showing an example. It needs to be provided fposA, fposB and an image
#  imB, and to import visualization function from subs.
#    
if __name__ == '__main__':
     
    v = transfpar(fposA,fposB,transform = 2)
    
    v0 = v[:2]
    v1 = v[2:4]
    v2 = array([[0,v[4]],[v[5],0]])
    fposAinB =  fposA*v1+v0+dot(v2,fposA.transpose()).transpose()
    visualization(imB,fposAinB,figname="Red channel 2",color='green')
