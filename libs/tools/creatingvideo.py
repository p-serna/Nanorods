import sys
sys.path.append("/users/bssn/serna/GitIBENS/Nanorods/")
from ROIS_analysis_all_speeding import fitROI

def ROIanimation(movie):
    St = zeros((6000,2))
    St[:,0] = arange(6000)*10e-3

    for i in range(6000):
        imt = movie[i,:,:]
        St[i,1] = (imt-imt.min()).sum()
    St2 = column_stack(((St[:,0].reshape(1500,4))[:,:2].flatten(),(St[:,1].reshape(1500,4))[:,:2].flatten()))

    fig, axf = subplots(nrows=2,ncols=1)
    ax = axf[0]

    ax2 = axf[1]

    x = randn(1000)
    y = randn(1000)

    foox = []
    fooy = []


    ax.set_xlim(0,1)
    ax.set_ylim(min(St[:,1]),max(St[:,1]))

    imt = movie[0,:,:]
    im = ax2.imshow(imt,cmap='gray')
    im.set_clim(movie[1:,:,:].min(),movie[1:,:,:].max())

    #sc = ax0.scatter(0,(imt-imt.min()).sum())

    line, = ax.plot([], [], '.-', lw=2,alpha=0.6)
    line2, = ax.plot([], [], 'C1.', lw=2,alpha=0.6)

    foo, = ax2.plot([], [], '.', lw=2)

    time_template = 'time = %.1fs'
    time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

    def init():
        line.set_data([], [])
        line2.set_data([], [])
        im.set_data(movie[0,:,:])
        time_text.set_text('')
        ax.set_xlim(0,1)
        return line, time_text, im,ax,line2


    def animate(i):
        thisx = St[max(0,i-80):(i+1),0]
        thisy = St[max(0,i-80):(i+1),1]

        imt = movie[i,:,:]
        im.set_array(imt)

        line.set_data(thisx, thisy)
        line2.set_data(St2[max(0,(i-80)//2):(i//2+1),0], St2[max(0,(i-80)//2):(i//2+1),1])

        time_text.set_text(time_template%(i*dt))
        ax.set_xlim(max(0,i-80)*10e-3,max(i+20,100)*10e-3)
        return line, time_text, im,ax,line2


    ani = animation.FuncAnimation(fig, animate, np.arange(6000),
        interval=2, blit=False, init_func=init)
    return(St,ani,fig)


%pylab
def ftvaugmentr(img0,magn=2):
    img = img0.transpose()
    sh = img.shape
    sh2 = array([sh[0],sh[1]])*(magn-1)
    sh2 = ((sh2[0])//2+1,(sh2[0])//2+1)
    fftim  = fftshift(fft2(img,axes=(0,1)),axes=(0,1))
    fftim  =  pad(fftim,(sh2[0],sh2[1]),'constant') 
    imgn = real(ifft2(ifftshift(fftim,axes=(0,1)),axes=(0,1)))
    imgn = imgn[:,:,sh2[0]:-sh2[0]]
    return(imgn.transpose())
    
rois = load("roi_moviewave.npy")     
roirs = rois.reshape(rois.shape[0],6000,5,5)  
#roia = map(ftvaugment,roirs[:4,:,:,:])





from scipy.optimize import minimize
def gaussd(x,par): return(exp(-(x-par[0])**2/2/par[1]**2)/sqrt(2*pi*par[1]**2) )
def dblgaussd(x,par): return(par[0]*gaussd(x,par[1:3])+(1.0-par[0])*gaussd(x,par[3:]))
def dblgausfit(x,y,wy=1.0,par0=array([0.5,-1,.5,1,.5])):
    def minf(par): return( sum(wy*(dblgaussd(x,par)-y)**2)/sum(wy))
    minx  = minimize(minf,par0)
    return minx

nx = int(sqrt(len(roirs)))
npix = 15
nsize = nx*npix


movie = zeros((6000,nsize,nsize))

for i,roi in enumerate(roirs):
    roia = ftvaugmentr(roi,magn=2)
      
    roia = (roia-roia.min())/roia.max()
    ix = i%nx
    iy = i//nx

    yl1,yl2 = (iy*npix+2,(iy+1)*npix-2) 
    xl1,xl2 = (ix*npix+2,(ix+1)*npix-2)
    movie[:,yl1:yl2,xl1:xl2] = roia

    
    
from matplotlib import animation

fig,ax = subplots(nrows=2,ncols=1)
ax0 = ax[0]
ax0.set_xlim(0,1)
ax0.set_ylim(340,1400)

imt = movie[0,:,:]
im = ax[1].imshow(imt,cmap='gray')
im.set_clim(movie[1:,:,:].min(),movie[1:,:,:].max())

St = zeros((6000,2))
St[:,0] = arange(6000)*10e-3
for i in range(6000):
    imt = movie[i,:,:]
    im.set_array(imt)
    St[i,1] = (imt-imt.min()).sum()
ax0.set_xlim(0,1)

ax0.set_ylim(min(St[:,1]),max(St:,1]))

def init():
    im.set_data(movie[0,:,:])
    return im,
    
# animation function.  This is called sequentially
def animate(i):
    imt = movie[i,:,:]
    im.set_array(imt)
    ax0.set_xlim(max(0,i-80)*10e-3,max(i+20,100)*10e-3)
    St[i,1] = (imt-imt.min()).sum()
    #sc.set_offsets(St[max(0,i-80):(i+1),:])
    sc.set_array(St[max(0,i-80):(i+1),:])

    #title("t = %.1f" % tt)
    return im,


anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=6000, interval=10, blit=True,repeat_delay = 10000)

# animation function.  This is called sequentially
def animate(i):
    im.set_array(movie[i,:,:])
    return [im]

anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=6000, interval=10, blit=True,repeat_delay = 10000)
# interval in ms : each cycle is 1 seconds now!

movie = zeros((6000,nsize,nsize))
for i in range(1600):
    #roia = ftvaugmentr(roi,magn=2)
    
    ix = i%nx
    iy = i//nx

    yl1,yl2 = (iy*npix+2,(iy+1)*npix-2) 
    xl1,xl2 = (ix*npix+2,(ix+1)*npix-2)
    roi = movie[:,yl1:yl2,xl1:xl2] 
    roir = (roi-roi.min())/roi.mean()
    movie[:,yl1:yl2,xl1:xl2] = roir
    
    
for i,roi in enumerate(roirs):
    roia = ftvaugmentr(roi,magn=2)
    
    ix = i%nx
    iy = i//nx

    yl1,yl2 = (iy*npix+2,(iy+1)*npix-2) 
    xl1,xl2 = (ix*npix+2,(ix+1)*npix-2)
    movie[:,yl1:yl2,xl1:xl2] = roia



    roia = roia.reshape(6000,121)
    
    
    rrmin = roi.reshape(6000,25).min(axis=1)
    amp = (roi.reshape(6000,25).sum(axis=1)-rrmin*25)/rrmin/25.0

    sel = isfinite(amp)     
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
        
    seli = ~sel
    x0 = zeros(6000)
    x0[sel] = 1
    x0 = x0.reshape(1500,4)
    x0[x0.sum(axis=1)!=4,:] = 0

    
