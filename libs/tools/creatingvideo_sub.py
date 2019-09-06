import sys
sys.path.append("/users/bssn/serna/GitIBENS/Nanorods/")
from ROIS_analysis_all_speeding import fitROI
from numpy import *
from matplotlib.pylab import *
from matplotlib import animation

def ROIanimation(movie,dt = 10e-3,nwindow = 1000):
    nframes = movie.shape[0]
    St = zeros((nframes,2))
    St[:,0] = arange(nframes)*dt

    for i in range(nframes):
        imt = movie[i,:,:]
        St[i,1] = (imt-imt.min()).sum()
    St2 = column_stack(((St[:,0].reshape(nframes//4,4))[:,:2].flatten(),(St[:,1].reshape(nframes//4,4))[:,:2].flatten()))

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
    im = ax2.imshow(imt,cmap='hot')
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

    nwindow2 = nwindow*8//10
    def animate(i):
        thisx = St[max(0,i-nwindow2):(i+1),0]
        thisy = St[max(0,i-nwindow2):(i+1),1]

        imt = movie[i,:,:]
        im.set_array(imt)

        line.set_data(thisx, thisy)
        line2.set_data(St2[max(0,(i-nwindow2)//2):(i//2+1),0], St2[max(0,(i-nwindow2)//2):(i//2+1),1])

        time_text.set_text(time_template%(i*dt))
        ax.set_xlim(max(0,i-nwindow2)*dt,max(i+nwindow*2//10,nwindow)*dt)
        return line, time_text, im,ax,line2


    ani = animation.FuncAnimation(fig, animate, np.arange(nframes),
        interval=10, blit=False, init_func=init)
    return(St,ani,fig)
