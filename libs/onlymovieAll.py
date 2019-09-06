from numpy import *
import matplotlib
matplotlib.use('Agg')
from matplotlib.pylab import *
import pytiff
from scipy.optimize import minimize,least_squares
import os
from sub.subs import readBigtifFile,visualization

# ~ def f0(x,par): return(exp(par[0])+par[1]*x+exp(par[2]-x/par[3]))

# ~ def funcfit(fun,x,y,par0,ey=1.0):
    # ~ def minf(par): return( sum((fun(x,par)-y)**2/ey**2))
    # ~ minx  = minimize(minf,par0)
    # ~ return minx

wdir = "/export/home1/users/bssn/serna/anastasia/ROIS_raw/FullMovies/"


if len(sys.argv)>1:
    ifile = sys.argv[1] 
    # ~ print("Trying to print args:",sys.argv[1],sys.argv[2])
    print(ifile)

    try:
        ifile = int(ifile)

    except:
        print("Did you provide the name of the file? We are assuming this")
        f = ifile
        
        wdir = ''
        for fs in f.split("/")[:-1]:
            wdir = wdir+fs+'/'
        
        posdir = f.split(".")[0]+'output/'
        outputdir = f.split(".")[0]+'output/gif/'
        if not os.path.isdir(outputdir):
            try:
                os.system("mkdir "+outputdir)
            except ValueError:
                print("I cannot create the folder for the output!")
                raise SystemExit(0)
 
    # ~ if len(sys.argv)>2:
        # ~ wdir = sys.argv[2]
    # ~ else:
        # ~ wdir = "./"
else:
    ifile = 0
    wdir = "./"
    print("No argument, we do it for file number 0")


# ----------------------------------------------------------------
# We read the file.     
    
# ~ f = ifile
movie = readBigtifFile(f)


pos = loadtxt(posdir+"FposA.dat")

nframes,h,w2 = movie.shape
w = w2//2
for il, i in enumerate(range(1,nframes,100)):
    imt = 1.0*movie[i,:,:w];
    visualization(imt,pos,widthr=5,heightr=5)
    savefig(outputdir+"imA_"+str(il).zfill(2)+".png",bbox_inches='tight',pad_inches = 0)
    close()
    if il >120:
        break

pos = loadtxt(posdir+"FposB.dat")

for il, i in enumerate(range(1,nframes,100)):
    imt = 1.0*movie[i,:,w:];
    visualization(imt,pos,widthr=5,heightr=5)
    savefig(outputdir+"imB_"+str(il).zfill(2)+".png",bbox_inches='tight',pad_inches = 0)
    close()
    if il >120:
        break
#cat movie/im*.png | ffmpeg -framerate 24 -i - movie/vid"+str(i).zfill(3)+".mp4

os.system('convert '+outputdir+'imA_* '+outputdir+'imA.gif')
os.system('convert '+outputdir+'imB''_* '+outputdir+'imB.gif')


os.system('rm '+outputdir+'imA_*png')
os.system('rm '+outputdir+'imB_*png')
os.system('convert '+outputdir+'imA.gif -fuzz 10% -layers Optimize '+outputdir+'optimA.gif')
os.system('convert '+outputdir+'imB.gif -fuzz 10% -layers Optimize '+outputdir+'optimB.gif')
