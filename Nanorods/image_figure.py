from numpy import *
import matplotlib
from matplotlib.pylab import *
import pytiff
import os
from sub.subs import readBigtifFile,visualization



if len(sys.argv)>1:
    cfile = sys.argv[1] 
else:
    cfile = 0
    wdir = "./"
    print("No argument, we do it for file 0")


print("We will do analysis for file "+cfile)

    
wdir = ''
for fs in cfile.split("/")[:-1]:
    wdir = wdir+fs+'/'

resdir = cfile.split(".")[0]+'output/'
sptdir = cfile.split(".")[0]+'output/sptrack/'
if not os.path.isdir(sptdir):
    print("No folder with sptrack, you probably need to run features ROIS")
    raise error
    
cellname = cfile.split("/")[-1].split("_")[0]

# ----------------------------------------------------------------
# We read the file.     
    
# ~ f = ifile
movie = readBigtifFile(cfile)
mask = pytiff.Tiff(wdir+cellname+"_Mask_CMOS.tif")[:,:]

data = loadtxt(sptdir+"/ROIsfeatures.dat")

posA = loadtxt(resdir+"FposA.dat")

posB = loadtxt(resdir+"FposB.dat")

#pos = row_stack((posA,posB+[512,0]))
pos = posA

projimg = movie[1:,:,:].sum(axis=0)

imt = 1.0*projimg[:,:512];
visualization(imt,pos,widthr=5,heightr=5)
imshow(mask,alpha=0.3,cmap='Blues_r')


figure()
contrastd=1.0;contrastu=2.5
imt = 1.0*projimg[:,:512]; me = mean(imt.flatten()); sd = std(imt.flatten());
imt = imt/me
immin = max(imt.min(),1.0-contrastd*sd/me)
immax = min(imt.max(),1.0+contrastu*sd/me)

imt[imt<immin] = 0.0;
imt[(imt<immax)*(imt>immin)] = (imt[(imt<immax)*(imt>immin)]-immin)/(immax-immin)
imt[imt>immax] = 1.0;

imshow(imt,'gray')
cmapgnu = get_cmap("autumn")

sel = isfinite(data[:674,0])
scatter(pos[sel,0],pos[sel,1],c=6.0+log10(data[:674,0][sel]),cmap=cmapgnu,alpha=0.5)


def visualization(imA,pos=None,feature = None,widthr=5,heightr=5,contrastd=1.0,contrastu=2.5,figname="",color='red'):
    fig = figure(figname,figsize = (7,4.5))
    ax = fig.add_subplot(111)
    imt = imA*1.0; me = mean(imt.flatten()); sd = std(imt.flatten());
    imt[imt>me+contrastu*sd] = me+contrastu*sd; imt[imt<me-contrastd*sd] = me-contrastd*sd;
    ax.imshow(imt,cmap='gray')
    if pos is None:
        pass
    else:
        for i,a  in enumerate(pos):
            a_x, a_y = a
            fval = feature[i]
            ax.add_patch(Rectangle(xy=(a_x-widthr/2, a_y-heightr/2) ,width=widthr, height=heightr, linewidth=1.5, color=color,alpha=0.5, fill=))
    axis('off')
    subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0)
    margins(0,0)
    gca().xaxis.set_major_locator(NullLocator())
    gca().yaxis.set_major_locator(NullLocator())

    return((fig,ax))

figure()


