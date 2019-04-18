import os, sys
sys.path.append("/export/home1/users/bssn/serna/GitIBENS/Nanorods")




wdirt = "selpeaks/popt/"
files = os.listdir(wdirt)
dfilesr = []
for f in files:
    if f[:5] == "popts":
        dfilesr.append(wdirt+f)


dfilesr.sort()

from numpy import loadtxt,arange,load,save


drift = loadtxt("../drifts0.dat")

wdirt2 = 'selpeaks/poptC/'

dfilesr2 = []
for i,fname in enumerate(dfilesr):
    hdrift = drift[i,:]
    t = arange(6000)
    popt = load(fname)
    popt[:,:,5] = popt[:,:,5]-t*hdrift[0]
    popt[:,:,6] = popt[:,:,6]-t*hdrift[1]
    save(wdirt2+'popt'+str(i).zfill(2)+'.npy',popt)
    dfilesr2.append(wdirt2+'popt'+str(i).zfill(2)+'.npy')

outputdir = wdirt+'movies/'
for i,fname in enumerate(dfilesr):
    popt = load(fname)
    sh2 = popt.shape[0]//2
    for j in range(sh2):
        try:
            fig = plotNR(popt[j,:,:],ns = 25)
            fig.savefig(outputdir+"finalB"+str(i).zfill(2)+"_"+str(j).zfill(3)+".png",bbox_inches='tight',pad_inches = 0,facecolor=((0.4,0.4,0.4)))
            fig.clear()
        except:
            pass
    for j in range(sh2):
        try:
            fig = plotNR(popt[j+sh2,:,:],ns = 25)
            fig.savefig(outputdir+"finalR"+str(i).zfill(2)+"_"+str(j).zfill(3)+".png",bbox_inches='tight',pad_inches = 0,facecolor=((0.4,0.4,0.4)))
            fig.clear()
        except:
            pass

        try:
            disp = popt[j,:,5:7].mean(axis=0)-popt[j+sh2,:,5:7].mean(axis=0)
            popt[j+sh2,:,5] += disp[0]
            popt[j+sh2,:,6] += disp[1]
            fig = plotNR(popt[j,:,:] + popt[j+sh2,:,:],ns = 25)
            fig.savefig(outputdir+"finalT"+str(i).zfill(2)+"_"+str(j).zfill(3)+".png",bbox_inches='tight',pad_inches = 0,facecolor=((0.4,0.4,0.4)))
            fig.clear()
        except: 
            pass
        
        
outputdir = wdirt2+'movies/'
for i,fname in enumerate(dfilesr):
    popt = load(fname)
    sh2 = popt.shape[0]//2
    for j in range(sh2):
        try:
            fig = plotNR(popt[j,:,:],ns = 25)
            fig.savefig(outputdir+"finalB"+str(i).zfill(2)+"_"+str(j).zfill(3)+".png",bbox_inches='tight',pad_inches = 0,facecolor=((0.4,0.4,0.4)))
            fig.clear()
        except:
            pass
    for j in range(sh2):
        try:
            fig = plotNR(popt[j+sh2,:,:],ns = 25)
            fig.savefig(outputdir+"finalR"+str(i).zfill(2)+"_"+str(j).zfill(3)+".png",bbox_inches='tight',pad_inches = 0,facecolor=((0.4,0.4,0.4)))
            fig.clear()
        except:
            pass

        try:
            disp = popt[j,:,5:7].mean(axis=0)-popt[j+sh2,:,5:7].mean(axis=0)
            popt[j+sh2,:,5] += disp[0]
            popt[j+sh2,:,6] += disp[1]
            fig = plotNR(popt[j,:,:] + popt[j+sh2,:,:],ns = 25)
            fig.savefig(outputdir+"finalT"+str(i).zfill(2)+"_"+str(j).zfill(3)+".png",bbox_inches='tight',pad_inches = 0,facecolor=((0.4,0.4,0.4)))
            fig.clear()
        except: 
            pass
