import sys
import os
from numpy import *
from matplotlib.pylab import *
import sys
import os

# ~ wdirs = ["/mnt/data/Anastasia/18_11_29_pd23_11_div6_25Hzsqwave/",
# ~ "/mnt/data/Anastasia/QDs_prelabsem1811/","/mnt/data/Anastasia/Initial/"]

# ~ dfiles = []
# ~ for dirt in wdirs:
    # ~ basedir = dirt
    # ~ files = os.listdir(basedir)
    # ~ if dirt[-1] != '/':
        # ~ dirt = dirt+'/'
    # ~ for f in files:
        # ~ if f[-4:]=='.tif': 
            # ~ try:
                # ~ i = int(f[-5])
                # ~ dfiles.append(dirt+f)
            # ~ except:
                # ~ pass
    

if len(sys.argv)>1:
    wdir = sys.argv[1]
    if wdir[-1] != '/':
        wdir = wdir+'/'
else:
    wdir = "./"
    print("No argument, we do it for current folder")

outputdir = wdir+'ROISgifs/'
if not os.path.isdir(outputdir):
    try:
        os.system("mkdir "+outputdir)
    except ValueError:
        print("I cannot create the folder for the output GIFs!")
        raise SystemExit(0)


dirt = wdir

basedir = dirt
files = os.listdir(basedir)

dfiles = []
for f in files:
    if f[:5]=='roi_s' and f[-3:] == 'npy': dfiles.append(basedir+f)
 
dfiles.sort()

for filepath in dfiles:
    print("Starting with "+filepath)
    namesp = (filepath.split(sep=".")[0]).split(sep="/")[-1]
    rois = load(filepath)
    sh = shape(rois)
    if len(sh)>2 :
        print("This is a different one: "+filepath)
    #ion(); fig = figure(1)
    if sh[1] == 25:
        rois = rois.reshape(sh[0],5,5)
    elif sh[1] == 9:
        rois = rois.reshape(sh[0],3,3)
    else:
        ssh = sqrt(sh[1])
        if abs(ssh-floor(ssh))<1.0e-8:
            ssh = int(ssh)
            rois.reshape(sh[0],ssh,ssh)
        else:
            print("Second dimension is not perfect square? Is it a square ROI?")
            break
            
    me = mean(rois.flatten()); sd = std(rois.flatten());
    #Monitor specific!!!?!?!?!
    my_dpi = 96
    fig = figure(1,figsize=(100/my_dpi, 100/my_dpi), dpi=my_dpi)
    contrastu = 2.5
    contrastd = 1.0
    for i in range(1,sh[0],100):
        ax = fig.add_subplot(111)
        imt = (rois[i:(i+100),:,:]*1.0).min(axis=0)
        #imt = rois[i,:,:]*1.0
        imt[imt>me+contrastu*sd] = me+contrastu*sd; imt[imt<me-contrastd*sd] = me-contrastd*sd;
        ax.imshow(imt,cmap='gray')
        axis('off')
        subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
                hspace = 0, wspace = 0)
        margins(0,0)
        gca().xaxis.set_major_locator(NullLocator())
        gca().yaxis.set_major_locator(NullLocator())
        #pause(0.01)
        savefig(outputdir+"imt_"+str(i).zfill(5)+".png",bbox_inches='tight',pad_inches = 0,dpi = my_dpi)
        fig.clear()
    
    close(1)
    
    os.system('convert '+outputdir+'imt_* '+outputdir+'imt.gif')
    #os.system('convert '+outputdir+'imt.gif -fuzz 20% -layers Optimize '+outputdir+'im'+namesp+'.gif')
    os.system('cp '+outputdir+'imt.gif '+outputdir+'im'+namesp+'.gif')

    os.system('rm '+outputdir+'imt*')
    
for i in range(0,len(dfiles),9):
    names = []
    for j in range(i,i+9):
        filepath = dfiles[j]
        namesp = (filepath.split(sep=".")[0]).split(sep="/")[-1]
        namesp = outputdir+'im'+namesp+'.gif'
        names.append(namesp)
    
    out = []
    for j in range(0,9,3):    
        namer0 = names[j]
        namer1 = names[j+1]
        namer2 = names[j+2]
        outfile = outputdir+'temp.gif'
        outfile2 = outputdir+'temp'+str(j//3+1)+'.gif'
        os.system('convert '+namer0+' -repage 214x107 -coalesce null: \( '+namer1+' -coalesce \) -geometry +107+0 -layers Composite '+outfile)
        os.system('convert '+outfile+' -repage 321x107 -coalesce null: \( '+namer2+' -coalesce \) -geometry +214+0 -layers Composite '+outfile2)
        out.append(outfile2)
        
    outfile2 = outputdir+'collage_'+str(i).zfill(3)+'.gif'
    os.system('convert '+out[0]+' -repage 321x214 -coalesce null: \( '+out[1]+' -coalesce \) -geometry +0+107 -layers Composite '+outfile)
    os.system('convert '+outfile+' -repage 321x321 -coalesce null: \( '+out[2]+' -coalesce \) -geometry +0+214 -layers Composite '+outfile2)
