import numpy as np
#import pytiff
import PIL.Image as Image
import PIL.ExifTags as Exiftags
from codecs import open as copen

def readBigTifFile(fname):
    img = Image.open(fname)
    imf = img.tag_v2
    keys = []
    for key in imf.keys():
        keys.append(key)

    width, height = imf[keys[0]],imf[keys[1]]
    nframes = int(imf[270].split()[1][7:])
    size = nframes*width*height
    
    offset = imf[273][0]
    with open(fname,"rb") as file:
        file.seek(offset)
        temp = fromfile(file,dtype=">u2")

    temp = temp[:size].reshape((nframes,height,width))
    return(temp)
    
def readtifImage(fname):    
    return( np.array(Image.open(fname)))
    
def readtifInfo(fname,verbose=True,maxlength = 500):
    img = Image.open(fname)
    imf = img.tag_v2
    #keys = []
    info = {}
    for key in imf.keys():
        if len(str(imf[key]))<maxlength:
            #keys.append(key)
            info[key] = imf[key]
            print(key,':',imf[key])
        else:
            print(key,': _too long to print and to store_')
    return(info)
    
def gettimes(fname,nframes=None):
    if nframes is None:
        try:
            tags = readtifInfo(fname,verbose=False)
            nframes = int([d for d in info[270].split('\n') \
                         if d.find('frames')>=0][0].split('=')[-1])

            # ~ with pytiff.Tiff(fname) as handle:
                # ~ tags = handle.read_tags()
                # ~ nframes = int(tags['image_description'].split()[2][7:])
        except:
            print("I could not get the number of frames, please provide it")
            return

    with copen(fname,"r","windows-1252") as f:
        j = 0
        times = zeros((nframes))
        while True:
            try:
                # ~ print('linea ',j)
                line = f.readline()
                # ~ print(line)
                linesp = line.replace('\x00','').strip().split()
                if len(linesp)==5:
                    if linesp[2]=='Time_From_Last':
                        k,t = linesp[1],linesp[-1]
                        # print(int(k),float(t))
                        times[int(k)-1] = float(t)
                        j = j+1
            except Exception as e:
                print(e)
                break
            if j>=nframes:
                break
    return(times)
