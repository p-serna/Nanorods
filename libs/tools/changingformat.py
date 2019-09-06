import h5py,os
from numpy import save,savetxt,array

dirs = ["./"]

for dirt in dirs:
    basedir = dirt
    files = os.listdir(basedir)

    dfiles = []
    for f in files:
        if f[-4:]=='.mat': dfiles.append(basedir+f)

arrays = {}
names = []
print("Data files:")
for filepath in dfiles:
    f = h5py.File(filepath)
    #filepath="stimulation.mat"

    name = filepath.split(sep=".")[1]
    for k,v in f.items():
        print(k,type(v))
        arrays[name+k] = array(v)
        names.append(name+k)
        print(name,k,v.shape)


for k in arrays.keys():
    #savetxt("."+k+".dat",arrays[k])
    save("."+k+".npy",arrays[k])


# ~ i = 0

# ~ for k in names:
    # ~ for j in range(arrays[k].shape[0]):
        # ~ savetxt("./temp/roi_sA"+str(i).zfill(4)+".dat",arrays[k][j,:,:])
        # ~ i += 1
