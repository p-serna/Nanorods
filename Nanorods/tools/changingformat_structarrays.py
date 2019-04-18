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
        if k=="Exp_Data":
            for k2,v2 in v.items():
                print("----",k2,type(v2))
                if type(v2) == h5py._hl.dataset.Dataset:
                    #print(2*"----",v2.shape)
                    for i,v2e in enumerate(v2):
                        if v2e.shape[0] == 1:
                            try:
                                arname = k+k2+str(i).zfill(2)
                                #arname = arname[1:]
                                d = array(f[v2e[0]])
                                names.append(arname)
                                print(2*"----",arname)
                                arrays[arname] = d
                            except Exception as e:
                                print(e)
                                break
                                

for k in arrays.keys():
    #savetxt("."+k+".dat",arrays[k])
    save("./data/"+k+".npy",arrays[k])


# ~ i = 0

# ~ for k in names:
    # ~ for j in range(arrays[k].shape[0]):
        # ~ savetxt("./temp/roi_sA"+str(i).zfill(4)+".dat",arrays[k][j,:,:])
        # ~ i += 1
