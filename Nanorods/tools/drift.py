%pylab
a = loadtxt("driftsc.dat")

fig = figure()
ax = fig.add_subplot(111,projection="polar")

r = sqrt(a[:,0]**2+a[:,2]**2)
phi = arctan2(a[:,2],a[:,0])

ax.scatter(phi,r)
ax.set_ylabel("nm/s")
savefig("driftallfields.png")
