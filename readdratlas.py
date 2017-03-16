# Functions to read in data from the North American Drought Atlas,
# Monsoon Asia Drought Atlas, and Old World Drought Atlas over a specified
# region, time period; filename can also be specified arbitrarily
# Note: all times are in units of years, longitudes are 0-360E
# Sam Stevenson
# March 2017

def readnada(minlat=-90,maxlat=90,minlon=0,maxlon=360,tstart=1,tstop=2006,nadafile="/glade/u/home/samantha/PDSI/NADA/NADAv2-2008.nc"):
	print nadafile
	import numpy as np
	import netCDF4 as nc
	
	# Read in coordinate data
	dat = nc.Dataset(nadafile)
	nadatime=dat.variables['time'][:]
	nadalat=dat.variables['lat'][:]
	nadalon=np.add(dat.variables['lon'][:],360)
	
	# Find needed index locations
	latind=np.where((nadalat >= minlat) & (nadalat <= maxlat))
	lonind=np.where((nadalon >= minlon) & (nadalon <= maxlon))
	tind=np.where((nadatime >= tstart) & (nadatime <= tstop))
	
	# Read in PDSI, filter time
	nadatime=nadatime[tind]
	nadapdsi=dat.variables['PDSI'][list(tind[0]),list(latind[0]),list(lonind[0])]

	# NADA: default time axis is in descending order. Flip to ascending for intuitive use.
        ind=np.argsort(nadatime,axis=0)
	return nadatime[ind],nadalat[latind],nadalon[lonind],nadapdsi[ind,:,:];


def readmada(minlat=-90,maxlat=90,minlon=0,maxlon=360,tstart=1300,tstop=2005,madafile="/glade/u/home/samantha/PDSI/MADA/MADApdsi.nc"):
        print madafile
        import numpy as np 
        import netCDF4 as nc

	# Read in coordinate data
        dat = nc.Dataset(madafile)
        madatime=dat.variables['time'][:]
        madalat=dat.variables['lat'][:]
        madalon=dat.variables['lon'][:]
        
	# Find needed index locations
	latind=np.where((madalat >= minlat) & (madalat <= maxlat))
        lonind=np.where((madalon >= minlon) & (madalon <= maxlon))
	tind=np.where((madatime >= tstart) & (madatime <= tstop))

	# Read in PDSI, filter time
	madatime=madatime[tind]
        madapdsi=dat.variables['PDSI'][list(tind[0]),list(latind[0]),list(lonind[0])]

        # MADA: default time axis is in descending order. Flip to ascending for intuitive use.
        ind=np.argsort(madatime,axis=0)
        return madatime[ind],madalat[latind],madalon[lonind],madapdsi[ind,:,:];


def readowda(minlat=-90,maxlat=90,minlon=0,maxlon=360,tstart=1,tstop=2012,owdafile="/glade/u/home/samantha/PDSI/OWDA/data.nc"):
        print owdafile
        import numpy as np
        import netCDF4 as nc

	# Read in coordinate data
        dat = nc.Dataset(owdafile)
        owdatime=dat.variables['T'][:]    # months since 1960-01-01
        owdalat=dat.variables['Y'][:]
        owdalon=dat.variables['X'][:]

	# Find needed index locations
        latind=np.where((owdalat >= minlat) & (owdalat <= maxlat))
        lonind=np.where((owdalon >= minlon) & (owdalon <= maxlon))
	tind=np.where((owdatime >= tstart) & (owdatime <= tstop))

	# Read in PDSI, filter time
	owdatime=owdatime[tind]
        owdapdsi=dat.variables['pdsi'][list(tind[0]),list(latind[0]),list(lonind[0])]

        # OWDA: default time axis is in descending order. Flip to ascending for intuitive use.
        ind=np.argsort(owdatime,axis=0)
        return owdatime[ind],owdalat[latind],owdalon[lonind],owdapdsi[ind,:,:];


def readanzda(minlat=-90,maxlat=90,minlon=0,maxlon=360,tstart=1500,tstop=2012,anzdadir="/glade/u/home/samantha/PDSI/ANZDA/"):
	print anzdadir
	import numpy as np
	
	# Coordinates of ANZDA time series
	crds = np.loadtxt((anzdadir+"/anzda-pdsi-xy.txt"),usecols=(0,1))  # col 0 = lon; col 1 = lat

	# Read full PDSI dataset
	anzdat = np.loadtxt((anzdadir+"/anzda-recon.txt"))
	anzdatime=anzdat[:,0]

	# Filter lat/lon, time
        lind=np.where((crds[:,1] >= minlat) & (crds[:,1] <= maxlat) & (crds[:,0] >= minlon) & (crds[:,0] <= maxlon))
	tind=np.where((anzdatime >= tstart) & (anzdatime <= tstop))

	# Filter PDSI data
	anzdapdsi=anzdat[:,1:]
	anzdapdsi=anzdapdsi[tind[0],:]
	anzdapdsi=anzdapdsi[:,lind[0]]

	return anzdatime[tind],crds[lind,1],crds[lind,0],anzdapdsi
