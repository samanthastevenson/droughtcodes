# Functions to read in data from the North American Drought Atlas,
# Monsoon Asia Drought Atlas, and Old World Drought Atlas over a specified
# region, time period; filename can also be specified arbitrarily
# Note: all times are in units of years, longitudes are 0-360E
# Sam Stevenson
# March 2017

def readnada(minlat=-90,maxlat=90,minlon=0,maxlon=360,tstart=1,tstop=2006,file="/glade/u/home/samantha/PDSI/NADA/NADAv2-2008.nc"):
	print file
	import numpy as np
	import netCDF4 as nc
	
	dat = nc.Dataset(file)
	nadatime=dat.variables['time'][:]
	nadalat=dat.variables['lat'][:]
	nadalon=np.add(dat.variables['lon'][:],360)
	latind=np.where((nadalat >= minlat) & (nadalat <= maxlat))
	lonind=np.where((nadalon >= minlon) & (nadalon <= maxlon))
	
	# NADA: default time axis is in descending order. Flip to ascending for intuitive use.
	ind=np.argsort(nadatime,axis=0)
	nadapdsi=dat.variables['PDSI'][:,list(latind[0]),list(lonind[0])]
	
	return nadatime[ind],nadalat[latind],nadalon[lonind],nadapdsi[ind,:,:];


def readmada(minlat=-90,maxlat=90,minlon=0,maxlon=360,tstart=1300,tstop=2005,file="/glade/u/home/samantha/PDSI/MADA/MADApdsi.nc"):
        print file
        import numpy as np 
        import netCDF4 as nc

        dat = nc.Dataset(file)
        madatime=dat.variables['time'][:]
        madalat=dat.variables['lat'][:]
        madalon=dat.variables['lon'][:]
        latind=np.where((madalat >= minlat) & (madalat <= maxlat))
        lonind=np.where((madalon >= minlon) & (madalon <= maxlon))

        # MADA: default time axis is in descending order. Flip to ascending for intuitive use.
        ind=np.argsort(madatime,axis=0)
        madapdsi=dat.variables['PDSI'][:,list(latind[0]),list(lonind[0])]

        return madatime[ind],madalat[latind],madalon[lonind],madapdsi[ind,:,:];


def readowda(minlat=-90,maxlat=90,minlon=0,maxlon=360,tstart=1,tstop=2012,file="/glade/u/home/samantha/PDSI/OWDA/data.nc"):
        print file
        import numpy as np
        import netCDF4 as nc

        dat = nc.Dataset(file)
        owdatime=dat.variables['T'][:]    # months since 1960-01-01
        owdalat=dat.variables['Y'][:]
        owdalon=dat.variables['X'][:]
        latind=np.where((owdalat >= minlat) & (owdalat <= maxlat))
        lonind=np.where((owdalon >= minlon) & (owdalon <= maxlon))

        # OWDA: default time axis is in descending order. Flip to ascending for intuitive use.
        ind=np.argsort(owdatime,axis=0)
        owdapdsi=dat.variables['pdsi'][:,list(latind[0]),list(lonind[0])]

        return owdatime[ind],owdalat[latind],owdalon[lonind],owdapdsi[ind,:,:];


