# -*- coding: utf-8 -*-
"""
Tools for handling the OSU tidal prediction software (OTPS) output data
(http://volkov.oce.orst.edu/tides/)

This software is based on the tide model driver (TMD) matlab code from here:
	http://polaris.esr.org/ptm_index.html

Matt Rayson
Stanford University
March 2013
"""

import os
import numpy as np
from scipy import interpolate

import collections
import datetime
import numbers

import pdb
import warnings

otis_constits = { 'M2':{'index':1,'omega':1.405189e-04,'v0u':1.731557546},\
     'S2':{'index':2,'omega':1.454441e-04,'v0u':0.000000000},\
     'N2':{'index':3,'omega':0.00013787970,'v0u':6.050721243},\
     'K2':{'index':4,'omega':0.0001458423,'v0u':3.487600001},\
     'K1':{'index':5,'omega':7.292117e-05,'v0u':0.173003674},\
     'O1':{'index':6,'omega':6.759774e-05,'v0u':1.558553872},\
     'P1':{'index':7,'omega':7.252295e-05,'v0u':6.110181633},\
     'Q1':{'index':8,'omega':6.495854e-05,'v0u':5.877717569}}


def _preprocess_time(time):
	"""
	See if we can massage time to be something useful
	"""
	print(time, type(time))
	if not isinstance(time, collections.Iterable):
		time = [time]
	if isinstance(time[0], np.datetime64):
		return time
	if isinstance(time[0], datetime.datetime):
		return np.array([np.datetime64(d) for d in time])
	if isinstance(time[0], numbers.Number):
		return np.array([np.datetime64(datetime.datetime.fromordinal(d))
							for d in time])
	raise TypeError('time must be np.datetime64, '
				    'list of datetime objects, or ordinal floats')

def tide_pred(modfile, lon, lat, time, conlist=None):
	"""
	Performs a tidal prediction at all points in [lon,lat] at times.

	Parameters
	----------
	modfile : string
		(Relative) path of the OSU model file on your file system
	lon, lat : array-like
		each is an n-length array of longitude and latitudes in units of
		degrees to perform predictions at. Usually a numpy array, though a
		list will work.  lat ranges from -90 to 90.  lon can range
		from -180 to 360.
  	time : array-like
		m-length array of times.  Acceptable formats are a list of `datetime`
		objects, a list or array of `numpy.datetime64` objects, or
		a list of floats, in which case they are assumed to be ordinal
		dates (See `.datetime.datetime.fromordinal`).
	conlist : list of strings (optional)
		If supplied, gives a list of tidal constituents to include in
		prediction. Available are 'M2', 'S2', 'N2', 'K2', 'K1', 'O1', 'P1',
		and 'Q1'.

	Returns
	-------
	h : m-by-n numpy array of tidal heights
		height is in meters, times are along the rows, and positions along
		the columns
	u : m-by-n numpy array of east-west tidal velocity [m/s]
	v : m-by-n numpy array of north tidal velocity [m/s]

	Examples
	--------

	dates = np.arange(np.datetime64('2001-04-03'),
	                  np.datetime64('2001-05-03'), dtype='datetime64[h]' )

	lon = np.array([198, 199])
	lat = np.array([21, 19])

	h, u, v = otp.tide_pred(modfile, lon, lat, dates,conlist=None)

	"""

	time = _preprocess_time(time)

	# Read and interpolate the constituents
	u_re, u_im, v_re, v_im, h_re, h_im, omega, conlist = extract_HC(modfile,
					lon, lat, conlist=conlist)

	# Initialise the output arrays
	sz = lon.shape
	nx = np.prod(sz)
	nt = time.shape[0]
	ncon = omega.shape[0]

	h_re = h_re.reshape((ncon,nx))
	h_im = h_im.reshape((ncon,nx))
	u_re = u_re.reshape((ncon,nx))
	u_im = u_im.reshape((ncon,nx))
	v_re = v_re.reshape((ncon,nx))
	v_im = v_im.reshape((ncon,nx))

	# Calculate nodal correction to amps and phases
	#t1992 = othertime.SecondsSince(time[0],basetime=datetime(1992,1,1))/86400.0

    #pu,pf,v0u = nodal(t1992+48622.0,conlist)

    # Calculate the time series
    #tsec = othertime.SecondsSince(time,basetime=datetime(1992,1,1))


	# nodal needs days since 1992:
	days = (time[0].astype('datetime64[D]') -
			np.datetime64('1992-01-01', 'D')).astype(float)

	pu,pf,v0u = nodal(days + 48622.0, conlist)

	# Calculate the time series
	tsec = (time.astype('datetime64[s]') -
			np.datetime64('1992-01-01', 's')).astype(float)

	h=np.zeros((nt,nx))
	u=np.zeros((nt,nx))
	v=np.zeros((nt,nx))
	for nn,om in enumerate(omega):
	    for ii in range(0,nx):
	        h[:,ii] += pf[nn]*h_re[nn,ii] * np.cos(om*tsec + v0u[nn] + pu[nn]) - \
	            pf[nn]*h_im[nn,ii] * np.sin(om*tsec + v0u[nn] + pu[nn])

	        u[:,ii] += pf[nn]*u_re[nn,ii] * np.cos(om*tsec + v0u[nn] + pu[nn]) - \
	            pf[nn]*u_im[nn,ii] * np.sin(om*tsec + v0u[nn] + pu[nn])

	        v[:,ii] += pf[nn]*v_re[nn,ii] * np.cos(om*tsec + v0u[nn] + pu[nn]) - \
	            pf[nn]*v_im[nn,ii] * np.sin(om*tsec + v0u[nn] + pu[nn])

	szo = (nt,)+sz
	return h.reshape(szo), u.reshape(szo), v.reshape(szo)


def extract_HC(modfile, lon, lat, z=None, conlist=None):
	"""
	Extract harmonic constituents from OTIS binary output and interpolate onto points in lon,lat

	set "z" to specifiy depth for transport to velocity conversion

	set "constituents" in conlist

	Returns:
	    u_re, u_im, v_re, v_im, h_re, h_im, omega, conlist

	"""

	###
	# Make sure the longitude is between 0 and 360
	lon = np.mod(lon,360.0)

	###
	# Read the filenames from the model file
	pathfile = os.path.split(modfile)
	path = pathfile[0]

	f = open(modfile,'r')
	hfile = path+'/' + f.readline().strip()
	uvfile = path+'/' + f.readline().strip()
	grdfile = path+'/' + f.readline().strip()
	f.close()

	###
	# Read the grid file
	X, Y, depth, mask = read_OTPS_grd(grdfile)
	#X[X>180.0] = 180.0 - X[X>180.0]
	mask = mask == 1
	def interpit(Z, lon, lat):
		spl = interpolate.RectBivariateSpline(X[0, :], Y[:, 0], Z.T)
		return spl(lon, lat, grid=False)

	# Create an interpolation object
	sz = lon.shape
	lon = lon.ravel()
	lat = lat.ravel()
	nx = lon.size
	z = interpit(depth, lon, lat)


	###
	# Check that the constituents are in the file
	conOTIS = get_OTPS_constits(hfile)

	if conlist == None:
	    conlist = conOTIS

	for vv in conlist:
	    if not vv in conOTIS:
	        warnings.warn('Constituent name: %s not present in OTIS file.'%vv)
	        conlist.remove(vv)

	###
	# Now go through and read the data for each

	# Initialse the arrays
	ncon = len(conlist)
	u_re = np.zeros((ncon,nx))
	u_im = np.zeros((ncon,nx))
	v_re = np.zeros((ncon,nx))
	v_im = np.zeros((ncon,nx))
	h_re = np.zeros((ncon,nx))
	h_im = np.zeros((ncon,nx))
	omega = np.zeros((ncon,))

	for ii, vv in enumerate(conlist):
		idx = otis_constits[vv]['index']
		omega[ii] = otis_constits[vv]['omega']
		print('Interpolating consituent: %s...'%vv)

		# Read and interpolate h
		X ,Y, tmp_h_re, tmp_h_im = read_OTPS_h(hfile,idx)

		h_re[ii, :] = interpit(tmp_h_re, lon, lat)
		h_im[ii, :] = interpit(tmp_h_im, lon, lat)

		# Read and interpolate u and v - Note the conversion from transport to velocity
		X ,Y, tmp_u_re, tmp_u_im, tmp_v_re, tmp_v_im = read_OTPS_UV(uvfile,idx)
		u_re[ii, :] = interpit(tmp_u_re, lon, lat) / z
		u_im[ii, :] = interpit(tmp_u_im, lon, lat) / z
		v_re[ii, :] = interpit(tmp_v_re, lon, lat) / z
		v_im[ii, :] = interpit(tmp_v_im, lon, lat) / z

	# Return the arrays in their original shape
	szout = (ncon,) + sz
	return u_re.reshape(szout), u_im.reshape(szout), v_re.reshape(szout), \
	    v_im.reshape(szout), h_re.reshape(szout), h_im.reshape(szout), omega, conlist


def nodal_correction(year,conlist,amp, phase):
        """
		### UNUSED ###
        Applies a lunar nodal correction to the amplitude and phase

        Code modified from Rusty Holleman's GET_COMPONENTS code below...
        #
        # GET_COMPONENTS
        #   [UAMP,UPHASE,VAMP,VPHASE,HAMP,HPHASE]=GET_COMPONENTS(YEAR,OMEGAT,LUN_NODE,V0U,AG)
        #   calculates the tidal amplitudes and phases from the interpolated OTIS
        #   data in the AG matrix.
        #
        #   This code has been adapted from Brian Dushaw's matlab scripts
        #   obtained from http://909ers.apl.washington.edu/~dushaw/tidegui/tidegui.html
        #
        #function [uamp,uphase,vamp,vphase,hamp,hphase]=get_components(YEAR,omegat,lun_node,v0u,AG)
        """

        import tide_consts as tc

        #oneday=np.array( [335.62,  0,  322.55,  1.97,   334.63,  0.99,  -0.99,  321.57])
        oneday = {'M2':335.62, 'S2':0, 'N2':322.55, 'K2':1.97, 'O1':334.63, 'K1':0.99, 'P1':-0.99, 'Q1':321.57}
        #if year < 1970 or year > 2037:
        #    print 'Constants for prediction year are not available'
            #return None

        # Find the index
        JJ=[]
        od=np.zeros((len(conlist),1))
        for ii,vv in enumerate(conlist):

            jj=[item for item in range(len(tc.const_names)) if tc.const_names[item] == vv]
            if len(jj) > 0:
                JJ.append(jj)

            if oneday.has_key(vv):
                od[ii]=(np.pi/180)*oneday[vv]

        I = int( np.where(year==tc.years)[0] )
        vou=tc.v0u[JJ,I]
        lunnod=tc.lun_nodes[JJ,I]

        vou=(np.pi/180)*vou


        #oneday=(np.pi/180)*oneday

        hamp = amp*lunnod
        #hphase = - oneday[JJ] + vou[JJ] - G
        hphase =  -od + vou - phase

        return hamp, hphase


def read_OTPS_UV(uvfile,ic):
    """
    Reads the tidal transport constituent data from an otis binary file

    ic = constituent number

    Returns: X, Y, h_re and h_im (Real and imaginary components)

    See this post on byte ordering
         http://stackoverflow.com/questions/1632673/python-file-slurp-w-endian-conversion
    """
    f = open(uvfile,'rb')
    #f = hfile
    # Try numpy
    ll = np.fromfile(f,dtype=np.int32,count=1)
    nm = np.fromfile(f,dtype=np.int32,count=3)
    th_lim = np.fromfile(f,dtype=np.float32,count=2)
    ph_lim = np.fromfile(f,dtype=np.float32,count=2)

    # Need to go from little endian to big endian
    ll.byteswap(True)
    nm.byteswap(True)
    th_lim.byteswap(True)
    ph_lim.byteswap(True)

    n = nm[0]
    m = nm[1]
    nc = nm[2]

    if ic < 1 or ic > nc:
        raise Exception('ic must be > 1 and < %d'%ic)

    # Read the actual data
    nskip = int((ic-1)*(nm[0]*nm[1]*16+8) + 8 + ll - 28)
    f.seek(nskip,1)

    htemp = np.fromfile(f,dtype=np.float32,count=4*n*m)
    htemp.byteswap(True)

    f.close()

    htemp = np.reshape(htemp,(m,4*n))
    U_re = htemp[:,0:4*n-3:4]
    U_im = htemp[:,1:4*n-2:4]
    V_re = htemp[:,2:4*n-1:4]
    V_im = htemp[:,3:4*n:4]

    X,Y = np.meshgrid(np.linspace(th_lim[0],th_lim[1],n),np.linspace(ph_lim[0],ph_lim[1],m))

    return X, Y, U_re, U_im, V_re, V_im


def read_OTPS_grd(grdfile):
	"""
	Reads the grid data from an otis binary file

	Returns: X, Y, hz, mask

	See this post on byte ordering
	     http://stackoverflow.com/questions/1632673/python-file-slurp-w-endian-conversion
	"""
	f = open(grdfile,'rb')
	#
	## Try numpy
	f.seek(4,0)
	n = np.fromfile(f,dtype=np.int32,count=1)
	m = np.fromfile(f,dtype=np.int32,count=1)
	lats = np.fromfile(f,dtype=np.float32,count=2)
	lons = np.fromfile(f,dtype=np.float32,count=2)
	dt = np.fromfile(f,dtype=np.float32,count=1)

	n.byteswap(True)
	m.byteswap(True)
	n = int(n)
	m = int(m)
	lats.byteswap(True)
	lons.byteswap(True)
	dt.byteswap(True)

	nob = np.fromfile(f,dtype=np.int32,count=1)
	nob.byteswap(True)
	if nob == 0:
		f.seek(20,1)
		iob = []
	else:
		f.seek(8,1)
		iob = np.fromfile(f, dtype=np.int32, count=int(2 * nob))
		iob.byteswap(True)
		iob = np.reshape(iob, (2, int(nob)))
		f.seek(8,1)

	hz = np.fromfile(f,dtype=np.float32,count=int(n * m))
	f.seek(8,1)
	mask = np.fromfile(f,dtype=np.int32,count=int(n * m))

	hz.byteswap(True)
	mask.byteswap(True)

	hz = np.reshape(hz,(m,n))
	mask = np.reshape(mask,(m,n))

	f.close()

	X,Y = np.meshgrid(np.linspace(lons[0],lons[1],n),np.linspace(lats[0],lats[1],m))

	return X, Y ,hz, mask


def read_OTPS_h(hfile,ic):
    """
    Reads the elevation constituent data from an otis binary file

    ic = constituent number

    Returns: X, Y, h_re and h_im (Real and imaginary components)

    See this post on byte ordering
         http://stackoverflow.com/questions/1632673/python-file-slurp-w-endian-conversion
    """
    f = open(hfile,'rb')
    #f = hfile
    # Try numpy
    ll = np.fromfile(f,dtype=np.int32,count=1)
    nm = np.fromfile(f,dtype=np.int32,count=3)
    th_lim = np.fromfile(f,dtype=np.float32,count=2)
    ph_lim = np.fromfile(f,dtype=np.float32,count=2)

    # Need to go from little endian to big endian
    ll.byteswap(True)
    nm.byteswap(True)
    th_lim.byteswap(True)
    ph_lim.byteswap(True)

    n = nm[0]
    m = nm[1]
    nc = nm[2]

    if ic < 1 or ic > nc:
        raise Exception('ic must be > 1 and < %d'%ic)
        #return -1

    # Read the actual data
    nskip = int((ic-1)*(nm[0]*nm[1]*8+8) + 8 + ll - 28)
    f.seek(nskip,1)

    htemp = np.fromfile(f,dtype=np.float32,count=2*n*m)
    htemp.byteswap(True)

    #
    f.close()

    htemp = np.reshape(htemp,(m,2*n))
    h_re = htemp[:,0:2*n-1:2]
    h_im = htemp[:,1:2*n:2]

    X,Y = np.meshgrid(np.linspace(th_lim[0],th_lim[1],n),np.linspace(ph_lim[0],ph_lim[1],m))

    return X ,Y, h_re, h_im


def get_OTPS_constits(hfile):
    """
    Returns the list of constituents in the file
    """
    f = open(hfile,'rb')
    ll = np.fromfile(f,dtype=np.int32,count=1)
    nm = np.fromfile(f,dtype=np.int32,count=3)

    ll.byteswap(True)
    nm.byteswap(True)

    f.close()

    ncon = nm[2]
    conList = []
    for ii in range(1,ncon+1):
        for vv in otis_constits:
            if otis_constits[vv]['index']==ii:
                conList.append(vv)

    return conList

def cart2pol(re,im):

    amp = np.abs(re + 1j*im)
    phs = np.angle(re + 1j*im)

    return amp, phs


def pol2cart(amp,phs):

    re = amp * np.cos(phs)
    im = amp * np.sin(phs)

    return re, im


def astrol(time):
    """
    %function  [s,h,p,N]=astrol(time);
    %  Computes the basic astronomical mean longitudes  s, h, p, N.
    %  Note N is not N', i.e. N is decreasing with time.
    %  These formulae are for the period 1990 - 2010, and were derived
    %  by David Cartwright (personal comm., Nov. 1990).
    %  time is UTC in decimal MJD.
    %  All longitudes returned in degrees.
    %  R. D. Ray    Dec. 1990
    %  Non-vectorized version. Re-make for matlab by Lana Erofeeva, 2003
    % usage: [s,h,p,N]=astrol(time)
    %        time, MJD
    circle=360;
    T = time - 51544.4993;
    % mean longitude of moon
    % ----------------------
    s = 218.3164 + 13.17639648 * T;
    % mean longitude of sun
    % ---------------------
    h = 280.4661 +  0.98564736 * T;
    % mean longitude of lunar perigee
    % -------------------------------
    p =  83.3535 +  0.11140353 * T;
    % mean longitude of ascending lunar node
    % --------------------------------------
    N = 125.0445D0 -  0.05295377D0 * T;
    %
    s = mod(s,circle);
    h = mod(h,circle);
    p = mod(p,circle);
    N = mod(N,circle);

    """

    circle=360;
    T = time - 51544.4993;
    # mean longitude of moon
    # ----------------------
    s = 218.3164 + 13.17639648 * T;
    # mean longitude of sun
    # ---------------------
    h = 280.4661 +  0.98564736 * T;
    # mean longitude of lunar perigee
    # -------------------------------
    p =  83.3535 +  0.11140353 * T;
    # mean longitude of ascending lunar node
    # --------------------------------------
    N = 125.0445 -  0.05295377 * T;
    #
    s = np.mod(s,circle);
    h = np.mod(h,circle);
    p = np.mod(p,circle);
    N = np.mod(N,circle);

    return s,h,p,N


def nodal(time,con):
	"""
	Nodal correction

	Derived from the tide model driver matlab scipt: nodal.m
	"""

	rad = np.pi/180.0

	s,h,p,omega=astrol(time)
	#
	#    omega =
	#
	#     determine nodal corrections f and u
	#     -----------------------------------
	sinn = np.sin(omega*rad);
	cosn = np.cos(omega*rad);
	sin2n = np.sin(2*omega*rad);
	cos2n = np.cos(2*omega*rad);
	sin3n = np.sin(3*omega*rad);

	ndict={'M2':{'f':np.sqrt((1.-.03731*cosn+.00052*cos2n)**2 + (.03731*sinn-.00052*sin2n)**2),\
	    'u':np.arctan((-.03731*sinn+.00052*sin2n)/(1.-.03731*cosn+.00052*cos2n))/rad},\
	    'S2':{'f':1.0, 'u':0.0},\
	    'K1':{'f':np.sqrt((1.+.1158*cosn-.0029*cos2n)**2 + (.1554*sinn-.0029*sin2n)**2),\
	        'u':np.arctan((-.1554*sinn+.0029*sin2n)/(1.+.1158*cosn-.0029*cos2n))/rad},\
	    'O1':{'f':np.sqrt((1.0+0.189*cosn-0.0058*cos2n)**2 + (0.189*sinn-0.0058*sin2n)**2),\
	        'u':10.8*sinn - 1.3*sin2n + 0.2*sin3n},\
	    'N2':{'f':np.sqrt((1.-.03731*cosn+.00052*cos2n)**2 + (.03731*sinn-.00052*sin2n)**2),\
	        'u':np.arctan((-.03731*sinn+.00052*sin2n)/(1.-.03731*cosn+.00052*cos2n))/rad},\
	    'P1':{'f':1.0, 'u':0.0},\
	    'K2':{'f':np.sqrt((1.+.2852*cosn+.0324*cos2n)**2 + (.3108*sinn+.0324*sin2n)**2),\
	        'u':np.arctan(-(.3108*sinn+.0324*sin2n) /(1.+.2852*cosn+.0324*cos2n))/rad},\
	    'Q1':{'f':np.sqrt((1.+.188*cosn)**2+(.188*sinn)**2),\
	        'u':np.arctan(.189*sinn / (1.+.189*cosn))/rad} }

	# Prepare the output data
	ncon = len(con)
	pu = np.zeros((ncon,1))
	pf = np.ones((ncon,1))
	v0u = np.zeros((ncon,1))
	for ii,vv in enumerate(con):
		if vv in ndict:
			pu[ii,:] = ndict[vv]['u']*rad
			pf[ii,:] = ndict[vv]['f']
		if vv in otis_constits:
			v0u[ii,:] = otis_constits[vv]['v0u']

	return pu, pf, v0u
