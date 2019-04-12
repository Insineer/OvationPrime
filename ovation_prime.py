import os
from scipy import interpolate

from ovation_utilities import *

#Determine where this module's source file is located
#to determine where to look for the tables
src_file_dir = os.path.dirname(os.path.realpath(__file__))
ovation_datadir = os.path.join(src_file_dir,'data')

class FluxEstimator(object):
	"""
	A class which estimates auroral flux
	based on the Ovation Prime regressions,
	at arbitrary locations and times.
	"""
	def __init__(self, fluxType):
		self.fluxType = fluxType

		seasons = ['spring','summer','fall','winter']
		self.seasonal_flux_estimators = {season : SeasonalFluxEstimator(season, fluxType) for season in seasons}

	def get_season_fluxes(self,dF):
		"""
		Extract the flux for each season and hemisphere and
		store them in a dictionary
		Return positive latitudes, since northern and southern
		latitude/localtime grids are the same
		"""
		seasonfluxesN, seasonfluxesS = {}, {}
		gridmlats, gridmlts = None, None
		for season, estimator in self.seasonal_flux_estimators.items():
			flux_outs = estimator.get_gridded_flux(dF)
			gridmlatsN, gridmltsN, gridfluxN = flux_outs[:3]
			gridmlatsS, gridmltsS, gridfluxS = flux_outs[3:]
			seasonfluxesN[season]=gridfluxN
			seasonfluxesS[season]=gridfluxS
			gridmlats = gridmlatsN
			gridmlts = gridmltsN
		return gridmlats, gridmlts, seasonfluxesN, seasonfluxesS

	def get_flux_for_time(self, dt, hemisphere='N', combine_hemispheres=True):
		yday = dt.timetuple().tm_yday

		if hemisphere == 'N':
			weights = self.season_weights(yday)
		elif hemisphere == 'S':
			weights = self.season_weights(365 - yday)
		else:
			raise ValueError('Invalid hemisphere %s (use N or S)' % hemisphere)

		solarwind = CalculateSolarWindCouplingFunc(dt)

		season_fluxes_outs = self.get_season_fluxes(solarwind)
		grid_mlats, grid_mlts, seasonfluxesN, seasonfluxesS = season_fluxes_outs

		gridflux = np.zeros_like(seasonfluxesN['summer'])
		for season,W in weights.items():
			gridfluxN, gridfluxS = seasonfluxesN[season], seasonfluxesS[season]
			if combine_hemispheres:
				gridflux += W*(gridfluxN+gridfluxS)/2
			elif hemisphere == 'N':
				gridflux += W*gridfluxN
			elif hemisphere == 'S':
				gridflux += W*gridfluxS

		if hemisphere == 'S':
			grid_mlats = -1. * grid_mlats # by default returns positive latitudes

		return grid_mlats,grid_mlts,gridflux


	def season_weights(self,doy):
		"""
		Determines the relative weighting of the
		model coeffecients for the various seasons for a particular
		day of year (doy). Nominally, weights the seasons
		based on the difference between the doy and the peak
		of the season (solstice/equinox)

		Returns:
			a dictionary with a key for each season.
			Each value in the dicionary is a float between 0 and 1
		"""
		weight = {'winter':0.,'spring':0.,'summer':0.,'fall':0.}
		winter_w,spring_w,summer_w,fall_w = 0.,0.,0.,0.

		if doy >= 79. and doy < 171:
			weight['summer'] = 1. - (171.-doy)/92.
			weight['spring'] = 1. - weight['summer']

		elif doy >= 171. and doy < 263.:
			weight['fall'] = 1. - (263.-doy)/92.
			weight['summer'] = 1. - weight['fall']

		elif doy >= 263. and doy < 354.:
			weight['winter'] = 1. - (354.-doy)/91.
			weight['fall'] = 1. - weight['winter']

		elif doy >= 354 or doy < 79:
			#For days of year > 354, subtract 365 to get negative
			#day of year values for computation
			doy0 = doy- 365. if doy >= 354 else doy
			weight['spring'] = 1. - (79.-doy0)/90.
			weight['winter'] = 1. - weight['spring']

		return weight

class SeasonalFluxEstimator(object):
	"""
	A class to hold and calculate predictions from the regression coeffecients
	which are tabulated in the data/premodel/{season}_{atype}_*.txt
	files.

	Given a particular season, type of aurora ( one of ['diff','mono','wave','ions'])
	and type of flux, returns
	"""
	def __init__(self, season, fluxType):
		"""

		season - str,['winter','spring','summer','fall']
			season for which to load regression coeffients

		fluxType - str, ['diff','mono','wave','ions']
			type of aurora for which to load regression coeffients

		"""
		nmlt = 96                           # number of mag local times in arrays (resolution of 15 minutes)
		nmlat = 160                         # number of mag latitudes in arrays (resolution of 1/4 of a degree (.25))
		ndF = 12							# number of coupling strength bins
		self.fluxType = fluxType

		self.n_mlt_bins, self.n_mlat_bins, self.n_dF_bins = nmlt, nmlat, ndF

		#The mlat bins are orgainized like -50:-dlat:-90,50:dlat:90
		self.mlats = np.concatenate([np.linspace(-90, -50 , self.n_mlat_bins // 2)[::-1],
									 np.linspace(50, 90 , self.n_mlat_bins // 2)])

		self.mlts = np.linspace(0.,24.,self.n_mlt_bins)

		self.afile = os.path.join(ovation_datadir,'premodel/%s_%s.txt' % (season, fluxType))
		self.pfile = os.path.join(ovation_datadir,'premodel/%s_prob_b_%s.txt' % (season, fluxType))

		#self.valid_types = ['diff','mono','wave','ions']

		with open(self.afile,'r') as f:
			f.readline() # pass header
			adata = np.genfromtxt(f,max_rows=nmlat*nmlt)

		# fill model parameters
		self.b1a, self.b2a = np.zeros((nmlt,nmlat)), np.zeros((nmlt,nmlat))
		self.b1a.fill(np.nan)
		self.b2a.fill(np.nan)
		mlt_bin_inds,mlat_bin_inds = adata[:,0].astype(int),adata[:,1].astype(int)
		self.b1a[mlt_bin_inds,mlat_bin_inds]=adata[:,2]
		self.b2a[mlt_bin_inds,mlat_bin_inds]=adata[:,3]

		self.b1p, self.b2p = np.zeros((nmlt,nmlat)), np.zeros((nmlt,nmlat))
		self.prob = np.zeros((nmlt,nmlat,ndF))
		self.b1p.fill(np.nan)
		self.b2p.fill(np.nan)
		self.prob.fill(np.nan)

		if fluxType in ['diff','mono','wave']:
			with open(self.pfile,'r') as f:
				f.readline() # pass header
				pdata_b = np.genfromtxt(f,max_rows=nmlt*nmlat) # 2 columns, b1 and b2
				pdata_p = np.genfromtxt(f,max_rows=nmlt*nmlat*ndF) # 1 column, pval

			#in the file the probability is stored with coupling strength bin
			#varying fastest (this is Fortran indexing order)
			pdata_p_column_dFbin = pdata_p.reshape((-1,ndF),order='F')

			#mlt is first dimension
			self.b1p[mlt_bin_inds,mlat_bin_inds]=pdata_b[:,0]
			self.b2p[mlt_bin_inds,mlat_bin_inds]=pdata_b[:,1]
			for idF in range(ndF):
				self.prob[mlt_bin_inds,mlat_bin_inds,idF]=pdata_p_column_dFbin[:,idF]
		else:
			print("Error: bad fluxType (%s)" % fluxType)


	def which_dF_bin(self,dF):
		"""

		Given a coupling strength value, finds the bin it falls into

		"""
		dFave = 4421. #Magic numbers!
		dFstep = dFave/8.
		i_dFbin = np.round(dF/dFstep)
		#Range check 0 <= i_dFbin <= n_dF_bins-1
		if i_dFbin < 0 or i_dFbin > self.n_dF_bins-1:
			i_dFbin = 0 if i_dFbin < 0 else self.n_dF_bins-1
		return int(i_dFbin)

	def prob_estimate(self,dF,i_mlt_bin,i_mlat_bin):
		"""

		Estimate probability of <something> by using tabulated
		linear regression coefficients ( from prob_b files )
		WRT coupling strength dF (which are different for each position bin)

		If p doesn't come out sensible by the initial regression,
		(i.e both regression coefficients are zero)
		then tries loading from the probability array. If the value
		in the array is zero, then estimates a value using adjacent
		coupling strength bins in the probability array.

		"""
		#Look up the regression coefficients
		b1,b2 = self.b1p[i_mlt_bin,i_mlat_bin],self.b2p[i_mlt_bin,i_mlat_bin]

		p = b1 + b2*dF #What is this the probability of?

		#range check 0<=p<=1
		if p < 0. or p > 1.:
			p = 1. if p > 1. else 0.

		if b1 == 0. and b2 == 0.:
			i_dFbin = self.which_dF_bin(dF)
			#Get the tabulated probability
			p = self.prob[i_mlt_bin,i_mlat_bin,i_dFbin]

			if p == 0.:
				#If no tabulated probability we must estimate by interpolating
				#between adjacent coupling strength bins
				i_dFbin_1 = i_dFbin - 1 if i_dFbin > 0 else i_dFbin+2 #one dF bin before by preference, two after in extremis
				i_dFbin_2 = i_dFbin + 1 if i_dFbin < self.n_dF_bins-1 else i_dFbin-2 #one dF bin after by preference, two before in extremis
				p = (self.prob[i_mlt_bin,i_mlat_bin,i_dFbin_1] + self.prob[i_mlt_bin,i_mlat_bin,i_dFbin_2])/2.

		return p

	def estimate_auroral_flux(self,dF,i_mlt_bin,i_mlat_bin):
		"""
		Does what it says on the tin,
		estimates the flux using the regression coeffecients in the 'a' files

		"""
		b1,b2 = self.b1a[i_mlt_bin,i_mlat_bin],self.b2a[i_mlt_bin,i_mlat_bin]
		p = self.prob_estimate(dF,i_mlt_bin,i_mlat_bin)
		#print p,b1,b2,dF
		flux = (b1+b2*dF)*p
		return self.correct_flux(flux)

	def correct_flux(self,flux):
		"""
		A series of magical (unexplained,unknown) corrections to flux given a particular
		type of flux
		"""

		if flux < 0.:
			flux = 0.

		if self.fluxType is not 'ions':
			#Electron Energy Flux
			if flux > 10.:
				flux = 0.5
			elif flux > 5.:
				flux = 5.
		else:
			#Ion Energy Flux
			if flux > 2.:
				flux = 2.
			elif flux > 4.:
				flux = 0.25

		return flux

	def get_gridded_flux(self, dF):
		"""
		Return the flux interpolated onto arbitary locations
		in mlats and mlts
		"""

		fluxgridN = np.zeros((self.n_mlat_bins // 2, self.n_mlt_bins))
		fluxgridN.fill(np.nan)
		#Make grid coordinates
		mlatgridN,mltgridN = np.meshgrid(self.mlats[self.n_mlat_bins // 2:],self.mlts,indexing='ij')

		fluxgridS = np.zeros((self.n_mlat_bins // 2,self.n_mlt_bins))
		fluxgridS.fill(np.nan)
		#Make grid coordinates
		mlatgridS,mltgridS = np.meshgrid(self.mlats[:self.n_mlat_bins // 2],self.mlts,indexing='ij')

		for i_mlt in range(self.n_mlt_bins):
			for j_mlat in range(self.n_mlat_bins // 2):
				#The mlat bins are orgainized like -50:-dlat:-90,50:dlat:90
				fluxgridN[j_mlat, i_mlt] = self.estimate_auroral_flux(dF, i_mlt, self.n_mlat_bins // 2 + j_mlat)
				fluxgridS[j_mlat, i_mlt] = self.estimate_auroral_flux(dF, i_mlt, j_mlat)

		# Interpolate flux linearly for each latitude ring in the wedge
		# of low coverage in northern hemisphere dawn/midnight region
		fluxgridN,inwedge = self.interp_wedge(mlatgridN, mltgridN, fluxgridN)
		self.inwedge = inwedge

		return mlatgridN, mltgridN, fluxgridN, mlatgridS, mltgridS, fluxgridS

	def interp_wedge(self,mlatgridN,mltgridN,fluxgridN):
		"""
		Interpolates across the wedge shaped data gap
		around 50 magnetic latitude and 23-4 MLT.
		Interpolation is performed individually
		across each magnetic latitude ring,
		only missing flux values are filled with the
		using the interpolant
		"""
		#Constants copied verbatim from IDL code
		x_mlt_min=-1.0   #minimum MLT for interpolation [hours] --change if desired
		x_mlt_max=4.0    #maximum MLT for interpolation [hours] --change if desired
		x_mlat_min=49.0  #minimum MLAT for interpolation [degrees]
		#x_mlat_max=67.0
		x_mlat_max=75.0  #maximum MLAT for interpolation [degrees] --change if desired (LMK increased this from 67->75)

		valid_interp_mlat_bins = np.logical_and(mlatgridN[:,0]>=x_mlat_min,mlatgridN[:,0]<=x_mlat_max).flatten()
		inwedge = np.zeros(fluxgridN.shape,dtype=bool) #Store where we did interpolation

		for i_mlat_bin in np.flatnonzero(valid_interp_mlat_bins).tolist():
			#Technically any row in the MLT grid would do, but for consistancy use the i_mlat_bin-th one
			this_mlat = mlatgridN[i_mlat_bin,0]
			this_mlt = mltgridN[i_mlat_bin,:]
			this_flux = fluxgridN[i_mlat_bin,:]

			#Change from 0-24 MLT to -12 to 12 MLT, so that there is no discontiunity at midnight
			#when we interpolate
			this_mlt[this_mlt>12.] = this_mlt[this_mlt>12.]-24.

			valid_interp_mlt_bins = np.logical_and(this_mlt>=x_mlt_min,this_mlt<=x_mlt_max).flatten()
			mlt_bins_missing_flux = np.logical_not(this_flux>0.).flatten()

			interp_bins_missing_flux = np.logical_and(valid_interp_mlt_bins,mlt_bins_missing_flux)

			inwedge[i_mlat_bin,:] = interp_bins_missing_flux

			if np.count_nonzero(interp_bins_missing_flux) > 0:

				#Bins right next to missing wedge probably have bad statistics, so
				#don't include them
				interp_bins_missing_flux_inds = np.flatnonzero(interp_bins_missing_flux)
				nedge=6
				for edge_offset in range(1,nedge+1):
					lower_edge_ind = interp_bins_missing_flux_inds[0]-edge_offset
					upper_edge_ind = np.mod(interp_bins_missing_flux_inds[-1]+edge_offset,len(interp_bins_missing_flux))
					interp_bins_missing_flux[lower_edge_ind] = interp_bins_missing_flux[interp_bins_missing_flux_inds[0]]
					interp_bins_missing_flux[upper_edge_ind] = interp_bins_missing_flux[interp_bins_missing_flux_inds[-1]]

				interp_source_bins = np.flatnonzero(np.logical_not(interp_bins_missing_flux))

				#flux_interp = interpolate.PchipInterpolator(this_mlt[interp_source_bins],this_flux[interp_source_bins])
				flux_interp = interpolate.interp1d(this_mlt[interp_source_bins],this_flux[interp_source_bins],kind='linear')
				fluxgridN[i_mlat_bin,interp_bins_missing_flux] = flux_interp(this_mlt[interp_bins_missing_flux])

		return fluxgridN,inwedge



