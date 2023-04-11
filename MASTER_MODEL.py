###########################################################################################
#################### SELF EXPLANATORY NAME ################################################
##################### TWO FUNCTIONS, FIRST MAKES THE FORCING, SECOND RUNS THE MODEL #######
###########################################################################################

import numpy as np
from netCDF4 import Dataset
from netCDF4 import MFDataset
import random as rand


def e_s(TK):
	
	T = TK - 273.15
	e_s = 6.11*10**(7.5*T/(237.5+T))  # TEMPERATURES IN CELSIUS!!!
	return(e_s) # Answer is in millibars!

###########################################################################
############ FORCING SET HAS F, P, e, SST, DT
############ PARAM SET HAS r, G_s
###########################################################################

def the_model(POINT_LON,POINT_LAT,YYY):

	FORCING_SET,PARAM_SET = Make_Forcing(POINT_LON,POINT_LAT,YYY)
	F_forcing = FORCING_SET[0]
	P_forcing = FORCING_SET[1]
	e_forcing = FORCING_SET[2]
	T_forcing = FORCING_SET[3]
	S_forcing = FORCING_SET[4]
	DT 	  = FORCING_SET[5]
	Nyears = FORCING_SET[6]
	N = len(F_forcing)

	Surf_P	= PARAM_SET[0]
	Jackson = PARAM_SET[1]
	v_c 	= PARAM_SET[2]
	r_v 	= PARAM_SET[3]
	r_s	= PARAM_SET[4]
	r_a 	= PARAM_SET[5]

	############ FOR EXPERIMENTS
	if np.isfinite(r_v) == True:
		r_v = np.linspace(r_v,1.43*r_v,N)
	else:
		r_v = np.ones(N)*np.inf
	##### OTHERWISE CONSTANT
#	r_v = np.ones(N)*r_v

	emis = 0.67*(e_forcing**0.08)
	q_forcing = e_forcing*0.622/(Surf_P - 0.37*e_forcing)

############ Geometry

	theta = 0.45
	h_s = 0.1			# surface layer depth [m]
	h_d = 0.9
	r = 1 - Jackson**(h_s*100) # Roots by Jackson

############### Physical Constants
	
	sig = 5.67e-08			# SB
	rho_a = 1.25			# density of air [kg/m^3]
	rho_l = 1000			# denisty of water [kg/m^3]
	rho_s = 1000			# density of dry soil [kg/m^3]
	c_ps  = 1000			# heat capacity of dry soil [J/kg/K]
	c_pa = 1003			# heat capacity of dry air [J/kg/K]
	c_pl  =  4182			# heat capacity of water [J/kg/K]
	L = 2257000			# Latent enthalpy of vaporization [J/kg]

################### Combined parameters

	mu_s = rho_l*h_s
	mu_d = rho_l*h_d		# storage capacity of root layer

################# State Variables

	T = np.zeros(N)
	m_s = np.zeros(N)
	m_d = np.zeros(N)
	LHF = np.zeros(N)
	TRANSP = np.zeros(N)
	VPD = np.zeros(N)

	T[0] = T_forcing[0]			# Initial condition for surface temperature
	m_s[0] = 0			# ""			surface moisture
	m_d[0] = 0			# ""			rood mosisture

	i = 0

	while i < N-1:		

		OLR = sig*(T[i]**4)
		DLR = sig*emis[i]*(T_forcing[i]**4)
		H = c_pa*rho_a*(T[i] - T_forcing[i])/r_a

		if S_forcing[i] > 0:
			Melt = S_forcing[i]
		if S_forcing[i] <= 0:
			Melt = 0

		es_sat = e_s(T[i])
		qs_sat = es_sat*0.622/(Surf_P - 0.37*es_sat)
		q_diff = qs_sat - q_forcing[i]			# Humidity difference near surface

		VPD[i] = es_sat - e_forcing[i] 

		if q_diff > 0:
			E_s 	= rho_a*m_s[i]*q_diff/(theta*r_s)		# Surface Transpiration
		else:
			E_s = 0

		if q_diff > 0:# and noTransp == False:

			beta_s = m_s[i]/theta
			beta_r = m_d[i]/theta

			T_s 	= rho_a*r*beta_s*q_diff/(r_v[i])			# Surface Evaporation [kg/m^2/s]
			T_r     = rho_a*(1-r)*beta_r*q_diff/(r_v[i])		# Root Transpiration  [kg/m^2/s]
		else:
			T_r = 0
			T_s = 0

		Prate = P_forcing[i]*rho_l/DT
		Mrate = Melt*rho_l/DT

		TRANSP[i] = T_r + T_s

		LHF[i] = L*(E_s + T_r + T_s)

		Cap = rho_l*v_c*(m_d[i] - m_s[i])

		C_eff = h_s*(c_ps*rho_s + c_pl*m_s[i]*rho_l)		# Effective heat capacity of storage [J/m^2/K]

		dT_dt = (F_forcing[i] - OLR + DLR - LHF[i] - H)/C_eff
		dms_dt = (Mrate + Prate + Cap - E_s - T_s)/mu_s
		dmd_dt = -(T_r + Cap)/mu_d

		T[i+1] 		= T[i] + dT_dt*DT
		m_s[i+1]	= m_s[i] + dms_dt*DT 
		m_d[i+1]	= m_d[i] + dmd_dt*DT

		if m_s[i+1] > m_d[i+1]:
			surface_ro = (m_s[i+1] - m_d[i+1])*h_s/h_d
			m_s[i+1] = m_d[i+1]
			m_d[i+1] = m_d[i+1] + surface_ro

		if m_d[i+1] > theta:
			m_d[i+1] = theta	
		if m_s[i+1] < 0:
			m_s[i+1] = 0
		if m_d[i+1] < 0:
			m_d[i+1] = 0
	
		i+=1

	mon_time = np.linspace(15,Nyears*365.25 - 15,Nyears*12)
	Time = np.linspace(1,Nyears*365.25,N)

######################################################33333

#	T_mon = np.interp(mon_time,Time,T)
	ms_mon = np.interp(mon_time,Time,m_s)
	VPD_mon = np.interp(mon_time,Time,VPD)

#	LHF_mon = np.interp(mon_time,Time,LHF)
#	md_mon = np.interp(mon_time,Time,m_d)

	return(VPD_mon,ms_mon,YYY)
#	return(T_mon,ms_mon,LHF_mon,md_mon,YYY,XXX)
#	return(ms_mon,T_mon,LHF_mon,YYY)

def Make_Forcing(POINT_LON,POINT_LAT,YYY):

	####################################
	####### NO LON VALUES LESS THAN 0
	####################################

	FORCING_SET = []
	PARAM_SET = []
	Nyears = 100

	Time = np.linspace(0,Nyears*12 -1 ,Nyears*12)

	data_path = '/n/huybers_lab/lzeppetello/RANDOM_DATA/'

#################### Radiation

	Rad_set = Dataset(data_path+'CERES_EBAF_Ed4.1_Subset_200003-202111.nc','r')
	Rad_lat = Rad_set['lat'][:]
	Rad_lon = Rad_set['lon'][:]

	A = abs(Rad_lat - POINT_LAT)
	B = abs(Rad_lon  - POINT_LON)

	Rlon = B.argmin()
	Rlat = A.argmin() 

	Rad_linear = Rad_set['sfc_net_sw_all_mon'][10:250,Rlat,Rlon]
	Rad = np.reshape(Rad_linear,(20,12))
	Rad_cycle = np.nanmean(Rad[:,:],axis=0)

############ CONTROL
#	Rad_forcing_mon = np.zeros(480)
#	Rad_forcing_mon[:240] = np.tile(Rad_cycle,20)
#	Rad_forcing_mon[240:] = Rad_linear
######### Constant
	Rad_forcing_mon = np.tile(Rad_cycle,Nyears)

########################## P forcing

	crupath = '/n/huybers_lab/lzeppetello/RANDOM_DATA/CRUv4/'
	pset = MFDataset(crupath+'cru*pre.dat.nc','r')
	p_lat = pset['lat'][:]
	p_lon = pset['lon'][:]
	P = pset['pre'][:,:,:] # in m/month
	if POINT_LON > 180:
		V_LON = POINT_LON - 360
	else:
		V_LON = POINT_LON
	A = abs(p_lat - POINT_LAT)
	B = abs(p_lon  - V_LON)
	Plon = B.argmin()
	Plat = A.argmin() 
############# CONTROL
#	P_forcing_mon = P[:,Plat,Plon]/1000.	
#	P_reshape = np.reshape(P[:,Plat,Plon]/1000.,(40,12))
#	P_cycle = np.nanmean(P_reshape,axis=0)
############# CONSTANT
#	P = np.reshape(P[:,Plat,Plon]/1000.,(40,12))
#	P_cycle = np.nanmean(P,axis=0)
#	P_forcing_mon = np.tile(P_forcing_cycle,Nyears)
############## STOCHASTIC
	P_forcing_mon = np.zeros(Nyears*12)
	P_reshape = np.reshape(P[:,Plat,Plon],(40,12))
	P_cycle = np.nanmean(P_reshape[:,:],axis=0)
	Trend_set = Dataset('Grass_Trends.nc','r')
	p_trend = Trend_set['trend'][YYY]*Nyears/10.
	FACT = np.linspace(1,1*(100 + p_trend)/100.,Nyears*12)

	j = 0
	while j < Nyears*12:
		mon = int(j%12)
		rand.seed()
#		P_mon = rand.gammavariate(P_cycle[mon]/10.,10)
		P_mon = rand.gammavariate(FACT[j]*P_cycle[mon]/10.,10)
		P_forcing_mon[j] = P_mon
		j+=1

	P_forcing_mon = P_forcing_mon/1000.

############## VAPOR
	big_vap_set = MFDataset(crupath+'cru*vap.dat.nc')
######## CONTROL
#	E_forcing_mon = big_vap_set['vap'][:,Plat,Plon]
######### CONSTANT	
#	E = np.reshape(big_vap_set['vap'][:,Plat,Plon],(40,12))
#	E_cycle = np.nanmean(E[:,:],axis=0)
#	E_forcing_mon = np.tile(E_cycle,Nyears)
######### WARMING
	E = np.reshape(big_vap_set['vap'][:,Plat,Plon],(40,12))
	E_cycle = np.nanmean(E[:,:],axis=0)
	F_warm = np.log10(E_cycle/6.11)
	T_o = np.nanmean((237.5*F_warm)/(7.5 - F_warm) + 273.15)
	T_rate = 4.5*Nyears/100.
	E_forcing_mon = np.tile(E_cycle,Nyears)
	EO = e_s(np.linspace(T_o,T_o + T_rate,Nyears*12))
	EO = EO - EO[0]
	E_forcing_mon = E_forcing_mon + EO
######### GET Dew Point Forcing
	fact_linear = np.log10(E_forcing_mon/6.11)
	T_forcing_mon = (237.5*fact_linear)/(7.5 - fact_linear) + 273.15
############## SNOW FORCING #######################################################
	snow_set = Dataset(data_path+'ERA5_global_snowdepth.nc')
	Slat = snow_set['latitude'][:]
	Slon = snow_set['longitude'][:]
	A = abs(Slat - POINT_LAT)
	B = abs(Slon - POINT_LON)
	slatdex = A.argmin()
	slondex = B.argmin()	
########### CONTROL
#	SD = snow_set['sd'][:,slatdex,slondex]
############# CONSTANT
	SD = snow_set['sd'][:,slatdex,slondex]
	Snow_raw = np.reshape(SD,(40,12))
	Snow_cycle = np.nanmean(Snow_raw[:,:],axis=0)
	SD = np.tile(Snow_cycle,Nyears)

############ REMOVING NON SUMMER HYDRO CYCLE
#	i = 0
#	while i < 40:
#		P_forcing_mon[i*12:(i*12)+5] = P_cycle[:5]
#		P_forcing_mon[(i*12)+5:(i*12)+8] = P_reshape[i,5:8]
#		P_forcing_mon[(i*12)+8:(i*12)+12] = P_cycle[8:]
#
#		SD[i*12:(i*12)+5] = Snow_cycle[:5]
#		SD[(i*12)+5:(i*12)+8] = Snow_raw[i,5:8]
#		SD[(i*12)+8:(i*12)+12] = Snow_cycle[8:]
#		i+=1
#######################################################################
############ INTERPOLATING TO DAYS ######################
######################################################################

	steps_per_day = 10
	DT = 86400/steps_per_day
	N = Nyears*365*steps_per_day

########### Adjusting monthly total rainfall to DT-ly total

	mon = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
	mon_lengths = np.tile(mon,Nyears)
	mon_time = np.linspace(15,Nyears*365.15-15,Nyears*12)
	Time = np.linspace(1,Nyears*365.15,N)
	Rain_rate = P_forcing_mon/(mon_lengths[:]*steps_per_day) # in m/DT

#####################################################################################

	Rad_forcing = np.interp(Time,mon_time,Rad_forcing_mon)
	T_forcing = np.interp(Time,mon_time,T_forcing_mon)
	e_forcing = np.interp(Time,mon_time,E_forcing_mon)
	P_forcing = np.interp(Time,mon_time,Rain_rate)
	SD_interp = np.interp(Time,np.linspace(30,365.25*Nyears,Nyears*12),SD)

	Melt = np.zeros(shape=SD_interp.shape)
	Melt[1:] = -(SD_interp[1:] - SD_interp[:-1]) # in m/DT
	S_forcing = Melt

####################################################################
################### FIGURING OUT PARAMETERS ########################
####################################################################

	biome_set = Dataset(data_path+'BIOMES.nc')
	blat = biome_set['latitude'][:]
	blon = biome_set['longitude'][:]
	biome = biome_set['layer'][:,:]
	np.putmask(biome,biome<1,np.nan)

	Indicies = ["SLC","SMC","SHC","TLC","TMC","THC",
		"SLD","SMD","SHD","TLD","TMD","THD","SLB","SMB","SHB","TLB","TMB","THB","SLN","SMN","SHN",
		"TLN","TMN","THN"]	

	if POINT_LON > 180:
		A = abs(blat - POINT_LAT)
		B = abs(blon - (POINT_LON-360))

	else:

		A = abs(blat - POINT_LAT)
		B = abs(blon - POINT_LON)

	Blon = B.argmin()
	Blat = A.argmin()
	B = biome[Blat,Blon]
	DEX = Indicies[int(B-1)]

	######## PRESSURE ###################################
	############# FROM ERA5 LAND ########################

	SP_set = Dataset('ERA5mean_surface_pressure.nc','r')
	SP_lon = SP_set['lon'][:]	 # min0
	SP_lat = SP_set['lat'][:]

	A = abs(SP_lat - POINT_LAT)
	B = abs(SP_lon - POINT_LON)

	latdex = A.argmin()
	londex = B.argmin()

	Pressure = SP_set['sp'][latdex,londex]

	if DEX[0] == "T":
		Jackson = 0.96 		# based on temperate forest
		v_c = 1e-06		# tuned to match CA forest 
		r_a = 200		# Raupach
		r_s = 1000		# Given from nevada case
		r_v = 300		# Tuned to match CA forest

	if DEX[1] == "H":
		r_v = 100

	if DEX[0] == "S":
		Jackson = 0.93  #Based on grasslands
		v_c = 1e-05	# Tuned to match South Dakota
		r_a = 750	# Raupach
		r_s = 1000	# Given from Nevada Case
		r_v = 52	# Tuned to match south dakota

	if DEX[0] == "S" and DEX[1] == "L": # DESERT
		Jackson = 0 	# No roots
		v_c = 5e-07	# tuned to match nevada desert case study
		r_v = np.inf	# No transpiration
		r_s = 125
		r_a = 1000

	PARAM_SET.append(Pressure) 
	PARAM_SET.append(Jackson)
	PARAM_SET.append(v_c)
	PARAM_SET.append(r_v)
	PARAM_SET.append(r_s)
	PARAM_SET.append(r_a)

################## FINAL #################################################

	FORCING_SET.append(Rad_forcing)
	FORCING_SET.append(P_forcing)
	FORCING_SET.append(e_forcing)
	FORCING_SET.append(T_forcing)
	FORCING_SET.append(S_forcing)
	FORCING_SET.append(DT)
	FORCING_SET.append(Nyears)

	return(FORCING_SET,PARAM_SET)


#import matplotlib.pyplot as plt
# DRY PLACE
#point_lat = 39.467
#point_lon = 360-117.617
## OK
#point_lat = 37
#point_lon = 360-97
# Wet place
#point_lat = 42.1
#point_lon = 360-75.4
#T,m,x = the_model(point_lon,point_lat,15)
#v = np.reshape(T,(100,12))
#vjja = np.nanmean(v[:,5:8],axis=1)
#m = np.reshape(m,(100,12))
#mjja = np.nanmean(m[:,5:8],axis=1)
#plt.scatter(mjja,vjja)
#plt.show()


