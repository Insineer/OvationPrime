import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import mysql.connector
import datetime

import config

conn = None

def ConnectToDB():
	try:
		global conn
		conn = mysql.connector.connect(host=config.host,
									   database=config.database,
									   user=config.user,
									   password=config.password)
		if conn.is_connected():
			print('Successfully connected to database!')
		else:
			print('Failed to connect database!')
			exit()
	except mysql.connector.Error as e:
		print(e)
		exit()


def PrintFluxesFields(fluxes, forecastTime):
	fluxesFileNamePrefix = forecastTime.strftime("%Y%m%d-%H:%M:%S")
	diffFileName = fluxesFileNamePrefix + "-diffusal.txt"
	monoFileName = fluxesFileNamePrefix + "-mono.txt"
	waveFileName = fluxesFileNamePrefix + "-wave.txt"

	np.savetxt(diffFileName, fluxes[0], fmt="%6.4f")
	np.savetxt(monoFileName, fluxes[1], fmt="%6.4f")
	np.savetxt(waveFileName, fluxes[2], fmt="%6.4f")

def DrawPlot(flux, poleward, equatorward, forecastTime):
	theta = np.linspace(0, 24, 96) / 24 * 2 * np.pi
	r = np.linspace(0, 40, 80)
	theta, r = np.meshgrid(theta, r)

	thetaForBorders = np.arange(0, 24, 1) / 24 * 2 * np.pi
	poleward = [90 - i for i in poleward]
	equatorward = [90 - i for i in equatorward]

	ax = plt.subplot(111, polar = True)

	ax.set_theta_zero_location("N")
	ax.set_theta_direction(-1)
	ax.set_rlabel_position(90)

	xticks = [i * 2 / 24 * 2 * np.pi for i in range(0, 12)]
	ax.set_xticks(xticks)
	ax.set_xticklabels([int(i / 2 / np.pi * 24) for i in xticks], color='black', fontsize=8)

	yticks = [0, 10, 20, 30]
	ax.set_yticks(yticks)
	ax.set_yticklabels([90 - i for i in yticks], color='#AAAAAA', fontsize=8)

	ctf = ax.pcolormesh(theta, r, flux, cmap = 'plasma')
	cbar = plt.colorbar(ctf, extend="max", label="aurorae power", shrink=0.5)
	cbar.minorticks_on()

	ax.plot(thetaForBorders, poleward)
	ax.plot(thetaForBorders, equatorward)

	ax.grid(True, color='#666666')

	fluxesFileNamePrefix = forecastTime.strftime("%Y%m%d-%H:%M:%S")
	plt.savefig(fluxesFileNamePrefix + '-acceleratedFluxes.png')

def SaveToDatabase(forecast_dt, processing_dt, poleward, equatorward, diffuse):
	print("Save to database:")
	print("Forecast time:   " + str(forecast_dt))
	print("Processing time: " + str(processing_dt))
	print("Poleward:        " + str(poleward))
	print("Equatorward:     " + str(equatorward))
	print("Diffuse:         " + str(diffuse))
	print("------")

	forecast_dt = forecast_dt.strftime("%Y-%m-%d %H:%M:%S")
	processing_dt = processing_dt.strftime("%Y-%m-%d %H:%M:%S")
	query = "REPLACE INTO `oval_forecast_ovation` (`forecast_dt`, `processing_dt`, `al_index`, `poleward`, `equatorward`, `diffuse`, `processed`) VALUES ('{}', '{}', '0', '{}', '{}', '{}', '0')".\
		format(forecast_dt, processing_dt, poleward, equatorward, diffuse)
	cursor = conn.cursor()
	cursor.execute(query)
	conn.commit()

	print("Successfully saved to DB!")


def __getMagneticData(dt):
	"""
		Return magnetic data for specified datetime

		input:
			dt - datetime
		output:
			bx_gsm, by_gsm, bz_gsm
	"""

	dt_start = dt - datetime.timedelta(hours=6)
	dt_end = dt

	dt_start = dt_start.strftime("%Y-%m-%d %H:%M:%S")
	dt_end = dt_end.strftime("%Y-%m-%d %H:%M:%S")

	cursor = conn.cursor()
	cursor.execute("SELECT * FROM `mag` WHERE `time_tag` BETWEEN ('{}') AND ('{}')".format(dt_start, dt_end))

	row = cursor.fetchone()
	mgt = [[], [], []]
	while row:
		mgt[0].append(row[1])
		mgt[1].append(row[2])
		mgt[2].append(row[3])
		row = cursor.fetchone()

	return sum(mgt[0]) / len(mgt[0]), sum(mgt[1]) / len(mgt[1]), sum(mgt[2]) / len(mgt[2])

def __getPlasmaData(dt):
	"""
		Return solar wind data for specified datetime

		input:
			dt - datetime
		output:
			density, speed
	"""
	dt_start = dt - datetime.timedelta(hours=6)
	dt_end = dt

	dt_start = dt_start.strftime("%Y-%m-%d %H:%M:%S")
	dt_end = dt_end.strftime("%Y-%m-%d %H:%M:%S")

	cursor = conn.cursor()
	cursor.execute("SELECT * FROM `plasma` WHERE `time_tag` BETWEEN ('{}') AND ('{}')".format(dt_start, dt_end))

	row = cursor.fetchone()
	speeds = []
	while row:
		if row[2] is not None:
			speeds.append(row[2])
		row = cursor.fetchone()

	return sum(speeds) / len(speeds)


def __couplingFunc(mag, V):
	"""
		Empirical coupling Func by Newell

		input:
			mag (x, y, z) - magnetic field
		output:
			V - solar wind speed
	"""
	Ec = np.zeros_like(mag[0])
	Ec.fill(np.nan)
	B = np.sqrt(mag[0] ** 2 + mag[1] ** 2 + mag[2] ** 2)
	BT = np.sqrt(mag[1] ** 2 + mag[2] ** 2)
	bztemp = mag[2]
	if bztemp == 0:
		bztemp = .001
	# Caculate clock angle (theta_c = t_c)
	tc = np.arctan2(mag[1], bztemp)
	neg_tc = BT * np.cos(tc) * mag[2] < 0
	if neg_tc:
		tc += np.pi
	sintc = np.abs(np.sin(tc / 2.))

	return (V ** 1.33333) * (sintc ** 2.66667) * (BT ** 0.66667)


def CalculateSolarWindCouplingFunc(dt):
	magX, magY, magZ = __getMagneticData(dt)
	plasmaSpeed = __getPlasmaData(dt)

	return __couplingFunc((magX, magY, magZ), plasmaSpeed)

def GetNextForecastDatetime():
	now = datetime.datetime.utcnow()
	previousHour = now.replace(minute=0, second=0, microsecond=0)
	return previousHour + datetime.timedelta(hours=1)
