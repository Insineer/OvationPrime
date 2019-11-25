import ovation_prime
from ovation_utilities import *
import datetime

def getNextForecastDatetime():
    now = datetime.datetime.utcnow()
    previousHour = now.replace(minute=0, second=0, microsecond=0)
    return previousHour + datetime.timedelta(hours=1)

def process(forecastDatetimes):
    ConnectToDB()

    estimators = [ovation_prime.FluxEstimator(tp) for tp in ["diff", "mono", "wave"]]

    for forecastTime in forecastDatetimes:
        print("Forecast for " + str(forecastTime))

        currentTime = forecastTime - datetime.timedelta(minutes=55)

        fluxes = []
        for estimator in estimators:
            mlatGrid, mltGrid, fluxGrid = estimator.get_flux_for_time(currentTime)
            fluxes.append(fluxGrid)

        poleward = []
        f = open('bound.txt', 'w')
        f.write('mlt mlat bound avg \n')
        f1 = open('avg.dat', 'w')
        f1.write('mlt mlat 90-mlat mlat_grid_res mlt_grid_res avg \n')
        #AvgMas=[][]
        for mlt in range(0, 24):
            for mlat in range(0, 40):
                mlat_grid_res = (40 - mlat) * 2 - 1 # from north to south
                mlt_grid_res = mlt * 4
                avg = 0
                for dmlt in range(0, 4):
                    for dmlat in range(0, 2):
                        avg += fluxes[1][mlat_grid_res - dmlat][mlt_grid_res + dmlt]
                        avg += fluxes[2][mlat_grid_res - dmlat][mlt_grid_res + dmlt]
                avg /= 8 * 2
                print(mlt, mlat, 90-mlat, mlat_grid_res, mlt_grid_res, avg, file=f1)

        lst_bound = None
        for mlt in range(0, 24):
            bound = None
            for mlat in range(0, 40):
                mlat_grid_res = (40 - mlat) * 2 - 1 # from north to south
                mlt_grid_res = mlt * 4
                avg = 0
                for dmlt in range(0, 4):
                    for dmlat in range(0, 2):
                        avg += fluxes[1][mlat_grid_res - dmlat][mlt_grid_res + dmlt]
                        avg += fluxes[2][mlat_grid_res - dmlat][mlt_grid_res + dmlt]
                avg /= 8 * 2
    #            print(mlt, mlat, avg,file=f1)
                if bound is None and avg >= 0.25:
                    bound = 90 - mlat
                    #f.write(mlt, mlat, bound, avg,'\n')
                    #f = open('bound.txt', 'a')
                    print(mlt, mlat, bound, avg,file=f)
                    #f.close()
                    break
                if mlat == 39:
                    bound = lst_bound
                    #f = open('bound.txt', 'a')
                    print(mlt, mlat, bound, avg,file=f)
                    #f.close()
            lst_bound = bound
            poleward.append(bound)
            

        f.close()
        f1.close()
        equatorward = []
        diff = []
        eq_lst_bound = None
        diff_lst_bound = None
        for mlt in range(0, 24):
            eq_bound = None
            diff_bound = None
            for mlat in range(0, 40):  # from south to north
                mlat_grid_res = mlat * 2
                mlt_grid_res = mlt * 4
                avg_accelerated = 0
                avg_diff = 0
                for dmlt in range(0, 4):
                    for dmlat in range(0, 2):
                        avg_accelerated += fluxes[1][mlat_grid_res + dmlat][mlt_grid_res + dmlt]
                        avg_accelerated += fluxes[2][mlat_grid_res + dmlat][mlt_grid_res + dmlt]
                        avg_diff += fluxes[0][mlat_grid_res + dmlat][mlt_grid_res + dmlt]
                avg_accelerated /= 8 * 2
                avg_diff /= 8
                if eq_bound is None and avg_accelerated >= 0.25:
                    eq_bound = mlat + 50
                if diff_bound is None and avg_diff >= 0.3:
                    diff_bound = mlat + 50
                if mlat == 39:
                    if diff_bound is None:
                        diff_bound = diff_lst_bound
                    if eq_bound is None:
                        eq_bound = eq_lst_bound
            eq_lst_bound = eq_bound
            diff_lst_bound = diff_bound
            equatorward.append(eq_bound)
            diff.append(diff_bound)

        poleward = ' '.join(str(a) for a in poleward)
        equatorward = ' '.join(str(a) for a in equatorward)
        diff = ' '.join(str(a) for a in diff)

        PrintFluxesFields(fluxes, forecastTime)
        SaveToDatabase(forecastTime, currentTime, poleward, equatorward, diff)


def main():
    print("Inited")
    # datetimes = []
    # for day in range(1, 2):
    #     for hour in range (0, 1):
    #         datetimes.append(datetime.datetime(2019,4,day,hour,0,0))
    # process(datetimes)
    process([getNextForecastDatetime()])

if __name__ == "__main__":
    main()