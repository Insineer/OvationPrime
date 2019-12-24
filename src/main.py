import ovation_prime
from ovation_utilities import *
from numpy import ndarray, flip
from config import *


def process(forecastDatetimes, threshold):
    ConnectToDB()

    estimators = [ovation_prime.FluxEstimator(tp) for tp in ["diff", "mono", "wave"]]

    for forecastTime in forecastDatetimes:
        processDatetime(estimators, forecastTime, threshold)


def countBorder(fluxGrid: ndarray, threshold):
    HOURS = 24
    LATITUDES = 40
    resolution = list(map(lambda a, b: int(a // b), fluxGrid.shape, (LATITUDES, HOURS)))

    absThreshold = fluxGrid.max() * threshold

    border = []

    lastBound = None
    for mlt in range(HOURS):
        bound = None
        for mlat in range(LATITUDES):  # from south to north
            mlatFirstPixel = mlat * resolution[0]
            mltFirstPixel = mlt * resolution[1]
            avg = 0
            for dmlt in range(resolution[1]):
                for dmlat in range(resolution[0]):
                    avg += fluxGrid[mlatFirstPixel + dmlat][mltFirstPixel + dmlt]
            avg /= resolution[0] * resolution[1]
            if avg >= absThreshold:
                bound = mlat
                break
            if mlat == LATITUDES - 1:
                if bound is None:
                    bound = lastBound
        lastBound = bound
        border.append(bound)

    return border

def processDatetime(estimators, forecastTime, threshold):
    print("Forecast for " + str(forecastTime))

    currentTime = forecastTime - datetime.timedelta(minutes=55)

    fluxes = []
    for estimator in estimators:
        mlatGrid, mltGrid, fluxGrid = estimator.get_flux_for_time(currentTime, CalculateSolarWindCouplingFunc(currentTime))
        fluxes.append(fluxGrid)

    acceleratedFlux = (fluxes[1] + fluxes[2]) / 2
    diffFlux = fluxes[0]

    poleward = countBorder(flip(acceleratedFlux, 0), threshold)
    poleward = [90 - b for b in poleward]

    equatorward = countBorder(acceleratedFlux, threshold)
    equatorward = [b + 50 for b in equatorward]

    diff = countBorder(diffFlux, threshold)
    diff = [b + 50 for b in diff]

    PrintFluxesFields(fluxes, forecastTime)
    DrawPlot(acceleratedFlux, poleward, equatorward, forecastTime)

    poleward = ' '.join(str(a) for a in poleward)
    equatorward = ' '.join(str(a) for a in equatorward)
    diff = ' '.join(str(a) for a in diff)

    SaveToDatabase(forecastTime, currentTime, poleward, equatorward, diff)

def main():
    print("Inited")

    datetimes = []

    if datetimeToProcess is None:
        datetimes.append(GetNextForecastDatetime())
    else:
        datetimes.append(datetimeToProcess)

    process(datetimes, threshold)

if __name__ == "__main__":
    main()