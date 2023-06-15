import csv

from core import *
from fbcm import FBCM

ap1 = AccessPoint("00:13:ce:95:e1:6f", Location(4.93, 25.81, 3.55))
ap2 = AccessPoint("00:13:ce:95:de:7e", Location(4.83, 10.88, 3.78))
ap3 = AccessPoint("00:13:ce:97:78:79", Location(20.05, 28.31, 3.74))
ap4 = AccessPoint("00:13:ce:8f:77:43", Location(4.13, 7.085, 0.80))
#ap5 = AccessPoint("00:13:ce:8f:78:d9", Location(5.74, 30.35, 2.04))
AP_list = [ap1, ap2, ap3, ap4]

calibration_database = read_database("data/td1_result.csv")
test_database = read_database("data/td2_test_data.csv")

fbcm = FBCM(AP_list, calibration_database)

sum_error_grid = 0
sum_error_least_squares = 0
nb_fingerprints = 0
# Iterate over all fingerprints to test
for fingerprint_to_evaluate in test_database:
    RSSI_to_evaluate = fingerprint_to_evaluate.sample
    
    # Estimate location using both grid and least squares multilateration solvers
    (location_grid, location_least_squares) = fbcm.estimate_fingerprint_location(RSSI_to_evaluate, AP_list, 6)
    
    # Compute error for both solvers
    error_grid = fingerprint_to_evaluate.location.dist_from_other(location_grid)
    error_least_squares = fingerprint_to_evaluate.location.dist_from_other(location_least_squares)
    
    # Verbose mode
    """
    print('Real location: \t\t\t\t' + str(fingerprint_to_evaluate.location))
    print('Estimated location (grid): \t\t' + str(location_grid))
    print('Estimated location (least squares): \t' + str(location_least_squares))
    print('Error (grid): \t\t\t\t%.3fm' % error_grid)
    print('Error (least squares): \t\t\t%.3fm' % error_least_squares)
    """

    nb_fingerprints += 1
    sum_error_grid += error_grid
    sum_error_least_squares += error_least_squares

    print("Fingerprint %d done" % nb_fingerprints)

print("\n%d fingerprints processed successfully" % nb_fingerprints)
print("Average error (grid): \t\t\t%.3fm" % (sum_error_grid/nb_fingerprints))
print("Average error (least squares): \t\t%.3fm" % (sum_error_least_squares/nb_fingerprints))
