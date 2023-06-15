


# load calibration database
from core import read_database
from SimpleFingerprint import simple_matching


calibration_database = read_database("data_files/data.csv")

# load test data (not filtered for histogram matching)
test_database = read_database("data_files/test_data_not_filtered.csv")


# Iterate over all fingerprints to test
for fingerprint_to_evaluate in test_database:
    RSSI_to_evaluate = fingerprint_to_evaluate.sample
    
    # Estimate location using both grid and least squares multilateration solvers
    location = simple_matching(calibration_database,RSSI_to_evaluate)
    
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
