import math
import csv
import numpy as np
from scipy.optimize import least_squares

class RSSISample:
    def __init__(self, mac_address: str, rssi: "list[float]") -> None:
        self.mac_address = mac_address
        self.rssi = rssi
    
    def get_average_rssi(self) -> float:
        rssi_milliwatts_sum = 0
        number_of_samples = 0
        for single_rssi in self.rssi:
            rssi_milliwatts_sum += 10 ** (single_rssi/10)
            number_of_samples += 1

        rssi_avg_milliwatts = rssi_milliwatts_sum / number_of_samples
        rssi_avg = 10 * math.log10(rssi_avg_milliwatts)
        return rssi_avg

class FingerprintSample:
    def __init__(self, samples: "list[RSSISample]") -> None:
        self.samples = samples
    
class SimpleLocation:
    def __init__(self, x: float, y: float, z: float) -> None:
        self.x = x
        self.y = y
        self.z = z

    def __eq__(self, other):
        x_eq = math.isclose(self.x, other.x, abs_tol=0.001)
        y_eq = math.isclose(self.y, other.y, abs_tol=0.001)
        z_eq = math.isclose(self.z, other.z, abs_tol=0.001)

        return (x_eq and y_eq and z_eq)
    
    def __str__(self):
        return f"({self.x}, {self.y}, {self.z})"
    
    def dist_from_other(self, other):
        delta_x_squared = math.pow(self.x - other.x, 2)
        delta_y_squared = math.pow(self.y - other.y, 2)
        delta_z_squared = math.pow(self.z - other.z, 2)
        distance = math.sqrt(delta_x_squared + delta_y_squared + delta_z_squared)
            
        return distance
    
class Fingerprint:
    def __init__(self, position: SimpleLocation, sample: FingerprintSample) -> None:
        self.position = position
        self.sample = sample
    
class FingerprintDatabase:
    def __init__(self) -> None:
        self.db = []

"""
database = FingerprintDatabase()
with open('data.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    
    while True:
        try:
            current_location = [next(reader), next(reader), next(reader), next(reader)]

        except StopIteration:
            break

        location = SimpleLocation(float(current_location[0][0]), float(current_location[0][1]), float(current_location[0][2]))
        samples = {}
        for row in current_location:
            row = row[4:] # Remove the first 4 elements of the row
            while row: # While row is not empty
                mac_address = row.pop(0)
                rssi = float(row.pop(0))

                if mac_address in samples:
                    samples[mac_address].append(rssi)
                else:
                    samples[mac_address] = [rssi]

        # Turn dictionary into a FingerprintSample object     
        rssi_samples = []
        for mac_address in samples.keys():
            rssi_samples.append(RSSISample(mac_address, samples[mac_address]))
        fingerprint_sample = FingerprintSample(rssi_samples)

        # Merge with location data to create a Fingerprint object, and add it to the database 
        database.db.append(Fingerprint(location, fingerprint_sample))

# Write results to an output csv file
with open("result.csv", "w", newline='') as outputfile:
    writer = csv.writer(outputfile, delimiter=",")
    for fingerprint in database.db:
        row = [fingerprint.position.x, fingerprint.position.y, fingerprint.position.z, 0]
        for fingerprint_sample in fingerprint.sample.samples:
            row.append(fingerprint_sample.mac_address)
            row.append(round(fingerprint_sample.get_average_rssi(), 2))
        writer.writerow(row)
"""

"""
print("(" + str(database.db[1].position.x) + ", " + str(database.db[1].position.y) + ", " + str(database.db[1].position.z) + ")")
for fingerprint_sample in database.db[1].sample.samples:
    print(fingerprint_sample.mac_address + " -> " + str(fingerprint_sample.get_average_rssi()))
"""

class AccessPoint:
    def __init__(self, mac: str, loc: SimpleLocation, p=20.0, a=5.0, f=2417000000):
        self.mac_address = mac
        self.location = loc
        self.output_power_dbm = p
        self.antenna_dbi = a
        self.output_frequency_hz = f

ap1 = AccessPoint("00:13:ce:95:e1:6f", SimpleLocation(4.93, 25.81, 3.55))
ap2 = AccessPoint("00:13:ce:95:de:7e", SimpleLocation(4.83, 10.88, 3.78))
ap3 = AccessPoint("00:13:ce:97:78:79", SimpleLocation(20.05, 28.31, 3.74))
ap4 = AccessPoint("00:13:ce:8f:77:43", SimpleLocation(4.13, 7.085, 0.80))
ap5 = AccessPoint("00:13:ce:8f:78:d9", SimpleLocation(5.74, 30.35, 2.04))
AP_list = {ap1.mac_address: ap1, ap2.mac_address: ap2, ap3.mac_address: ap3, ap4.mac_address: ap4, ap5.mac_address: ap5}

def compute_FBCM_index(distance: float, rssi_sample: RSSISample, ap: AccessPoint) -> float:
    """
    Computes a FBCM index based on the distance between transmitter and receiver,
    and the AP parameters. We consider the mobile device's antenna gain is 2.1 dBi.
    :param distance: the distance between AP and device
    :param rssi_values: the RSSI values associated to the AP for current calibration point. Use their average value.
    :return: one value for the FBCM index
    """

    # If we are too close to the access point, avoid division by zero
    if distance <= 0:
        distance = 0.0001
    
    Pr = (rssi_sample.get_average_rssi())
    Pt = (ap.output_power_dbm)
    Gr = (2.1)
    Gt = (ap.antenna_dbi)
    
    Lambda = 299792458/ap.output_frequency_hz
    
    # Friis formula computation
    index = (Pt - Pr + Gt + Gr + 20*math.log(Lambda/(4*math.pi), 10)) / (10*math.log(distance))
    return index

def estimate_distance(rssi_avg: float, fbcm_index: float, ap: AccessPoint) -> float:
    """
    Estimates the distance between an access point and a test point based on its rssi sample.
    :param rssi: average RSSI value for test point
    :param fbcm_index: index to use
    :param ap: access points parameters used in FBCM
    :return: the distance (meters)
    """

    Pr = (rssi_avg)
    Pt = (ap.output_power_dbm)
    Gr = (2.1)
    Gt = (ap.antenna_dbi)

    Lambda = 299792458/ap.output_frequency_hz
    
    # Friis formula computation
    estimated_dist = pow(10, (Pt - Pr + Gt + Gr + 20*math.log(Lambda/(4*math.pi))) / (10*fbcm_index))
    return estimated_dist

class Grid():
    """
    A volume defined by two points, where (xa < xb), (ya < yb) and (za < zb). 
    - lower_bound = coordinates of the lower bounding point A -> SimpleLocation
    - upper_bound = coordinates of the upper bounding point B -> SimpleLocation

    Methods
    - subpoints(self, step): Generator to get regularly-spaced points within the grid
    """

    def __init__(self, lower_bound, upper_bound):
        self.lower_bound = lower_bound
        self.upper_bound = upper_bound

    def subpoints(self, step):
        """
        Generator to get regularly-spaced points within the grid
        """
        # We use while loops instead of for loops because python's range() function doesn't support floats
        xi = self.lower_bound.x
        while xi <= self.upper_bound.x:
            yi = self.lower_bound.y
            while yi <= self.upper_bound.y:
                zi = self.lower_bound.z
                while zi <= self.upper_bound.z:
                    yield SimpleLocation(xi, yi, zi)
                    zi += step
                yi += step
            xi += step

def total_sphere_distance(location, distances, ap_locations):
    """
    For a given location, compute the sum of its distances to all AP spheres
    """
    dist_sum = 0
    for ap_mac in ap_locations:
        # If there is no distance for this access point, skip it
        if ap_mac not in distances:
            continue

        ap_location = ap_locations[ap_mac]
        ap_theoric_distance = distances[ap_mac]
        dist_from_ap = location.dist_from_other(ap_location)
        dist_from_sphere = abs(dist_from_ap - ap_theoric_distance)
        
        dist_sum += dist_from_sphere
    
    return dist_sum

def find_closest_location(distances, ap_locations, grid, step):
    """
    Iterate over all grid subpoints to find out which one is closest to all circles
    """
    closest_location = grid.lower_bound
    closest_location_dist = total_sphere_distance(grid.lower_bound, distances, ap_locations)
    for point in grid.subpoints(step):
        cur_location_dist = total_sphere_distance(point, distances, ap_locations)
        
        # If current point is closer than the previous closest, keep it instead of the old one
        if cur_location_dist < closest_location_dist:
            closest_location = point
            closest_location_dist = cur_location_dist
    
    return closest_location

def multilateration(distances: "dict[str, float]", ap_locations: "dict[str, SimpleLocation]") -> SimpleLocation:
    """
    Function multilateration computes a location based on its distances towards at least 3 access points
    :param distances: the distances associated to the related AP MAC addresses as a string
    :param ap_locations: the access points locations, indexed by AP MAC address as strings
    :return: a location
    """

    # Find out the grid boundaries by looking for the min and max coordinates amongst all APs
    lower_bound_x = min(map(lambda loc: loc.x, ap_locations.values()))
    lower_bound_y = min(map(lambda loc: loc.y, ap_locations.values()))
    lower_bound_z = min(map(lambda loc: loc.z, ap_locations.values()))
    upper_bound_x = max(map(lambda loc: loc.x, ap_locations.values()))
    upper_bound_y = max(map(lambda loc: loc.y, ap_locations.values()))
    upper_bound_z = max(map(lambda loc: loc.z, ap_locations.values()))

    step = 1.0
    target_precision = 3
    for current_precision in range(target_precision+1): # Here we refine the grid to a precision of "target_precision"
        print("Current precision: " + str(current_precision) + 
            " digits after decimal point (target: " + str(target_precision) + ")")
        
        # Prepare grid
        step = 1/(10**current_precision)
        grid = Grid(
            SimpleLocation(lower_bound_x, lower_bound_y, lower_bound_z), 
            SimpleLocation(upper_bound_x, upper_bound_y, upper_bound_z)
        )

        # Find closest point for the current grid
        closest_location = find_closest_location(distances, ap_locations, grid, step)

        # Round closest point to the current precision to avoid cumulative floating-point errors
        closest_location.x = round(closest_location.x, current_precision)
        closest_location.y = round(closest_location.y, current_precision)
        closest_location.z = round(closest_location.z, current_precision)
        """
        # Compute distance once more to display it to the user (this could be removed without affecting the algorithm)
        closest_location_distance = round(total_sphere_distance(closest_location, distances, ap_locations), target_precision)

        # Print current closest point for the user
        print("Closest point: ("
            + str(closest_location.x) + ", "
            + str(closest_location.y) + ", "
            + str(closest_location.z) + ") with a total distance of "
            + str(closest_location_distance) + "\n")
        """
        # Refine grid for next pass
        lower_bound_x = closest_location.x - step
        lower_bound_y = closest_location.y - step
        lower_bound_z = closest_location.z - step
        upper_bound_x = closest_location.x + step
        upper_bound_y = closest_location.y + step
        upper_bound_z = closest_location.z + step

    return closest_location

calibration_database = FingerprintDatabase()
with open('data/td1_result.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    
    for location_row in reader:
        location = SimpleLocation(float(location_row[0]), float(location_row[1]), float(location_row[2]))
        location_row = location_row[4:] # Remove the first 4 elements of the row
        
        # Parse row samples into an rssi samples list
        rssi_samples = []
        while location_row: # While row is not empty
            mac_address = location_row.pop(0)
            rssi = float(location_row.pop(0))
            rssi_samples.append(RSSISample(mac_address, [rssi]))

        # Use samples list to create a FingerprintSample object
        fingerprint_sample = FingerprintSample(rssi_samples)

        # Merge with location data to create a Fingerprint object, and add it to the database 
        calibration_database.db.append(Fingerprint(location, fingerprint_sample))

# list of fingerprint (used only for calibration)
fbcm_indexes = {}

# iterate over all access points
for ap in AP_list.values():

    # temporary list of indexes
    AP_fbcm_index = []

    # iterate over all fingerprints
    for fingerprint in calibration_database.db:

        # evaluate the real distance between the fingerprint and the AP
        distance = fingerprint.position.dist_from_other(ap.location)

        # run through the RSSI Sample list
        for rssiSample in fingerprint.sample.samples:
            
            #look for RSSI Samples that concern the AP (corresponding mac address)
            if rssiSample.mac_address == ap.mac_address:
                # compute index with this sample and knowing the distance. 
                # I firstly tried to filter long distances when rssi is high
                if rssiSample.rssi[0] > -80:
                    AP_POS_fbcm_index = compute_FBCM_index(distance, rssiSample, ap)
                    break

        # adding index to the list
        # I choose to filter indexes between 3 and 3.5
        if AP_POS_fbcm_index > 2 and AP_POS_fbcm_index < 3.5:
            AP_fbcm_index.append(AP_POS_fbcm_index)
    
    # averaging indexes
    fbcm_indexes[ap.mac_address] = sum(AP_fbcm_index)/len(AP_fbcm_index)

#print(fbcm_indexes)

"""
========================================================================================
"""

def evaluate(fingerprint_sample_to_evaluate, fbcm_indexes, AP_list, invp: int = 1) -> SimpleLocation:
    # distances with each access point
    AP_distances = {}

    # run through the AP list
    for ap in AP_list.values():
        # run through the list of RSSI values
        for rssi_sample in fingerprint_sample_to_evaluate.samples:
            
            # search for RSSI value that correspond to the AP
            if rssi_sample.mac_address == ap.mac_address:

                # estimate the distance using Friis index
                dist = estimate_distance(rssi_sample.rssi[0], fbcm_indexes[ap.mac_address], ap)

                # append to the list of distances
                AP_distances[ap.mac_address] = dist
                break

    #print(AP_distances)

    # extracts locations of APs
    locations_list = {}
    for ap in AP_list.values():
        locations_list[ap.mac_address] = ap.location

    # search for the location of the sample giving AP location and estimated distances with them
    pos = multilateration(AP_distances, locations_list)
    pos_better = beta_multilateration(AP_distances, locations_list)
    return (pos, pos_better)

"""
========================================================================================
"""

def beta_multilateration(distances: "dict[str, float]", ap_locations: "dict[str, SimpleLocation]") -> SimpleLocation:
        """
        Perform multilateration to compute the location based on distances and access point locations.
        :param distances: List of distances between the sample and each access point
        :param ap_locations: List of access point locations
        :return: Computed location as a SimpleLocation object
        """

        def objective_function(location: np.ndarray) -> np.ndarray:
            """
            Objective function for the optimization process.
            Calculates the differences between estimated distances and actual distances.
            :param location: Input location to evaluate as an np.ndarray
            :return: Differences between estimated distances and actual distances as an np.ndarray
            """
            location_obj = SimpleLocation(location[0].real, location[1].real, location[2].real)
            return np.array([
                location_obj.dist_from_other(loc) - dist
                for dist, loc in zip(distances.values(), ap_locations.values())
            ])

        initial_guess = SimpleLocation(0, 0, 0)  # Starting point for optimization

        # Optimize the objective function using least squares (Levenberg-Marquardt algorithm)
        result = least_squares(objective_function, [initial_guess.x, initial_guess.y, initial_guess.z])

        # Extract the optimized location
        optimized_location = result.x

        return SimpleLocation(*optimized_location)

"""
========================================================================================
"""

test_database = FingerprintDatabase()
with open('data/td2_test_data.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',')
    
    for row in reader:
        location = SimpleLocation(float(row[0]), float(row[1]), float(row[2]))
        row = row[4:] # Remove the first 4 elements of the row
        
        # Parse row samples into an rssi samples list
        rssi_samples = []
        while row: # While row is not empty
            mac_address = row.pop(0)
            rssi = float(row.pop(0))
            rssi_samples.append(RSSISample(mac_address, [rssi]))

        # Use samples list to create a FingerprintSample object
        fingerprint_sample = FingerprintSample(rssi_samples)

        # Merge with location data to create a Fingerprint object, and add it to the database 
        test_database.db.append(Fingerprint(location, fingerprint_sample))

sum = 0
sum_better = 0
N = 0
for fingerprint_to_evaluate in test_database.db:
    # real location
    RSSI_to_evaluate = fingerprint_to_evaluate.sample
    print('Real location: \t' + str(fingerprint_to_evaluate.position))
    
    # evaluate location
    (pos, pos_better) = evaluate(RSSI_to_evaluate, fbcm_indexes, AP_list, 2)
    print('Estimated location: \t' + str(pos))
    print('Estimated location better: \t' + str(pos_better))
    
    # compute distance between real location and the estimated one
    dist = fingerprint_to_evaluate.position.dist_from_other(pos)
    print('Error distance: \t' + str(dist))

    # compute distance between real location and the estimated one
    dist_better = fingerprint_to_evaluate.position.dist_from_other(pos_better)
    print('Error distance better: \t' + str(dist_better))

    N += 1
    sum += dist
    sum_better += dist_better
    #print()
    
print("Average: " + str(sum/N))
print("Average better: " + str(sum_better/N))

"""
#fingerprint database
_fingerprint_list = _db_calibration.fingerprint_list

#precompute FBCM indexes
_FBCM = FBCM(_AP_list,_fingerprint_list)

#load fingerprint to localize
_db_test = CSVParser("data_files/test_data.csv")

# run over the fingerprint list to evaluate
sum = 0;N = 0
for fingerprint_to_evaluate in _db_test.fingerprint_list:

    # real location
    RSSI_to_evaluate = fingerprint_to_evaluate._RSSISample_list
    print('Real location : \t' + str(fingerprint_to_evaluate.location))
    
    # evaluate location
    _pos = _FBCM.evaluate(RSSI_to_evaluate, 2)
    print('Estimated location : \t' + str(_pos))
    
    # compute distance between real location and the estimated one
    dist = evaluateDistance(_pos, fingerprint_to_evaluate.location)
    print('Error distance : \t' + str(dist))

    N += 1
    sum += dist
    #print()
    
print("Average : " + str(sum/N))
"""
