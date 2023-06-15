import math
import csv

class Location:
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
    
class Fingerprint:
    def __init__(self, location: Location, sample: FingerprintSample) -> None:
        self.location = location
        self.sample = sample

class AccessPoint:
    def __init__(self, mac: str, loc: Location, p=20.0, a=5.0, f=2417000000):
        self.mac_address = mac
        self.location = loc
        self.output_power_dbm = p
        self.antenna_dbi = a
        self.output_frequency_hz = f

def read_database(database_path: str) -> "list[Fingerprint]":
    database = []
    with open(database_path, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        
        for row in reader:
            location = Location(float(row[0]), float(row[1]), float(row[2]))
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
            database.append(Fingerprint(location, fingerprint_sample))
        
    return database