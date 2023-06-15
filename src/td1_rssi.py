import math
import csv

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
    
class Fingerprint:
    def __init__(self, position: SimpleLocation, sample: FingerprintSample) -> None:
        self.position = position
        self.sample = sample
    
class FingerprintDatabase:
    def __init__(self) -> None:
        self.db = []

database = FingerprintDatabase()
with open('data/td1_data.csv', newline='') as csvfile:
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
with open("data/td1_result.csv", "w", newline='') as outputfile:
    writer = csv.writer(outputfile, delimiter=",")
    for fingerprint in database.db:
        row = [fingerprint.position.x, fingerprint.position.y, fingerprint.position.z, 0]
        for fingerprint_sample in fingerprint.sample.samples:
            row.append(fingerprint_sample.mac_address)
            row.append(round(fingerprint_sample.get_average_rssi(), 2))
        writer.writerow(row)

"""
print("(" + str(database.db[1].position.x) + ", " + str(database.db[1].position.y) + ", " + str(database.db[1].position.z) + ")")
for fingerprint_sample in database.db[1].sample.samples:
    print(fingerprint_sample.mac_address + " -> " + str(fingerprint_sample.get_average_rssi()))
"""