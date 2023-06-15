from typing import Tuple

from core import *
from multilateration import multilateration_grid, multilateration_least_squares

class FBCM():
    def __init__(self, AP_list: "list[AccessPoint]", calibration_database: "list[Fingerprint]"):
        # Dictionary mapping mac addresses to friis indexes
        self.indexes = {}

        # Iterate over all access points
        for ap in AP_list:
            # List of all indexes computed for that AP
            all_AP_indexes = []

            # Iterate over all fingerprints
            for fingerprint in calibration_database:
                # Since this is the calibration data set, we do know 
                # the real distance between the fingerprint and the AP
                distance = fingerprint.location.dist_from_other(ap.location)

                # Iterate over all RSSI samples
                for rssiSample in fingerprint.sample.samples:
                    # Proceed only if the RSSI sample has the AP's mac address
                    if rssiSample.mac_address == ap.mac_address:
                        # Compute index if the RSSI isn't too weak
                        AP_POS_fbcm_index = -1
                        if rssiSample.rssi[0] > -80:
                            AP_POS_fbcm_index = FBCM.compute_index(distance, rssiSample, ap)

                        break
                
                # Add index to the list only when it's relevant
                if AP_POS_fbcm_index > 1.0 and AP_POS_fbcm_index < 3.0:
                    all_AP_indexes.append(AP_POS_fbcm_index)
            
            # Average all indexes and store the result in the dictionary
            self.indexes[ap.mac_address] = sum(all_AP_indexes)/len(all_AP_indexes)

    def estimate_fingerprint_location(
        self, fingerprint_sample_to_evaluate, AP_list, grid_target_precision
    ) -> "Tuple[Location, Location]":
        # Distances with each access point
        AP_distances = {}

        # Iterate over every access point
        for ap in AP_list:
            # Interate over every RSSI sample
            for rssi_sample in fingerprint_sample_to_evaluate.samples:
                # Proceed only if the RSSI sample has the AP's mac address 
                if rssi_sample.mac_address == ap.mac_address:

                    # Estimate the distance using Friis index and append it to the list
                    dist = self.estimate_distance(rssi_sample.rssi[0], ap)
                    AP_distances[ap.mac_address] = dist

                    break

        # Extract location of all access points
        locations_list = {}
        for ap in AP_list:
            locations_list[ap.mac_address] = ap.location

        # Solve location using two different multilateration algorithms
        location_grid = multilateration_grid(AP_distances, locations_list, grid_target_precision)
        location_least_squares = multilateration_least_squares(AP_distances, locations_list)
        return (location_grid, location_least_squares)

    def estimate_distance(self, rssi: float, ap: AccessPoint) -> float:
        """
        Estimates the distance between an access point and a test point based on its rssi sample.
        :param rssi: RSSI value for test point
        :param ap: access point corresponding to the measured RSSI value
        :return: the distance (meters)
        """

        Pr = (rssi)
        Pt = (ap.output_power_dbm)
        Gr = (2.1)
        Gt = (ap.antenna_dbi)

        Lambda = 299792458/ap.output_frequency_hz
        fbcm_index = self.indexes[ap.mac_address]

        # Friis formula computation
        estimated_dist = pow(10, (Pt - Pr + Gt + Gr + 20*math.log(Lambda/(4*math.pi))) / (10*fbcm_index))
        return estimated_dist

    def compute_index(distance: float, rssi_sample: RSSISample, ap: AccessPoint) -> float:
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
    