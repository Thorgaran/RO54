from core import Location, Fingerprint


def rssi_distance(sample1: Fingerprint, sample2: Fingerprint) -> float:
    distance = 0.0
    for rssisample1 in sample1.sample.samples:
        mac_address1 = rssisample1.mac_address
        rssi1 = rssisample1.get_average_rssi()
        for rssisample2 in sample2.sample.samples:
            if rssisample2.mac_address == mac_address1:
                rssi2 = rssisample2.get_average_rssi()
                distance += (rssi1 - rssi2) ** 2
                break
        else:
            distance += rssi1 ** 2  # Treat missing MAC addresses as maximum distance
    for rssisample2 in sample2.sample.samples:
        mac_address2 = rssisample2.mac_address
        for rssisample1 in sample1.sample.samples:
            if rssisample1.mac_address == mac_address2:
                break
        else:
            rssi2 = rssisample2.get_average_rssi()
            distance += rssi2 ** 2  # Treat missing MAC addresses as maximum distance
    return math.sqrt(distance)


def simple_matching(db: list[Fingerprint], sample: Fingerprint) -> Location:
    min_distance = float('inf')
    closest_location = None
    for fingerprint in db:
        distance = rssi_distance(sample, fingerprint)
        if distance < min_distance:
            min_distance = distance
            closest_location = fingerprint.location
    return closest_location