from core import Fingerprint, Location


class NormHisto:
    def __init__(self, histogram: dict[int, float]):
        self.histogram = histogram


def histogram_matching(db: list[Fingerprint], sample: list[NormHisto]) -> Location:
  pass



def probability(histo1: NormHisto, histo2: NormHisto) -> float:
    common_rssi = set(histo1.histogram.keys()) & set(histo2.histogram.keys())
    similarity = 0.0
    for rssi in common_rssi:
        similarity += min(histo1.histogram[rssi], histo2.histogram[rssi])
    return similarity