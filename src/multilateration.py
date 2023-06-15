import numpy as np
from scipy.optimize import least_squares

from core import *

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
                    yield Location(xi, yi, zi)
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

def multilateration_grid(distances: "dict[str, float]", ap_locations: "dict[str, Location]", target_precision: int = 3) -> Location:
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
    for current_precision in range(target_precision+1): # Here we refine the grid to a precision of "target_precision"
        """
        print("Current precision: " + str(current_precision) + 
            " digits after decimal point (target: " + str(target_precision) + ")")
        """

        # Prepare grid
        step = 1/(10**current_precision)
        grid = Grid(
            Location(lower_bound_x, lower_bound_y, lower_bound_z), 
            Location(upper_bound_x, upper_bound_y, upper_bound_z)
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

def multilateration_least_squares(distances: "dict[str, float]", ap_locations: "dict[str, Location]") -> Location:
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
            location_obj = Location(location[0], location[1], location[2])
            return np.array([
                location_obj.dist_from_other(loc) - dist
                for dist, loc in zip(distances.values(), ap_locations.values())
            ])

        initial_guess = Location(0, 0, 0)  # Starting point for optimization

        # Optimize the objective function using least squares (Levenberg-Marquardt algorithm)
        result = least_squares(objective_function, [initial_guess.x, initial_guess.y, initial_guess.z])

        # Extract the optimized location
        optimized_location = result.x

        return Location(*optimized_location)
