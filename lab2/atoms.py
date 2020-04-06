from math import sqrt
from typing import Tuple
from vectors import Vector


class Cell:
    mass = 1
    sigma = 1  # radius of interaction
    epsilon = 1  # constant of interaction force

    __t = 1 + sqrt(1 - 4 * (- 1 / 400000)) / 2
    r_max = sigma / (__t * __t * __t * __t * __t * __t)

    def __init__(self, xyz: Tuple[float], sxyz: Tuple[float]):
        self.coordinates: Vector = Vector(*xyz)
        self.speed: Vector = Vector(*sxyz)
        self.force: Vector
        self.analytic_force: Vector

    def __copy__(self):
        return Cell(self.coordinates.copy(), self.speed.copy())

    def __eq__(self, other):
        return self.coordinates == other.coordinates

    def __hash__(self):
        return hash(self.coordinates)

    def distance(self, other):
        return abs(other.coordinates - self.coordinates)
