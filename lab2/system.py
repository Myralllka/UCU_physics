from typing import Set

from atoms import Cell
from find_forces import find_force, find_velo, find_pos
from vectors import Vector


class System:
    def __init__(self, cells: Set[Cell], step: float, size: Vector):
        self.cells: Set[Cell] = cells
        self.step = step

    def update_coordinates(self, elem: Cell):
        """
        """
        elem.coordinates = find_pos(elem, self.step)

    def update_volosity(self, elem: Cell):
        """
        """
        elem.speed = find_velo(elem, self.step)

    def update_forces(self):
        """
        update all forces for cells
        """
        element: Cell
        element1: Cell
        element2: Cell
        for element in self.cells:
            element.force = Vector(.0, .0, .0)

        tmp_systems = self.cells.copy()
        to_delete = set()
        for element1 in self.cells:
            for element2 in self.cells:
                if element1 == element2:
                    continue
                if element2.distance(element1) <= element2.r_max:
                    force = find_force(element1, element2)
                    element1.force += force
                    element2.force += -force
                    to_delete.add(element2)
            tmp_systems -= to_delete
            to_delete.clear()

    def next_period_step(self):
        """
        Update the forces projections in cells
        """
        self.update_forces()
        element: Cell
        for element in self.cells:
            self.update_coordinates(element)
            self.update_volosity(element)
