#!/bin/env python3

from random import random

from atoms import Cell
from system import System
from vectors import Vector

if __name__ == '__main__':
    n = 10
    dt = 0.000001
    size = Vector(0.1, 0.1, 0.1)
    simulation = System({Cell(
        (float(random() * size.x), float(random() * size.y),
         float(random() * size.z)), (0, 0, 0)) for _ in range(n)}, dt, size)

    simulation.next_period_step()

#     F = ForceProjections();
# http://www.fizika.unios.hr/rf/wp-content/uploads/sites/67/2011/02/CPwP.pdfc
# https://kitchingroup.cheme.cmu.edu/blog/2017/11/14/Forces-by-automatic-differentiation-in-molecular-simulation/
# https://wiki.fysik.dtu.dk/ase/_modules/ase/calculators/lj.html
# http://espressomd.org/html/tutorials_html/01-lennard_jones/01-lennard_jones.html
# https://pythoninchemistry.org/sim_and_scat/molecular_dynamics/intro.html
