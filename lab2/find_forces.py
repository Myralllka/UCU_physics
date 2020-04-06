#!/bin/env python

import numpy as np
from typing import List, Set
from sympy import symbols, diff
from math import sqrt
from atoms import Cell
from vectors import Vector


def projection(val: float, axes: float, vector: Vector) -> float:
    """
    Calculate the projection on
    :param val - absolute value for projection
    :param axes - the axes to project on it
    :param vector - vector of all coordinates x, y, z
    """
    return val * (
            axes / np.sqrt(vector.x ** 2 + vector.y ** 2 + vector.z ** 2))


def lj_force(r, epsilon, sigma):
    """
    Implementation of the Lennard-Jones potential 
    to calculate the force of the interaction.
    
    Parameters
    ----------
    r: float
        Distance between two particles (Å)
    epsilon: float 
        Potential energy at the equilibrium bond 
        length (eV)
    sigma: float 
        Distance at which the potential energy is 
        zero (Å)
    
    Returns
    -------
    float
        Force of the van der Waals interaction (eV/Å)
    """
    return 48 * epsilon * np.power(sigma, 12) / np.power(r, 13) \
           - 24 * epsilon * np.power(sigma, 6) / np.power(r, 7)


def find_pos(cell: Cell, dt: float):
    """
    Update the particle positions.
    """
    return cell.coordinates + cell.speed * dt + cell.force * (dt * dt * 0.5)


def find_velo(cell: Cell, dt: float):
    """
    Update the particle velocities.
    """
    return cell.speed + cell.force / cell.mass * dt


def find_force(cell_1: Cell, cell_2: Cell):
    """
    Find the projections of forces two cells
    :return: [Fx, Fy, Fz]
    """
    vr = cell_1 - cell_2  # vector from cell_1 to cell_2 
    r = abs(vr)  # absolute distanse between cells
    abs_force = lj_force(r, Cell.epsilon, Cell.sigma)
    return Vector(
            projection(abs_force, vr.x, vr),
            projection(abs_force, vr.y, vr),
            projection(abs_force, vr.z, vr)
            )
