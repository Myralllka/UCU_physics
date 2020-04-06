#!/bin/env python

import numpy as np

from atoms import Cell
from vectors import Vector


def projection(val: float, axes: float, vector: Vector) -> float:
    """
    Calculate the projection on
    :param val - absolute value for projection
    :param axes - the axes to project on it
    :param vector - vector of all coordinates x, y, z
    """
    return val * (axes / np.sqrt(vector.x * vector.x +
                                 vector.y * vector.y +
                                 vector.z * vector.z))


def lj_force(r, epsilon, sigma) -> float:
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


def find_pos(cell: Cell, dt: float) -> Vector:
    """
    Update the particle positions.
    """
    return cell.coordinates + cell.speed * dt + cell.force / cell.mass * (
            dt * dt * 0.5)


def find_velo(cell: Cell, dt: float) -> Vector:
    """
    Update the particle velocities.
    """
    return cell.speed + cell.force / cell.mass * dt


def find_force(cell_1: Cell, cell_2: Cell) -> Vector:
    """
    Find the projections of forces two cells
    :return: [Fx, Fy, Fz]
    """
    # vector from cell_1 to cell_2
    vr = cell_1.coordinates - cell_2.coordinates
    r = abs(vr)  # absolute distance between cells
    abs_force = lj_force(r, Cell.epsilon, Cell.sigma)
    return Vector(
        projection(abs_force, vr.x, vr),
        projection(abs_force, vr.y, vr),
        projection(abs_force, vr.z, vr)
    )
