from random import randrange, random
from spin_interaction import h
from math import exp


def monte_carlo(particles_mtrx, n, l):
    to_work_n = len(particles_mtrx) // n
    for _ in range(n):
        to_work_with = list()
        for __ in range(to_work_n):
            x = randrange(l)
            y = randrange(l)
            z = randrange(l)
            to_work_with.append(particles_mtrx[x][y][z])
        to_work_with *= 2

        for particle in to_work_with:
            e_old = h(particles_mtrx, x, y, z, particle)
            e_new = h(particles_mtrx, x, y, z, -1 * particle)
            delta_e = e_new - e_old
            if delta_e < 0:
                particles_mtrx[x][y][z] *= -1
            else:
                probability = exp(-1 * delta_e)
                real_prob = random()
                if real_prob <= probability:
                    particles_mtrx[x][y][z] *= -1


def get_M(particles_mtrx):
    n = len(particles_mtrx)
    m = 0
    for x in range(n):
        for y in range(n):
            for z in range(n):
                m += particles_mtrx[x][y][z]
    return m
