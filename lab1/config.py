import random


def config(l):
    return [[[random.choice([-1, 1]) for _ in range(l)]
             for __ in range(l)] for ___ in range(l)]
