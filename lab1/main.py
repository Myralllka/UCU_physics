import random
from config import config
from monte_carlo import monte_carlo, get_M
from vizualization import visualization

L = 5
N = 10
N_max = 5

if __name__ == "__main__":
    matrix = config(L)
    Ms = list()
    for _ in range(N_max):
        monte_carlo(matrix, N, L)
        Ms.append(get_M(matrix))
    visualize_matrix(matrix)
