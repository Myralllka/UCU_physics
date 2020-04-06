"""
==============
3D scatter-plot
==============

Demonstration of a basic scatter-plot in 3D.
"""

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np


def visualize_matrix(array):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(len(array)):
        for j in range(len(array)):
            for k in range(len(array)):
                c = 'r' if array[i][j][k] == -1 else 'b'
                m = 'o' if array[i][j][k] == -1 else '^'
                ax.scatter(i, j, k, c=c, marker=m)
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')
    plt.show()

def visualize_M():
    with open