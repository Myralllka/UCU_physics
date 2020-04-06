from atoms import Cell


def lj_potential_cells(cell_1: Cell, cell_2: Cell):
    """
    LJ potential calculation for two cells
    :return: value of the potential between two cells
    """
    r_12 = abs((cell_1.coordinates - cell_2.coordinates))
    return lj_potential(r_12, cell_1)


def lj_potential(r_12: float, cell: Cell):
    """
    LJ potential calculation from the cell on distance r_12 from it

    :param r_12: distance from the cell
    :param cell: the cell
    :return: LJ potential from the cell on distance r_12
    """
    tmp = (cell.sigma / r_12) ** 6
    return cell.epsilon * 4 * tmp * (tmp - 1)
