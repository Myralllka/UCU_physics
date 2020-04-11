def h(matrix, x, y, z, new_spin):
    def _get_val(x1, y1, z1):
        if (0 <= x1 < len(matrix) and
                0 <= y1 < len(matrix[0]) and
                0 <= z1 < len(matrix[0][0])):
            return matrix[x1][y1][z1]
        else:
            return 0

    res = 0
    res += _get_val(x + 1, y, z)
    res += _get_val(x - 1, y, z)
    res += _get_val(x, y + 1, z)
    res += _get_val(x, y - 1, z)
    res += _get_val(x, y, z + 1)
    res += _get_val(x, y, z - 1)

    return res * -new_spin
