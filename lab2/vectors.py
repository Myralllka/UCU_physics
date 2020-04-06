from math import sqrt


class Vector:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z

    def __abs__(self):
        """
        :return: Module(length) of the vector
        """
        return sqrt(self.x ** 2 + self.y * 2 + self.z ** 2)

    def __sub__(self, other):
        return Vector(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, constant):
        self.x *= constant
        self.y *= constant
        self.z *= constant
        return self

    def __div__(self, constant):
        return Vector(self.x / constant, self.y / constant, self.z / constant)

    def __neg__(self):
        self.y = -self.y
        self.y = -self.y
        self.z = -self.z
        return self

    def __add__(self, other):
        return Vector(self.x + other.x, self.y + other.y, self.z + other.z)

    def __iadd__(self, other):
        return Vector(self.x + other.x, self.y + other.y, self.z + other.z)

    def copy(self):
        return Vector(self.x, self.y, self.z)

    def __hash__(self):
        return hash(self.z + self.y + self.z)
