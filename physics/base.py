from abc import ABC, abstractmethod
import numpy as np

class DynamicalSystem(ABC):
    """
    Abstract base class for any dynamical system.

    State vector:
        y = [position, velocity]
    """

    def __init__(self, mass=1.0):
        self.mass = mass

    @abstractmethod
    def acceleration(self, x, v, t):
        pass

    def derivatives(self, t, y):
        x, v = y
        return np.array([
            v,
            self.acceleration(x, v, t)
        ])

    def kinetic_energy(self, x, v):
        return 0.5 * self.mass * v**2

    @abstractmethod
    def potential_energy(self, x):
        pass

    def total_energy(self, x, v):
        return self.kinetic_energy(x, v) + self.potential_energy(x)from abc import ABC, abstractmethod
import numpy as np

class DynamicalSystem(ABC):
    """
    Abstract base class for any dynamical system.

    State vector:
        y = [position, velocity]
    """

    def __init__(self, mass=1.0):
        self.mass = mass

    @abstractmethod
    def acceleration(self, x, v, t):
        pass

    def derivatives(self, t, y):
        x, v = y
        return np.array([
            v,
            self.acceleration(x, v, t)
        ])

    def kinetic_energy(self, x, v):
        return 0.5 * self.mass * v**2

    @abstractmethod
    def potential_energy(self, x):
        pass

    def total_energy(self, x, v):
        return self.kinetic_energy(x, v) + self.potential_energy(x)
