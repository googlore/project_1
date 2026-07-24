from physics.base import DynamicalSystem

class SHMOscillator(DynamicalSystem):

    def __init__(self, k=1.0, mass=1.0):
        super().__init__(mass)
        self.k = k

    def acceleration(self, x, v, t):
        return -(self.k/self.mass)*x

    def potential_energy(self, x):
        return 0.5*self.k*x*x

    @property
    def angular_frequency(self):
        return (self.k/self.mass)**0.5

    @property
    def period(self):
        from math import pi
        return 2*pi/self.angular_frequency
