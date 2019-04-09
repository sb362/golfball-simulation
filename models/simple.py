import numpy as np

# Gravitational acceleration assuming no air resistance, in m/s/s
g = 9.81


class Golfball:
	def __init__(self):
		# Properties (in SI units)
		self._mass = 0.045
		self._radius = 0.0425 / 2

		# Launch properties
		self._vi = 50.0
		self._loft = np.deg2rad(45.0)

		self.update()

	def mass(self):
		return self._mass

	def setmass(self, newmass):
		self._mass = newmass
		self.update()

	def radius(self):
		return self._radius

	def setradius(self, newradius):
		self._radius = newradius
		self.update()

	def vi(self):
		return self._vi

	def setvi(self, newvi):
		self._vi = newvi

	def loft(self):
		return self._loft

	def setloft(self, newloft):
		self._loft = newloft

	def area(self):
		return self._area

	def volume(self):
		return self._volume

	# Updates properties like area and volume
	def update(self):
		self._area = 4 * np.pi * (self.radius() ** 2)
		self._volume = 0.75 * np.pi * (self.radius() ** 3)

	# Positions at given time (seconds)
	def pos_x(self, time):
		return self.vi() * np.cos(self.loft()) * time

	def pos_y(self, time):
		return self.vi() * np.sin(self.loft()) * time - (0.5 * g * (time ** 2))

	# Time of flight (pos_y == 0)
	# TODO: scipy.optimize.fsolve could be useful here for complex models
	def tof(self):
		return 2 * self.vi() * np.sin(self.loft()) / g

	# Returns position data over a given interval
	def pos_data(self, t0, t1):
		interval = np.linspace(t0, t1, 20)
		return self.pos_x(interval), self.pos_y(interval)
