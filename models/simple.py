import numpy as np
import sys

# Gravitational acceleration assuming no air resistance, in m/s/s
g = 9.81


# printf equivalent
def printf(fmt, *args):
	sys.stdout.write(fmt % args)


# extract components from velocity
def horizontal(velocity, theta):
	return np.cos(theta) * velocity


def vertical(velocity, theta):
	return np.sin(theta) * velocity


# force, dt, mass of ball to initial velocity
# vi ~ F_avg * dt / mass
def force_to_vi(force, dt, mass):
	return force * dt / mass


class Golfball:
	def __init__(self):
		# Properties (in SI units)
		self._mass = 0.045
		self._radius = 0.0425 / 2

		# Launch properties (in SI units)
		self._vi = 100.0
		self._loft = np.pi / 6

		# Calculated in Golfball.update()
		self._area = 0.0

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

	def set_vi(self, force, dt):
		self._vi = force_to_vi(force, dt, self.mass())

	def area(self):
		return self._area

	def print(self):
		printf("\tMass: %.2f g\n\tRadius: %.2f mm\n\tLoft: %.2f deg\n\tv_0: %.2f m/s\n\tArea: %.2f mm^2\n\n",
			  self.mass() * 1000,
			  self.radius() * 1000,
			  self.loft(),
			  self.vi(),
			  self.area() * 1000
			   )

	# Updates properties like area and volume
	def update(self):
		self._area = 4 * np.pi * (self.radius() ** 2)

	# Positions at given time (seconds)
	def pos_x(self, time):
		return horizontal(self.vi(), self.loft()) * time

	def pos_y(self, time):
		return vertical(self.vi(), self.loft()) * time - (0.5 * g * (time ** 2))

	# Time of flight (pos_y == 0)
	# TODO: scipy.optimize.fsolve could be useful here for complex models
	def tof(self):
		return 2 * self.vi() * np.sin(self.loft()) / g

	# Returns position data over a given interval
	def pos_data(self, t0, t1):
		interval = np.linspace(t0, t1, 20)
		return self.pos_x(interval), self.pos_y(interval)
