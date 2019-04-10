import numpy as np
from scipy.integrate import odeint

g = 9.81

class Golfball:
	def __init__(self):
		self.x = 0
		self.y = 0
		self.vx = 0
		self.vy = 0

		# 45 g, 42.67 mm
		self.mass = 0.045
		self.radius = 0.04267

	def area(self):
		return 4 * np.pi * self.radius ** 2

	def print(self):
		print("mass: {:.1f} kg\nradius: {:.1f} m\narea: {:.1f} m^2\nx: {:.1f} m\ny: {:.1f} m\nvx: {:.1f} m/s\nvy: {:.1f} m/s\n".format(self.mass, self.radius, self.area(), self.x, self.y, self.vx, self.vy))

	def set_coords(self, coords):
		self.x, self.y, self.vx, self.vy = coords

	def coords(self):
		return np.array([self.x, self.y, self.vx, self.vy])

	def positions(self):
		return np.array([self.x, self.y])

	def velocities(self):
		return np.array([self.vx, self.vy])

	def accelerations(self):
		# [x accel, y accel]
		fg = np.array([0, -g])

		return fg / self.mass

	def derivatives(self):
		derivatives = np.zeros(4)
		derivatives[0:2] = self.velocities()
		derivatives[2:4] = self.accelerations()

		return derivatives

	def eqns(self, coords, t):
		self.set_coords(coords)
		if self.y < 0:
			return np.zeros_like(coords)
		else:
			return self.derivatives()

	def solve(self, t0, t1, dt=0.01):
		interval = np.linspace(t0, t1, (t1 - t0) / dt)
		return interval, odeint(self.eqns, self.coords(), interval)[:, :2]

