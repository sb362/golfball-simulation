import numpy as np
from models import simple2

g = simple2.g


# F_d = 1/2 * air density * (flow velocity)^2 * coefficient of drag * surface area
def drag(area, cd, velocity, density=1.225):
	return 0.5 * density * velocity * np.abs(velocity) * cd * area


# F_l = 1/2 * air density * velocity * spin * coefficient of lift * area
def lift(area, cl, velocity, spin, density=1.225):
	return 0.5 * density * velocity * spin * cl * area


class Golfball(simple2.Golfball):
	def __init__(self):
		simple2.Golfball.__init__(self)

		self.cd = 0.3

		# rotational properties
		self.cl = 0.02
		self.spin = 50

	def set_spin(self, spin):
		self.spin = spin

	def accelerations(self):
		fg = np.array([0, -g * self.mass])
		fd = -drag(self.area(), self.cd, self.velocities())
		fl = np.array([0, lift(self.area(), self.cl, self.velocities()[0], self.spin)])

		return (fg + fd + fl) / self.mass

	def eqns(self, coords, t):
		self.set_coords(coords)

		if self.y < 0:
			return np.zeros_like(coords)
		else:
			return self.derivatives()
