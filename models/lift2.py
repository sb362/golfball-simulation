import numpy as np
from models import drag2

g = drag2.g
drag = drag2.drag


def lift(radius, cl, velocity, density=1.225):
	return cl * 4.0/3 * (4 * np.pi**2 * density * radius**3 * velocity * 200.0/60)


class Golfball(drag2.Golfball):
	def __init__(self):
		drag2.Golfball.__init__(self)

		self.cl = 0.25

	def accelerations(self):
		fg = np.array([0, -g * self.mass])
		fd = -drag(self.area(), self.cd, self.velocities())
		fl = np.array([0, lift(self.radius, self.cl, self.velocities()[0])])

		return (fg + fd + fl) / self.mass

