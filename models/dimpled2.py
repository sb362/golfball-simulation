import numpy as np
from models import lift2

g = lift2.g
drag = lift2.drag

angularvelocity = np.array([0, 200/60])


def lift(radius, cl, velocity, angular, density=1.225):
	return cl * 4.0/3 * (4 * np.pi**2 * density * radius**3 * np.cross(velocity, angular))


class Golfball(lift2.Golfball):
	def __init__(self):
		lift2.Golfball.__init__(self)

		self.cd = 0.15

	def accelerations(self):
		fg = np.array([0, -g * self.mass])
		fd = -drag(self.area(), self.cd, self.velocities())
		fl = lift(self.radius, self.cl, self.velocities(), angularvelocity)

		return (fg + fd + fl) / self.mass

