import numpy as np
import sys
import models.simple as simple


# F_d = 1/2 * air density * (flow velocity)^2 * coefficient of drag * surface area
def drag(area, cd, flowvelocity, density=1.225):
	return 0.5 * density * (flowvelocity ** 2) * cd * area


class Golfball(simple.Golfball):
	def __init__(self):
		simple.Golfball.__init__(self)

		# drag coefficient http://scienceworld.wolfram.com/physics/DragCoefficient.html
		self._cd = 0.3

	def cd(self):
		return self._cd

	def setcd(self, cd):
		self._cd = cd

	def pos_data(self, t0, t1):
		x_out = []
		y_out = []
		dt = 0.01

		# displacement
		s_x = 0
		s_y = 0

		# velocity
		v_x = simple.horizontal(self.vi(), self.loft())
		v_y = simple.vertical(self.vi(), self.loft())

		# acceleration
		a_x = -drag(self.area() / 2, self.cd(), v_x) / self.mass()
		a_y = -drag(self.area() / 2, self.cd(), v_y) / self.mass() - simple.g

		# calculate everything over an interval (todo: figure out how to numpy-ify this)
		i = 0
		while s_y >= 0:
			s_x += dt * (v_x + a_x * dt)
			s_y += dt * (v_y + a_y * dt)

			v_x += a_x * dt
			v_y += a_y * dt

			a_x = -drag(self.area() / 2, self.cd(), v_x) / self.mass()
			a_y = -drag(self.area() / 2, self.cd(), v_y) / self.mass() - simple.g

			x_out.append(s_x)
			y_out.append(s_y)

		return x_out, y_out

	# Time of flight (pos_y == 0)
	# TODO: scipy.optimize.fsolve could be useful here for complex models
	#def tof(self):
	#	return 0
