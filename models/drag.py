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
		x_raw, y_raw = self.__pos_data(t0, t1)
		x_out = np.array([0])
		y_out = np.array([0])

		# make adjustments here

		return x_out, y_out

	# Time of flight (pos_y == 0)
	# TODO: scipy.optimize.fsolve could be useful here for complex models
	def tof(self):
		return 0
