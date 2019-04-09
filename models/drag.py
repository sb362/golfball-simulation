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

	# Calculates the effect of drag on x-y position over time interval dt
	# F = ma, a =
	def pos_x_adjusted(self, time, dt):
		return self.pos_x(time) - (drag(self.area(), self.cd(), 10) / self.mass()) * (dt**2)

	# Returns position data over a given interval
	def pos_data(self, t0, t1):
		interval = np.linspace(t0, t1, 20)
		return self.pos_x_adjusted(interval, 0), self.pos_y(interval)

