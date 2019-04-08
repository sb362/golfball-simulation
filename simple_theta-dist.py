import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-y", "--height", type=float, default=0.0, help="Initial height (y_0) in meters")
parser.add_argument("-m", "--mass", type=float, default=0.045, help="Initial mass in kilograms")
parser.add_argument("-g", "--gravity", type=float, default=9.81, help="g in m/s^2")

# https://hypertextbook.com/facts/2001/EmilyAccamando.shtml
# dt ~ 0.5 milliseconds, F ~ 9 kN

parser.add_argument("-f", "--force", type=float, default=9000,help="Force appled by driver in newtons")
parser.add_argument("-t", "--dt", type=float, default=0.0005, help="Contact time in seconds")

params = parser.parse_args()

import math
import numpy
import matplotlib.pylab as plot

# v_0 = F * dt / m_ball
vi = params.force * params.dt / params.mass


# d = (v*v / 2g) sin(2 theta) * ( 1 + sqrt(1 + 2g y_0 / (v sin(theta))^2) )
def calc_distance(angle):
	dist = vi / (2 * params.gravity) * numpy.sin(2 * angle)

	if params.height != 0:
		dist *= (1 + numpy.sqrt(1 + (2 * params.gravity * params.height) / (vi**2 * numpy.sin(2 * angle)**2)))
	else:
		dist *= 2

	return dist


theta_range = numpy.arange(math.pi / 6, math.pi / 3, math.pi / 36)

xres = []
yres = []
for theta in theta_range:
	d = calc_distance(theta)
	xres.append(theta)
	yres.append(d)
	print(theta, d)


plot.plot(xres, yres)
plot.show()
