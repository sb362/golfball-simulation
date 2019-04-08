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


def calc_x(t, angle):
	return vi * numpy.cos(angle) * t


def calc_y(t, angle):
	return ((vi * t) * numpy.sin(angle)) - (0.5 * params.gravity * (time ** 2)) + params.height


for theta in numpy.arange(math.pi / 12, math.pi / 4, math.pi / 36):
	xres = []
	yres = []

	for time in numpy.linspace(0, 20, num=100):
		x = calc_x(time, theta)
		y = calc_y(time, theta)
		if y < 0:
			break

		xres.append(x)
		yres.append(y)

	# remove any y < 0
	p = [i for i, j in enumerate(yres) if j < 0]
	for i in sorted(p, reverse=True):
		del xres[i]
		del yres[i]

	plot.plot(xres, yres, label=format(numpy.rad2deg(theta), ".1f") + " degrees")

plot.legend()
plot.show()
