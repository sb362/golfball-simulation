from models import drag2, simple2
from matplotlib import pyplot as plot
import numpy as np

loftAngle = np.deg2rad(45)


def components_of(v, theta):
	return v * np.cos(theta), v * np.sin(theta)


def est_tof(v, theta):
	return 2 * v * np.sin(theta) / simple2.g


for vi in np.arange(20, 60, 10):
	ball = drag2.Golfball()

	vx, vy = components_of(vi, loftAngle)
	ball.set_coords([0, 0, vx, vy])

	time, res = ball.solve(0, est_tof(vi, loftAngle), 0.01)
	x, y = res.T
	plot.plot(x, y, label=format(vi, ".0f") + " m/s, with drag")

	ball.print()


plot.legend()
plot.show()
