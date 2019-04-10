
from models import drag2, simple2
from matplotlib import pyplot as plot
import numpy as np

initialVelocity = 50

def components_of(v, theta):
	return v * np.cos(theta), v * np.sin(theta)


def est_tof(v, theta):
	return 2 * v * np.sin(theta) / 9.81


for theta in np.arange(np.deg2rad(40), np.deg2rad(50), np.deg2rad(5)):
	ball = drag2.Golfball()

	vx, vy = components_of(initialVelocity, theta)
	ball.set_coords([0, 0, vx, vy])

	time, res = ball.solve(0, est_tof(initialVelocity, theta))
	x, y = res.T
	plot.plot(x, y, label=format(np.rad2deg(theta), ".1f")+" deg, w/ drag")

for theta in np.arange(np.deg2rad(40), np.deg2rad(50), np.deg2rad(5)):
	ball = simple2.Golfball()

	vx, vy = components_of(initialVelocity, theta)
	ball.set_coords([0, 0, vx, vy])

	time, res = ball.solve(0, est_tof(initialVelocity, theta))
	x, y = res.T
	plot.plot(x, y, label=format(np.rad2deg(theta), ".1f")+" deg, w/o drag")

plot.legend()
plot.show()
