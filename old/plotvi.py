from models import simple
from matplotlib import pyplot as plot
import numpy as np
import math as m


for vi in np.arange(75, 155, 5):
	ball = simple.Golfball()
	ball.setvi(vi)
	ball.print()

	x, y = ball.pos_data(0, ball.tof())
	plot.plot(x, y, label=format(vi, ".1f") + " m/s")

plot.legend()
plot.xlabel("x-position of ball (m)")
plot.ylabel("y-position of ball (m)")
plot.show()