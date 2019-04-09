from models import simple
from matplotlib import pyplot as plot
import numpy as np
import math as m


for theta in np.arange(m.pi / 6, m.pi / 3, m.pi / 36):
	ball = simple.Golfball()
	ball.setloft(np.deg2rad(theta))

	print(format(np.rad2deg(theta), ".1f") + " degrees, tof=", format(ball.tof(), ".4f"))
	x, y = ball.pos_data(0, ball.tof())
	plot.plot(x, y, label=format(np.rad2deg(theta), ".1f") + " degrees")

plot.legend()
plot.show()