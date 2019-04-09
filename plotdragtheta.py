from models import drag
from matplotlib import pyplot as plot
import numpy as np
import math as m


for theta in np.arange(m.pi / 6, m.pi / 3, m.pi / 36):
	ball = drag.Golfball()
	ball.setloft(theta)
	ball.print()

	x, y = ball.pos_data(0, ball.tof())
	plot.plot(x, y, label=format(np.rad2deg(theta), ".0f") + " degrees")

plot.legend()
plot.show()