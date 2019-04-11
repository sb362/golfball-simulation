from models import drag, simple
from matplotlib import pyplot as plot
import numpy as np
import math as m

for theta in np.arange(m.pi / 6, m.pi / 3, m.pi / 36):
    ball = drag.Golfball()
    ball.setloft(theta)
    x, y = ball.pos_data(0, ball.tof())
    plot.plot(x, y, label=format(np.rad2deg(theta), ".1f") + " deg")

plot.legend()
plot.xlabel("x-position of ball (m)")
plot.ylabel("y-position of ball (m)")
plot.show()