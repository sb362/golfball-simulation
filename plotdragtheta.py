from models import drag, simple
from matplotlib import pyplot as plot
import numpy as np
import math as m

ball = drag.Golfball()
ball.setloft(m.pi / 6)
ball.print()

x, y = ball.pos_data(0, ball.tof())
plot.plot(x, y, label=format(np.rad2deg(m.pi / 6), ".0f") + " degrees, w/ drag")

ball = simple.Golfball()
ball.setloft(m.pi / 6)
ball.print()

x, y = ball.pos_data(0, ball.tof())
plot.plot(x, y, label=format(np.rad2deg(m.pi / 6), ".0f") + " degrees, w/o drag")


plot.legend()
plot.xlabel("x-position of ball (m)")
plot.xlabel("y-position of ball (m)")
plot.show()