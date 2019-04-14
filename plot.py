
from models import drag2, simple2, lift3
from matplotlib import pyplot as plot
import numpy as np
from scipy.optimize import minimize

# parse command line arguments
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-li", "--loftinitial", type=float, default=40, help="Minimum loft in degrees")
parser.add_argument("-lf", "--loftfinal", type=float, default=50, help="Maximum loft in degrees ")
parser.add_argument("-st", "--stepsize", type=float, default=5, help="Step size")
parser.add_argument("-dt", "--dt", type=float, default=0.01, help="Time step - decrease this value if you see lines rather than curves")
parser.add_argument("-vi", "--velocity", type=float, default=50, help="Initial velocity to use")
parser.add_argument("-y0", "--height", type=float, default=0, help="Initial height to use")
parser.add_argument("-tb", "--timebonus", type=float, default=1, help="Useful for debugging, multiplies estimated time of flight")
parser.add_argument("-sp", "--spin", type=float, default=200, help="Spin parameter")
parser.add_argument("-cd", "--drag", type=float, default=0.3, help="Drag coefficient")
parser.add_argument("-cl", "--lift", type=float, default=0.02, help="Lift coefficient")

args = parser.parse_args()

assert args.loftinitial < args.loftfinal
assert args.stepsize <= (args.loftfinal - args.loftinitial)
assert args.stepsize != 0

# initial velocity = 50 m/s
initialVelocity = args.velocity

# extract components of velocity
def components_of(v, theta):
	return v * np.cos(theta), v * np.sin(theta)

# estimate time of flight (assumes no air resistance)
def est_tof(v, theta):
	vy = v * np.sin(theta)
	return args.timebonus * (2 * vy + np.sqrt(vy**2 + 2*simple2.g*args.height)) / simple2.g

# spin from launch angle
def spin_from_theta(spin, theta):
	return spin + theta/14

maxdist = 0
maxtheta = 0
withlift = False

# Plot for a range of loft angles
plot.figure(1)
for theta in np.arange(np.deg2rad(args.loftinitial), np.deg2rad(args.loftfinal), np.deg2rad(args.stepsize)):
	# Drag model
	ball = drag2.Golfball()
	ball.cd = args.drag

	# Set initial velocity
	vx, vy = components_of(initialVelocity, theta)
	ball.set_coords([0, args.height, vx, vy])

	# Solve eqn and plot
	time, res = ball.solve(0, est_tof(initialVelocity, theta), args.dt)
	x, y = res.T
	plot.plot(x, y, label=format(np.rad2deg(theta), ".1f")+" deg, w/ drag")

	# update maximum value
	if x[-1] >= maxdist:
		maxtheta = theta
		maxdist = x[-1]
		withlift = False

	# Lift model
	ball = lift3.Golfball()
	ball.cd = args.drag
	ball.cl = args.lift

	# Set initial velocity
	vx, vy = components_of(initialVelocity, theta)
	ball.set_coords([0, args.height, vx, vy])
	ball.set_spin(spin_from_theta(args.spin, theta))

	# Solve eqn and plot
	time, res = ball.solve(0, est_tof(initialVelocity, theta), args.dt)
	x, y = res.T
	plot.plot(x, y, label=format(np.rad2deg(theta), ".1f") + " deg, w/ drag and lift")

	# update maximum value
	if x[-1] >= maxdist:
		maxtheta = theta
		maxdist = x[-1]
		withlift = True


print("Maximum distance: " + format(maxdist, ".2f") + " m @ " + format(np.rad2deg(maxtheta), ".1f") + " deg")
print(withlift)

plot.grid(True)
plot.xlabel("x-position (m)")
plot.ylabel("y-position (m)")
plot.title("Ballistic trajectory of golf ball (v_i = " + format(initialVelocity, ".1f") + " m/s)")
plot.legend()
plot.show()