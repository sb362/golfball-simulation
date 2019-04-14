
from models import drag2, simple2, lift2, lift3, dimpled2
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
parser.add_argument("-sp", "--spin", type=float, default=200, help="Spin parameter")

args = parser.parse_args()

assert args.loftinitial < args.loftfinal
assert args.stepsize <= (args.loftfinal - args.loftinitial)
assert args.stepsize != 0

# initial velocity = 50 m/s
initialVelocity = args.velocity

# extract components of velocity
def components_of(v, theta):
	return v * np.cos(theta), v * np.sin(theta)

# spin from launch angle
def spin_from_theta(spin, theta):
	return spin + theta/14

# estimate time of flight (assumes no air resistance / lift)
def est_tof(v, theta):
	vy = v * np.sin(theta)
	return (2 * vy + np.sqrt(vy**2 + 2*simple2.g*args.height)) / simple2.g

"""
# Plot for a range of loft angles
plot.figure(1)
for theta in np.arange(np.deg2rad(args.loftinitial), np.deg2rad(args.loftfinal), np.deg2rad(args.stepsize)):
	ball = lift3.Golfball()

	# Set initial velocity
	vx, vy = components_of(initialVelocity, theta)
	ball.set_coords([0, args.height, vx, vy])
	ball.set_spin(spin_from_theta(args.spin, theta))

	# Solve eqn and plot
	time, res = ball.solve(0, est_tof(initialVelocity, theta) * 30, args.dt)
	x, y = res.T
	plot.plot(x, y, label=format(np.rad2deg(theta), ".1f")+" deg, w/ lift")

	# Drag model
	ball = drag2.Golfball()

	# Set initial velocity
	vx, vy = components_of(initialVelocity, theta)
	ball.set_coords([0, args.height, vx, vy])

	# Solve eqn and plot
	time, res = ball.solve(0, est_tof(initialVelocity, theta), args.dt)
	x, y = res.T
	plot.plot(x, y, label=format(np.rad2deg(theta), ".1f") + " deg, w/o lift")


"""
for theta in np.arange(np.deg2rad(args.loftinitial), np.deg2rad(args.loftfinal), np.deg2rad(args.stepsize)):
	ball = lift3.Golfball()

	vx, vy = components_of(initialVelocity, theta)
	ball.set_coords([0, args.height, vx, vy])
	ball.set_spin(spin_from_theta(args.spin, theta))

	time, res = ball.solve(0, 10)
	x, y = res.T

	plot.plot(np.rad2deg(theta), [np.amax(x)], 'ro')
	print(theta, np.amax(x))


plot.grid(True)
plot.xlabel("loft angle (degrees)")
plot.ylabel("range (m)")
plot.title("v_i = " + format(initialVelocity, ".1f") + " m/s)")
plot.show()