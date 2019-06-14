#!/usr/bin/env python3

import numpy as np
import argparse

from scipy.integrate import odeint as integrate
from matplotlib import pyplot as plot
from numpy.linalg import norm
from mpl_toolkits.mplot3d import Axes3D

parser = argparse.ArgumentParser()

# Ball parameters
constants = parser.add_argument_group("Constants")
constants.add_argument("-m", "--mass", default=0.04593, help="Mass of ball (kg)")
constants.add_argument("-r", "--radius", default=0.04267/2, help="Radius of ball (m)")
constants.add_argument("-i", "--inertia", type=float, default=9.145e-6, help="Inertia of golf ball")

constants.add_argument("--clubmass", type=float, default=0.2, help="Mass of club head (kg)")

constants.add_argument("-g", "--gravity", type=float, default=9.81, help="For when we get a Mars base (m/s/s)")
constants.add_argument("-d", "--density", type=float, default=1.225, help="Density of air (kg m^-3)")
constants.add_argument("--viscosity", type=float, default=1.46e-5, help="Kinematic viscosity of air")

# Initial parameters
initialparams = parser.add_argument_group("Initial parameters")
initialparams.add_argument("-yi", "--height", type=float, default=0, help="Initial height (m)")
initialparams.add_argument("--vclub", type=float, default=51.4, help="Club speed (m/s)")

# Loft angle
loftangleparams = parser.add_argument_group("Loft angle parameters")
loftangleparams.add_argument("-li", "--loftinitial", type=float, default=10, help="Loft angle (initial)")
loftangleparams.add_argument("-lf", "--loftfinal", type=float, default=35, help="Loft angle (final)")
loftangleparams.add_argument("-st", "--step", type=float, default=5, help="Loft angle (step)")

# Debugging
parser.add_argument("-v", "--verbose", action="store_true")

# Parse args
args = parser.parse_args()

# Input validation
assert args.loftfinal > args.loftinitial, "Final loft angle must be gretaer than initial loft angle!"
assert args.step != 0, "Step must be non-zero!"
assert ((args.loftfinal - args.loftinitial) / args.step).is_integer(), "Step size must divide the change in loft angle!"

assert args.mass != 0, "Mass must be non-zero."
assert args.radius != 0, "Radius must be non-zero."
assert args.viscosity != 0, "Kinematic viscosity must be non-zero."
assert args.density != 0, "Density of air must be non-zero."

g = args.gravity
density = args.density

# Ball speed from club speed and loft angle
def ball_speed(theta):
	theta = np.radians(theta)
	e = 0.86 - 0.0029 * args.vclub * np.cos(theta)

	bfn = (1 + e) * args.vclub * np.cos(theta) / (1 + args.mass / args.clubmass)
	bfp = args.vclub * np.sin(theta) / (1 + args.mass / args.clubmass + (args.mass * args.radius**2 / args.inertia))
	return np.sqrt(bfn**2 + bfp**2)

# Spin
def ball_spin(theta):
	theta = np.radians(theta)
	bfp = args.vclub * np.sin(theta) / (1 + args.mass / args.clubmass + (args.mass * args.radius**2 / args.inertia))

	return args.mass * bfp * args.radius / args.inertia

# Coefficient of drag from Reynolds number, based on degree four polynomial.
def re_to_cd(re):
	# Clamp output value as it is only an approximation
	if re > 120000:
		return 0.370
	elif re < 53000:
		return 0.8

	# Array of coefficients
	coeffs = np.array([
			  9.46410458e-20, -3.80736984e-14,
			  5.72048806e-09, -3.81337408e-04,
			  9.92620188e+00
			])

	# Return value of polynomial approximation
	return np.polyval(coeffs, re)


# Linear velocity to Reynolds number (Re = velocity * diameter / k. viscosity)
def reynolds(velocity, radius):
	return 2 * radius * velocity / args.viscosity


# Linear velocity to drag coefficient
def sphere_cd(velocity, radius):
	cd = re_to_cd(reynolds(velocity, radius))
	return cd


# Drag equation
# F_d = 1/2 * air density * ref. area * coefficient * |velocity| * v
def drag(density, area, cd, velocity):
	return -0.5 * density * area * cd * norm(velocity) * velocity


# Lift equation
# F_l = 1/2 * air density * ref. area * coefficient * |v|^2 * (what x vhat)
def lift(density, area, cl, velocity, rvelocity):
	if cl == 0:
		return np.array([0, 0, 0])

	S = 0.5 * density * area * cl

	# Cross product of angular velocity and linear velocity, for direction of spin
	rxv = np.cross(rvelocity, velocity)
	rxv /= norm(rxv)

	# Magnitude of spin is considered in coefficient of lift
	return S * norm(velocity)**2 * rxv

# Simple golfball, no drag, no lift, smooth
class BasicGolfball:
	def __init__(self):
		# Properties
		self.mass = args.mass
		self.radius = args.radius

		# Position
		self.x = 0
		self.y = args.height
		self.z = 0

		# Velocity
		self.vx = 0
		self.vy = 0
		self.vz = 0

		# Rotational velocities
		self.rvx = 0
		self.rvy = 0
		self.rvz = 0

	# Reference area, for a sphere this is the cross-section.
	def area(self):
		return np.pi * self.radius**2

	# Set initial velocity
	def set_velocity(self, v, theta):
		self.vx = v * np.cos(np.radians(theta))
		self.vy = v * np.sin(np.radians(theta))

	# Set spin
	def set_spin(self, spin):
		self.rvx, self.rvy, self.rvz = spin

	# Get all coordinates
	def coords(self):
		return np.array([self.x, self.y, self.z, self.vx, self.vy, self.vz, self.rvx, self.rvy, self.rvz])

	# Set all coordinates [x, y, z, vx, vy, vz, rvx, rvy, rvz]
	def set_coords(self, coords):
		self.x, self.y, self.z, self.vx, self.vy, self.vz, self.rvx, self.rvy, self.rvz = coords

	# Returns numpy array of position coordinates
	def position(self):
		return np.array([self.x, self.y, self.z])

	# Returns numpy array of velocity at the current position
	def velocity(self):
		return np.array([self.vx, self.vy, self.vz])

	# Returns numpy array of acceleration at the current position
	def acceleration(self):
		return np.array([0, -g, 0])

	# Returns numpy array of rotational velocity (spin) at the current position
	def rvelocity(self):
		return np.array([self.rvx, self.rvy, self.rvz])

	# Returns numpy array of rotational acceleration at the current position
	def racceleration(self):
		return np.array([0, 0, 0])

	# Returns numpy array of differential eqns to be solved by odeint
	def differentials(self):
		d = np.zeros(9)

		d[0:3] = self.velocity()
		d[3:6] = self.acceleration()

		d[6:9] = self.racceleration()

		return d

	# (Internal) Updates coordinates and returns list of equations to solve (for odeint)
	def __eqns(self, t, coords):
		self.set_coords(coords)

		if args.verbose:
			print(t, self.velocity(), self.rvelocity(), self.acceleration(), self.racceleration())

		return self.differentials()

	# Solve for trajectory over given interval
	def solve(self, t0, t1, dt=0.01):
		interval = np.linspace(t0, t1, int((t1 - t0) / dt))
		res = integrate(self.__eqns, self.coords(), interval, tfirst=True)[:, :3]

		out = np.array([e for e in res if e[1] >= 0])
		return out

# Simple golf ball but with drag
class DragGolfball(BasicGolfball):
	def __init__(self):
		BasicGolfball.__init__(self)

	# Coefficient of drag from velocity & radius
	def cd(self):
		return sphere_cd(norm(self.velocity()), self.radius)

	def acceleration(self):
		fd = drag(density, self.area(), self.cd(), self.velocity())
		return BasicGolfball.acceleration(self) + fd / self.mass


# Golfball with lift and drag
class LiftGolfball(DragGolfball):
	def __init__(self):
		DragGolfball.__init__(self)

	# Returns spin factor
	def spinf(self):
		v = norm(self.velocity())
		w = self.radius * norm(self.rvelocity())
		return w / v

	# Returns coefficient of lift based on spin factor
	def cl(self):
		s = self.spinf()
		return -3.25 * s**2 + 1.99 * s

	def acceleration(self):
		fl = lift(density, self.area(), self.cl(), self.velocity(), self.rvelocity())
		return DragGolfball.acceleration(self) + fl / self.mass

	# Spin decreases by about 1% every second
	def racceleration(self):
		return -0.01 * self.rvelocity()


# Initial conditions
density = 1.128
plot.figure()
for theta in np.arange(args.loftinitial, args.loftfinal, args.step):
	ball = LiftGolfball()
	ball.set_velocity(ball_speed(theta), theta)
	ball.set_spin([0, 0, ball_spin(theta)])

	res = ball.solve(0, 10)
	x, y, z = res.T

	plot.plot(x, y, label="Loft angle: " + format(theta, ".1f"))

plot.legend()
plot.grid(True)
plot.xlabel("Distance (m)")
plot.ylabel("Height (m)")
plot.title("Ballistic trajectory for air density " + format(density, ".3f") + " kg/m^3")

plot.figure()
xdata = []
ydata = []
for theta in np.arange(5, 30, 0.5):
	ball = LiftGolfball()
	ball.set_velocity(ball_speed(theta), theta)
	ball.set_spin([0, 0, ball_spin(theta)])

	res = ball.solve(0, 10)
	x, y, z = res.T

	xdata.append(theta)
	ydata.append(x[-1])

plot.plot(xdata, ydata, 'o', label="Air density: " + format(density, ".3f"))
plot.legend()
plot.grid(True)
plot.xlabel("Loft angle (m)")
plot.ylabel("Carry distance (m)")

plot.show()
