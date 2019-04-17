import numpy as np
import argparse

from scipy.integrate import odeint as integrate
from matplotlib import pyplot as plot
from numpy.linalg import norm
from mpl_toolkits.mplot3d import Axes3D

parser = argparse.ArgumentParser()

# Ball parameters
constants = parser.add_argument_group("Constants")
constants.add_argument("-m", "--mass", default=0.045, help="Mass of ball (kg)")
constants.add_argument("-r", "--radius", default=0.04267/2, help="Radius of ball (m)")

constants.add_argument("-g", "--gravity", type=float, default=9.81, help="For when we get a Mars base (m/s/s)")
constants.add_argument("-d", "--density", type=float, default=1.225, help="Density of air (kg m^-3)")

constants.add_argument("-cd", "--drag", type=float, default=0.5, help="Coefficient of drag")
constants.add_argument("-cl", "--lift", type=float, default=0.22, help="Coefficient of lift")

# Initial parameters
initialparams = parser.add_argument_group("Initial parameters")
initialparams.add_argument("-vi", "--velocity", type=float, default=50, help="Initial velocity (m/s)")
initialparams.add_argument("-yi", "--height", type=float, default=0, help="Initial height (m)")

initialparams.add_argument("-sp", "--spin", type=float, default=0, help="Spin (z)")
initialparams.add_argument("-spy", "--spiny", type=float, default=0, help="Spin (y)")
initialparams.add_argument("-spx", "--spinx", type=float, default=0, help="Spin (x)")

# Loft angle
parser.add_argument("-loft", "--loft", type=float, default=30, help="Loft angle")

# Debugging
parser.add_argument("-v", "--verbose", action="store_true")

# Parse args
args = parser.parse_args()

g = args.gravity
density = args.density


# Drag equation
def drag(density, area, cd, velocity):
	return -0.5 * density * area * cd * norm(velocity) * velocity


# Lift equation
def lift(density, area, cl, velocity, rvelocity):
	S = 0.5 * density * area * cl

	# Get cross product of angular velocity & linear velocity and make it a unit vector
	vxr = np.cross(rvelocity, velocity)

	return S * norm(velocity)**2 / norm(vxr) * vxr


# Simple golfball, no drag, no lift, smooth
class BasicGolfball:
	def __init__(self):
		# Properties
		self.mass = args.mass
		self.radius = args.radius

		# Coordinates
		self.x = 0
		self.y = args.height
		self.z = 0

		self.vx = 0
		self.vy = 0
		self.vz = 0

		self.rvx = 0
		self.rvy = 0
		self.rvz = 0

	# Reference area, for a sphere this is the cross-section.
	def area(self):
		return np.pi * self.radius**2

	# Set initial velocity
	def set_velocity(self, v, theta):
		self.vx = v * np.cos(theta)
		self.vy = v * np.sin(theta)

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

		# Ensure golf ball doesn't dig through the Earth
		if self.y < 0:
			return np.zeros_like(coords)

		if args.verbose:
			print(t, self.velocity(), self.rvelocity(), self.acceleration(), self.racceleration())

		return self.differentials()

	# Solve for trajectory over given interval
	def solve(self, t0, t1, dt=0.01):
		interval = np.linspace(t0, t1, (t1 - t0) / dt)
		return integrate(self.__eqns, self.coords(), interval, tfirst=True)[:, :3]


class DragGolfball(BasicGolfball):
	def __init__(self):
		BasicGolfball.__init__(self)

	def cd(self):
		return args.drag

	def acceleration(self):
		fd = drag(density, self.area(), self.cd(), self.velocity())
		return BasicGolfball.acceleration(self) + fd / self.mass


class LiftGolfball(DragGolfball):
	def __init__(self):
		DragGolfball.__init__(self)

	def spinf(self):
		v = (self.velocity()**2).sum()**0.5
		w = (self.rvelocity()**2).sum()**0.5
		return w / v

	def cl(self):
		return args.lift * self.spinf()

	def acceleration(self):
		fl = lift(density, self.area(), self.cl(), self.velocity(), self.rvelocity())

		# Assumption: x-component of lift is zero.
		fl[0] = 0

		return DragGolfball.acceleration(self) + fl / self.mass

	def racceleration(self):
		return 0


theta = np.radians(args.loft)

# Ball without drag
ball = BasicGolfball()
ball.set_velocity(args.velocity, theta)

# Ball with drag
ball2 = DragGolfball()
ball2.set_velocity(args.velocity, theta)

# Ball with drag and lift
ball3 = LiftGolfball()
ball3.set_velocity(args.velocity, theta)
ball3.set_spin([args.spinx, args.spiny, args.spin])

# Solve
res = ball.solve(0, 10)
x, y, z = res.T

res2 = ball2.solve(0, 10)
x2, y2, z2 = res2.T

res3 = ball3.solve(0, 100)
x3, y3, z3 = res3.T

# Plot 2d
fig = plot.figure(figsize=plot.figaspect(2))

ax = fig.add_subplot(2, 1, 1)
ax.plot(x, y, label="w/o drag")
ax.plot(x2, y2, label="w/ drag")
ax.plot(x3, y3, label="w/ drag & lift")
ax.grid(True)
ax.set_xlabel("Distance (m)")
ax.set_ylabel("Height (m)")
ax.legend()

# Plot 3d
ax = fig.add_subplot(2, 1, 2, projection="3d")
ax.plot(x, z, y, label="w/o drag")
ax.plot(x2, z2, y2, label="w/ drag")
ax.plot(x3, z3, y3, label="w/ drag & lift")
ax.set_xlabel("Distance (m)")
ax.set_zlabel("Height (m)")
ax.legend()

plot.show()


