import argparse
import numpy as np

from scipy.integrate import solve_ivp
from matplotlib import pyplot as plot

parser = argparse.ArgumentParser()

# Ball parameters
parser.add_argument("-m", "--mass", default=0.045, help="Mass of ball (kg)")
parser.add_argument("-r", "--radius", default=0.04267, help="Radius of ball (m)")

parser.add_argument("-cd", "--drag", type=float, default=0, help="Coefficient of drag")
parser.add_argument("-cl", "--lift", type=float, default=0, help="Coefficient of drag")

parser.add_argument("-g", "--gravity", type=float, default=9.81, help="For when we get a Mars base (m/s/s)")
parser.add_argument("-d", "--density", type=float, default=1.225, help="Density of air (kg m^-3)")

# Initial parameters
parser.add_argument("-vi", "--velocity", type=float, default=50, help="Initial velocity (m/s)")
parser.add_argument("-yi", "--height", type=float, default=0, help="Initial height (m)")

parser.add_argument("-sp", "--spin", type=float, default=0, help="Spin")
parser.add_argument("--decay", type=float, default=0, help="Spin decay rate")

# Range of loft angles
parser.add_argument("-li", "--loftinitial", type=float, default=10, help="Minimum loft angle (degrees)")
parser.add_argument("-lf", "--loftfinal", type=float, default=20, help="Maximum loft angle (degrees)")
parser.add_argument("-st", "--step", type=float, default=2, help="Loft angle step (degrees)")

# Debugging/experimental
parser.add_argument("--samples", type=int, default=100, help="Number of time samples to use")
parser.add_argument("--continuous", action="store_true")
parser.add_argument("--odemethod", type=str, default="RK45")

# Parse arguments
args = parser.parse_args()
g = args.g

# Ensure arguments are within reality
assert args.loftfinal > args.loftinitial, "Final loft angle must be gretaer than initial loft angle!"
assert args.step != 0, "Step must be non-zero!"
assert ((args.loftfinal - args.loftinital) / args.step).is_integer(), "Step size must divide the change in loft angle!"

assert args.mass != 0, "Mass must be non-zero."

assert args.decay == 0, "Spin decay rate is not implemented."


# Drag equation (F_d = 1/2 * rho * ref area * coefficient of drag * v|v|)
def drag(area, cd, velocity, density):
	return 0.5 * density * area * cd * (velocity * np.abs(velocity))


# Lift equation
def lift():
	return 0


# Basic model (no drag, no lift, smooth)
class BasicGolfball:
	def __init__(self, loft):
		# Properties
		self.mass = args.mass
		self.radius = args.radius

		# Position & velocity values
		self.x = 0
		self.y = args.height
		self.vx = args.velocity * np.cos(np.deg2rad(loft))
		self.vy = args.velocity * np.sin(np.deg2rad(loft))

	# Reference area (for a sphere this is cross-sectional area)
	def area(self):
		return np.pi * self.radius**2

	# Set initial position & velocity [x, y, v_x, v_y]
	def set_coords(self, coords):
		self.x, self.y, self.vx, self.vy = coords

	def coords(self):
		return np.array([self.x, self.y, self.vx, self.vy])

	def position(self):
		return np.array([self.x, self.y])

	def velocity(self):
		return np.array([self.vx, self.vy])

	def acceleration(self):
		fg = np.array([0, -g * self.mass])

		return fg / self.mass

	__acceleration = acceleration

	# Returns system of differential equations to solve
	def differentials(self):
		derivatives = np.zeros(4)
		derivatives[0:2] = self.velocity()
		derivatives[2:4] = self.acceleration()

		return derivatives

	# Finds trajectory of golf ball over given interval
	# https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html
	def solve(self, t0=0, t1=10):
		interval = np.linspace(t0, t1, args.samples)
		return solve_ivp(self.__eqns, t_span=interval, y0=self.coords(), method=args.odemethod, dense_output=args.continuous)

	# Internal, do not call - returns set of eqns to be solved by SciPy
	def __eqns(self, t, coords):
		self.set_coords(coords)

		# Ensure golf ball doesn't dig through the Earth
		if self.y < 0:
			return np.zeros_like(coords)

		return self.differentials()


# Drag model (drag and no lift, smooth)
class DragGolfball(BasicGolfball):
	def __init__(self, loft):
		BasicGolfball.__init__(self, loft)

		self.cd = args.drag

	def acceleration(self):
		fd = -drag(self.area(), self.cd, self.velocity(), args.density)

		return BasicGolfball.acceleration(self) + fd / self.mass


# Lift model (drag, lift, smooth)
class LiftGolfball(DragGolfball):
	def __init__(self, loft):
		DragGolfball.__init__(self, loft)

		self.cl = args.lift
		self.spin = args.spin

	def acceleration(self):
		fl = lift()

		return DragGolfball.acceleration(self) + fl / self.mass



# Plot for a range of loft angles
for theta in np.arange(args.loftinitial, args.loftfinal, args.step):
	ball = LiftGolfball(theta)
	print(ball.solve(0, 10).status)



plot.grid(True)
plot.xlabel("Distance (m)")
plot.ylabel("Height (m)")
plot.title("Ballistic trajectory of golf ball with initial velocity " + format(args.velocity, ".0f") + " m/s")
plot.show()







