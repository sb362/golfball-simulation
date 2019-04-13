# Golf Ball Simulator 5000
## Physics 1B Group Discovery Project

### Requirements
#### Python
3.7
#### Packages
`numpy` for all the maths\
`scipy` for solving ODEs\
`matplotlib` for plotting

To install the above: \
`python3 -m pip install --user numpy scipy matplotlib`

### Installation
#### Linux
Open a terminal
- Install Python with\
`apt-get install python3`
- Install the packages listed above

May also require:\
``apt-get install python3-tk``

### Usage

Generate a plot with\
`python3 plot.py <optional parameters, --help for list of commands>`

### Code example
See `plot.py` for a more interesting plot
```python
# Import model & plotting stuff
from models import simple2
from matplotlib import pyplot as plot

# Create new golfball (without drag)
ball = simple2.Golfball()

# [x pos, y pos, x velocity, y velocity]
ball.set_coords([0, 0, 50, 30])

# Simulate from t = 0 to 10 sec
time, res = ball.solve(0, 10)
x, y = res.T
plot.plot(x, y)

# Render plot
plot.show()
```
If the output appears dodgy (i.e. lots of straight lines rather than curves) try decreasing `dt` like so\
`ball.solve(0, 10, dt=0.001)`(default: 0.01)

