# Golf Ball Simulator 5000
## Physics 1B Group Discovery Project

### Requirements
#### Python
Only tested on Python 3/3.5 - might work on earlier versions
#### Packages
`numpy` for all the maths\
`scipy` for solving ODEs\
`matplotlib` for plotting

To install the above: \
`pip install --user numpy scipy matplotlib`\
or: \
`pip3 install --user numpy scipy matplotlib`

### Installation
#### Linux
Open a terminal
- Install Python with\
`apt-get install python3 pip3`
- Install the packages listed above

May also require:\
``apt-get install python3-tk``
#### Mac
Get Python here (3 or 3.5 should work) and follow instructions: https://www.python.org/downloads/mac-osx/

Open a terminal
- Get pip3 with `easy_install pip3` or `sudo easy_install pip3`
- Install the packages listed above

### Usage
Set current directory to wherever `plot.py` is with\
`cd <location>`

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

