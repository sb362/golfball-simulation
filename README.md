# Models for Physics 1B group discovery project

### Requirements
##### Python
Only tested on Python 3/3.5 - might work on earlier versions
##### Packages
`numpy` for all the maths\
`scipy` for solving ODEs\
`matplotlib` for plotting

To install the above: \
`pip install --user numpy scipy matplotlib`\
or: \
`pip3 install --user numpy scipy matplotlib`

May also require:\
``apt-get install python3-tk``

### Usage
```python
# Import model & plotting stuff
from models import simple2
from matplotlib import pyplot as plot

# Create new golfball
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
If the output appears dodgy try decreasing `dt` like so\
`ball.solve(0, 10, dt=0.001)`(default: 0.01)

