## Physics 1B Group Discovery Project

### Examples
![Alt text](images/results1128.png?raw=true "Trajectories for a variety of loft angles, default parameters (vclub = 51.4 m/s)")
![Alt text](images/results1128range.png?raw=true "Carry distance as a function of loft angle, same conditions/parameters")

### Files
`model.py` : deprecated 2D model, does not consider variable drag coefficient, etc.\
`model3d.py` : 3D model with drag based on Reynolds number; lift, backspin, etc.\
`images/*` : assorted images, mostly for the report\
`report/*` : a 2000-word report on the model (PDF and LaTeX documents)

### Requirements
#### Python
3.7 or 3.5 should work
#### Packages
`numpy` for all the maths\
`scipy` for solving ODEs\
`matplotlib` for plotting
