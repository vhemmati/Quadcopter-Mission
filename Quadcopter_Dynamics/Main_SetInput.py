import numpy as np
import Traj_Check

# Please, provide intial and final positions and velocities
vi= np.array([0, 0, 10])
vf= np.array([10, 0, 0])

xi= np.array([0.0, 0.0, 0.0])
xf= np.array([-10, 10, 0])

Traj_Check.dynamics(xi,vi,xf,vf)