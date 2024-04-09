# Quadcopter-Mission
Dynamics of Quadcopter based on mission 
# Quadcopter Trajectory and Control Command Package

This Python package calculates a feasible 3D trajectory and the control commands (thrust and torques) for a quadcopter, based on the specified mission inputs. The inputs include the initial and final positions and velocities, and the mission is assumed to be completed within 5 seconds.

## Inputs

Set the initial (Xi and Vi) and final (Xf and Vf) positions (X) and velocities (V) in `Main_SetInput.py` as shown below:

![code](https://github.com/vhemmati/Quadcopter-Mission/assets/93438814/db0bb244-8412-4f41-8e35-4c2dc53f4196)

## Outputs

The package generates the following outputs:

- **Feasible Trajectory**: A trajectory that respects mission limits.
- **Angular Position vs. Time**: Outputs are provided at each 0.1 second interval.
- **Thrust vs. Time**: Detailed thrust information over time.
- **Torque vs. Time**: Detailed torque information over time.

## Results Visualization

There are two scripts included for visualizing the results:

- `Plot_Path_Freez.py`: Shows the entire motion in one plot.
- `Plot_Path.py`: Provides an animation of the mission.

## How to Use

1. Download the code package.
2. Set the input values in `Main_SetInput.py`.
3. Run `Traj_Check.py` as the main library.
4. Execute `Main_SetInput.py` to start the trajectory and control command computation.

## Example Plots

[Comment: You should include the plots generated by the visualization scripts here. You can add images to your GitHub repository and link them in the README using relative paths or hosted URLs.]

For more detailed information on usage and contributions, refer to the specific documentation sections or the code comments within the package.
