# License
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.

# Manuscript
The primary purpose in distributing this source code is to enable readers to reproduce the numerical results reported in the manuscript "Optimal Path Planning with Random Breakdowns" by Marissa Gee and Alexander Vladimirsky.


# Contributions
* The manuscript for this project was written and edited by Marissa Gee and Alexander Valdimirsky
* The problem statement was developed by Alexander Vladimirsky and Marissa Gee
* The algorithm design was performed by Marissa Gee and Alexander Vladimirsky
* The implementation is by Marissa Gee, and makes use of existing code developed by Elliot Cartee, Zachary Clawson, and Marc Gilles


# Instructions
## Requirements
The C++ Code requires two external libraries:
* [Boost](http://www.boost.org/), which is used for implementation of multidimensional arrays and heaps.
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page), which is used for solving linear systems.

Currently the C++ code is run through a Makefile that assumes the libraries are installed in the `/usr/local/include` directory.
If those libraries are installed elsewhere, you will have to modify the [Makefile](https://github.com/eikonal-equation/Random_Breakdowns/blob/main/Makefile) to make sure both libraries are properly linked.

The code uses the C++17 standard, and can be compiled using both gcc and icpc.

The plotting scripts are in [jupyter notebooks](https://jupyter.org/).
The Python code in the notebooks requires [numpy](https://numpy.org/) and [matplotlib](https://matplotlib.org/).

### Data
Example 4 relies on terrain data from Mars, found in the `data` directory. The speed and breakdown rates can be computed using the file [DataProcessing](https://github.com/eikonal-equation/Random_Breakdowns/blob/main/src/DataProcessing.m). The code assumes that the speed and breakdown rate data is stored in the `/usr/local/data` directory in `.csv` format and formatted as an `nx` by `ny` array. The region of interest was selected by hand and the data was exported using [JMARS](https://jmars.asu.edu/), a Mars GIS. The region used in Example 4 is approximately square, spanning latitudes  18.41 to 18.58 degrees North, and longitudes 77.07 to 77.46 degrees East.

## Running the Code
Assuming the libraries are appropriately linked, you should be able to compile the code and run the test cases from the manuscript.

The four examples from the manuscript can be reproduced with the following commands
* To compile & run Example 1: Convergence of iterative scheme (radially symmetric, constant environment, target and depot at located at the center of the domain):
` make run TEST=Example1 `
* To compile & run Example 2: Inhomogenous partial breakdown rate (inhomogenous partial breakdown rate, includes cases with and without total breakdowns, single target and depot located at the center of the domain):
` make run TEST=Example2 `
* To compile & run Example 3: Changing partial breakdown rate (includes multiple partial breakdown rates, no total breakdowns, three depots, one target, all at distinct locations):
` make run TEST=Example3 `
* To compile & run Example 4: Real-world terrain (speed and partial breakdown rate based on Martian terrain, no total breakdowns, single target and depot at the same location)):
` make run TEST=Example4 `


## Visualizing Output
You can visualize the value functions for each example using the jupyter notebook [RandomBreakdowns](https://github.com/eikonal-equation/Random_Breakdowns/blob/main/plotting/RandomBreakdowns.ipynb). In the file, you will need to specify the example that was run, the version of the model ('TwoBreakdownTypes', 'OnlyPartialBreakdowns', or 'OnlyTotalBreakdowns'), and the version of the iterative algorithm used ('V' for value iterations or 'VP' for value-policy iterations). All examples in the manuscript use value-policy iterations by default. For example, to visualize Example 1 of the paper, set `test = Example1`, `model = TwoBreakdownTypes`, and `mode = VP`.

### Visualizing optimal trajectories
For Examples 2, 3, and 4 optimal trajectories can be visualized using the jupyter notebook [OptimalTrajectories](https://github.com/eikonal-equation/Random_Breakdowns/blob/main/plotting/OptimalTrajectories.ipynb). The first section computes and plots optimal trajectories for Example 2 or Example 4 (with the example specified as outlined above). The second section computes and visualizes optimal trajectories for a range of partial breakdown rates for Example 3. The optimal trajectories are visualized over a single value function.

### Visualizing convergence results
For Example 1, the convergence under grid refinement can be visualized using the jupyter notebook [ConvergencePlotting](https://github.com/eikonal-equation/Random_Breakdowns/blob/main/plotting/ConvergencePlotting.ipynb). The notebook also plots 1-dimensional representations of the radially symmetric value functions, and 1-dimensional representations of the location dependent error between the numerical and analytic solutions.