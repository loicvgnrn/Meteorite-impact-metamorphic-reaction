# Meteorite-impact-metamorphic-reaction
2D Meteorite Impact Simulation with Improved Physics
This MATLAB project simulates a 2D meteorite impact event, including dynamic wave propagation with improved Perfectly Matched Layer (PML) absorption, heat conduction, and the mapping of metamorphic facies. The simulation also computes lithostatic pressure effects and provides evolving visualization and progress feedback via an enhanced waitbar.
Table of Contents
* Overview
* Features
* Installation & Prerequisites
* Usage
* Input Parameters
* Code Structure
* Visualization & Output
* Saving Results
* License
Overview
The code simulates the dynamic response of geological media to a meteorite impact. It combines multiple physics:
1. Dynamic Simulation: A wave propagation model using an improved split-field formulation with PML for absorbing boundaries.
2. Heat Conduction: A post-impact heat diffusion model where the initial temperature is perturbed by the impact.
3. Metamorphic Facies Mapping: Based on total pressure (dynamic + lithostatic) and temperature, the code classifies regions into different metamorphic facies.
A multi-stage waitbar informs the user of progress throughout the simulation, progressing from the dynamic simulation stage (0–60%), through the heat conduction stage (60–80%), and finally during metamorphic facies mapping (80–100%).
Features
* Dynamic Wavefield Simulation: Uses a Ricker wavelet source and spatially varying geological properties.
* Perfectly Matched Layer (PML): Implements a high-order PML to minimize boundary reflections.
* Heat Conduction Modeling: Simulates the diffusion of heat following the impact.
* Lithostatic Pressure Calculation: Computes cumulative pressure due to overburden weight.
* Metamorphic Facies Mapping: Classifies regions based on computed pressure and temperature.
* Interactive User Input: Prompts the user for simulation parameters such as impact velocity, layer properties, grid dimensions, etc.
* Visualization: Displays dynamic fields for displacement and pressure, and final plots for temperature, pressure, and metamorphic facies.
* Result Saving: Offers the option to save simulation outputs for further analysis.
Installation & Prerequisites
* MATLAB: Ensure MATLAB (R2016b or later is recommended) is installed.
* Toolboxes: No additional toolboxes are required. The script uses built-in MATLAB functions (e.g., conv2, inputdlg, waitbar, etc.).
To install or use the code, simply place the provided .m file (e.g., dep_march_improved_with_lithostatic_evolvedWaitbar.m) in your MATLAB working directory.
Usage
1. Launch MATLAB and set the working directory to where the code file is located.
2. Run the simulation by typing the following command in the MATLAB Command Window:
matlab
Copy
dep_march_improved_with_lithostatic_evolvedWaitbar
3. Input Parameters: A dialog box will appear requesting various simulation parameters (e.g., meteorite radius, impact velocity, geological layer properties, simulation grid dimensions, etc.). Default values are provided for convenience.
4. Progress Feedback: A waitbar will display the simulation’s progress through the three main stages:
   * Dynamic Simulation (0%–60%)
   * Heat Conduction (60%–80%)
   * Metamorphic Facies Mapping (80%–100%)
   5. Visualization: Real-time visualization is shown during the dynamic simulation stage. After the simulation, final plots for temperature, pressure, and metamorphic facies mapping are displayed.
   6. Saving Results: You will be prompted whether to save the results. If you choose "Yes," a dialog will let you select a filename and location for saving a MAT-file containing the simulation data.
Input Parameters
When prompted, the following parameters are required:
   * Meteorite Parameters:
   * Meteorite radius (m)
   * Impact velocity (m/s)
   * Geological Layers:
   * Number of layers
   * Thickness of each layer (comma-separated, in meters)
   * Densities of layers (comma-separated, in kg/m³)
   * Compressibilities (inverse of bulk modulus, comma-separated, in 1/Pa)
   * Thermal conductivities (comma-separated, in W/(m·K))
   * Simulation Domain & Grid:
   * Simulation width and depth (m)
   * Number of grid points in the x- and y-directions
   * Total simulation time (s)
   * PML thickness (m)
   * Impact Location:
   * Impact coordinates (m from left and top)
Additional parameters (e.g., Ricker wavelet frequency, PML order, and reflection coefficient) are set internally within the code.
Code Structure
The code is divided into several sections:
   1. User Input & Validation:
   * getSimulationParameters(): Prompts the user for simulation parameters.
   * validateParameters(): Checks that inputs are valid and consistent.
   2. Grid and Geological Properties Setup:
   * setupGrid(): Creates the spatial grid.
   * mapProperties(): Maps density, compressibility, and thermal conductivity onto the grid.
   3. Dynamic Simulation:
   * Time step calculation based on the CFL condition (calculateTimeStep()).
   * Wavefield initialization and source definition (defineSource(), injectSource()).
   * Main simulation loop using a split-field formulation with PML.
   4. Heat Conduction Modeling:
   * Computes initial temperature distribution using initialTemperatureMap().
   * Advances temperature using a finite-difference method.
   5. Pressure & Lithostatic Pressure Calculation:
   * Computes dynamic pressure from the displacement field.
   * lithostaticPressure(): Calculates the lithostatic pressure based on the density map.
   6. Metamorphic Facies Mapping & Final Visualization:
   * Processes final pressure and temperature maps to classify metamorphic facies.
   * Displays final results using finalizeVisualization().
   7. Helper Functions:
   * PML profile creation (buildPML_1D()), pressure computation, and visualization updates.
Visualization & Output
   * Dynamic Visualization: Displays two subplots in real time:
   * Displacement field (wave propagation)
   * Pressure field (computed from the displacement field)
   * Final Visualization: Presents three subplots showing:
   * Temperature distribution (°C)
   * Total pressure distribution (Pa)
   * Metamorphic facies map with overlaid labels
Saving Results
After final visualization, the script prompts the user with an option to save the results. If "Yes" is chosen, the temperature map, pressure map, metamorphic facies indices, and spatial grids are saved into a MATLAB .mat file.


Gpt o-1mini completion from original readme (20.02.25)
