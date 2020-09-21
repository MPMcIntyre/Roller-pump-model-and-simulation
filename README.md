# Roller-pump-model-and-simulation
A project dedicated to modelling the pulsatile flow and simulating the pressure response of roller-type peristaltic pumps with various numbers of rollers.

# Instructions for running simulation
Open "Main.m" to insert the operating conditions, pump dimentions and parameter values. Once all user options and variables are set to the desired values, run the "Main.m" file to populate the workspace. To simulate the pressure response, open the "Lumped_parameter_sim.slx" file and run it.
## User options:
* **plot_enabled == true/false:** Plots flow variables pertaining to the variable inputs for the simulation.
* **measured_roller_values == true/false:** Uses measured roller volume displacement values stored in "Measured_roller_volume_data.mat" with the use of "Measured_roller_flow_solver.m". If false, the "Main.m" program uses the integration approximation method with the pump dimentions provided by the user to approximate the roller induced flow.
* **s_time:** Total simulation time - linked to the Simulink model.
* **h:** Time step (increment) value - Linked to the ode15 solver step size of the Simulink model.
