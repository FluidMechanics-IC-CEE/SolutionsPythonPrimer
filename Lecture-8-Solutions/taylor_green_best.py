import numpy as np
import matplotlib.pyplot as plt


def compute_taylor_green_theory(x: np.ndarray, y: np.ndarray, time: float, viscosity: float):
    """
    Compute the velocity field of the Taylor-Green vortex using theoretical expressions.

    Args:
        x (np.ndarray): X-coordinates (2D grid).
        y (np.ndarray): Y-coordinates (2D grid).
        time (float): Time at which to evaluate the velocity field (s).
        viscosity (float): Kinematic viscosity of the fluid (m²/s).

    Returns:
        tuple: u (x-velocity field), v (y-velocity field), vel_magnitude (velocity magnitude field), 
               vorticity (vorticity field), psi (stream function).
    """
    if viscosity <= 0:
        raise ValueError("Viscosity must be a positive value.")
    if time < 0:
        raise ValueError("Time must be non-negative.")

    # Velocity of Taylor Green system 
    u = np.sin(x) * np.cos(y) * np.exp(-2.0 * viscosity * time)
    v = - np.cos(x) * np.sin(y) * np.exp(-2.0 * viscosity * time)
    vel_magnitude = np.sqrt(u**2 + v**2)
    
    # Stream function
    psi = np.sin(x) * np.sin(y) * np.exp(-2.0 * viscosity * time)
    
    # Vorticity field of Taylor Green system 
    vorticity = 2.0 * np.sin(x) * np.sin(y) * np.exp(-2.0 * viscosity * time)

    return u, v, vel_magnitude, vorticity, psi


def compute_taylor_green_from_stream_stream(x: np.ndarray, y: np.ndarray, psi: np.ndarray):
    """
    Compute the velocity field from the stream function of the Taylor-Green vortex.

    Args:
        x (np.ndarray): X-coordinates (2D grid).
        y (np.ndarray): Y-coordinates (2D grid).
        psi (np.ndarray): Stream function defined on the 2D grid.

    Returns:
        tuple: u (x-velocity field), v (y-velocity field), 
               vel_magnitude (velocity magnitude field), vorticity (vorticity field).
    """

    # Velocity components from stream function
    u = np.gradient(psi, y[0,:], axis=1)  # del_psi/del_y
    v = -np.gradient(psi, x[:,0], axis=0) # -del_psi/del_x
    vel_magnitude = np.sqrt(u**2 + v**2)
    
    # Vorticity = del_v/del_x - del_u/del_y
    vorticity = np.gradient(v, x[:,0], axis=0) - np.gradient(u, y[0,:], axis=1)

    return u, v, vel_magnitude, vorticity


def compute_error(var_theory: np.ndarray, var: np.ndarray, var_name):
    """
    Compute the absolute error.

    Args:
        var_theory: Variable calculated from theoretical expressions.
        var: Variable calculated from numerical methods.
        var_name: Variable name as a string.

    Returns:
        tuple: Absolute differences between the fields.
    """
    
    error = np.abs(var_theory - var)
    print("Maximum error in " + var_name + f" computation: {np.max(error):.6e}")
    
    return error
    
    
def plot_field_with_quiver(x, y, u, v, F, title, climits, lquiver):
    """
    Plot the velocity field.

    Args:
        x, y: Coordinate grids.
        u, v: Velocity components.
        F: Field.
        title: Title for the plot.
        climits: A vetor of two elements containinf the color map limits
        lquiver: logicital switch to plot (or not plot) the quiver vectors
    """
    plt.figure(figsize=(8, 6))
    cm = plt.contourf(x, y, F, levels=20, cmap='jet',vmin=climits[0], vmax=climits[1])
    plt.colorbar(cm)
    if lquiver:
        plt.quiver(x, y, u, v, scale=40, pivot="middle", color="black")
    plt.title(title)
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.axis("equal")
    plt.show()


#%% Parameter setting
# Physical parameters
time = 1.0       # Time [s]
viscosity = 0.1  # Kinematic viscosity [m^2/s]

# Grid parameters
grid_resolution_x = 50  # Number of points along x axis
grid_resolution_y = 40  # Number of points along y axis
x_bounds = [0.0, 2.0 * np.pi]   # Domain bound along x direction [m]
y_bounds = [0.0, 2.0 * np.pi]   # Domain bound along y direction [m]

# Plotting parameters
clim_max = 0.9                         # maximum color map limit of the velocity contour plots
clim = np.array([-clim_max, clim_max]) # upper and lower limits of the color map of the velocity contour plots
clim_error = np.array([0.0, 0.003])    # upper and lower limits of the color map of the error contour plots
lquiver = True      # Logicital switch to plot (or not plot) the quiver vectors


#%% Main

# Grid  generation
x = np.linspace(x_bounds[0], x_bounds[1], grid_resolution_x)
y = np.linspace(y_bounds[0], y_bounds[1], grid_resolution_y)
x, y = np.meshgrid(x, y, indexing='ij')


# Compute velocity from theoretical expressions
( 
   u_theory, 
   v_theory, 
   vel_mag_theory, 
   vorticity_theory, 
   psi_theory 
) = compute_taylor_green_theory(x, y, time, viscosity)


# Compute velocity from stream function
u_psi, v_psi, vel_mag_psi, vorticity_psi = compute_taylor_green_from_stream_stream(x, y, psi_theory)


# Compute errors
u_error = compute_error(u_theory, u_psi, "u velocity")
v_error = compute_error(v_theory, v_psi, "v velocity")
vel_mag_error = compute_error(vel_mag_theory, vel_mag_psi, "velocity magnitude")
vorticity_error = compute_error(vorticity_theory, vorticity_psi, "vorticity")


# Plot fields
plot_field_with_quiver(x, y, u_theory, v_theory, u_theory, "Velocity Field u: theoretical [m/s]", clim, lquiver)
plot_field_with_quiver(x, y, u_psi, v_psi, u_psi, "Velocity Field u: from stream function [m/s]", clim, lquiver)
plot_field_with_quiver(x, y, u_theory, v_theory, u_error, "Error in computation of u [m/s]", clim_error, lquiver)

plot_field_with_quiver(x, y, u_theory, v_theory, v_theory, "Velocity Field v: theoretical [m/s]", clim, lquiver)
plot_field_with_quiver(x, y, u_psi, v_psi, v_psi, "Velocity Field v: from stream function [m/s]", clim, lquiver)
plot_field_with_quiver(x, y, u_theory, v_theory, v_error, "Error in computation of v [m/s]", clim_error, lquiver)

plot_field_with_quiver(x, y, u_theory, v_theory, vel_mag_theory, "Velocity magnitude: theoretical [m/s]", [0.0, clim_max], lquiver)
plot_field_with_quiver(x, y, u_psi, v_psi, vel_mag_psi, "Velocity magnitude: from stream function [m/s]", [0.0, clim_max], lquiver)
plot_field_with_quiver(x, y, u_theory, v_theory, vel_mag_error, "Error in computation of velocity magnitude [m/s]", clim_error, lquiver)

plot_field_with_quiver(x, y, u_theory, v_theory, vorticity_theory, "Vorticity field of the theoretical velocity [s^-1]", 2.0*clim, lquiver)
plot_field_with_quiver(x, y, u_psi, v_psi, vorticity_psi, "Vorticity field of the velocity obtained from stream function [s^-1]", 2.0*clim, lquiver)
plot_field_with_quiver(x, y, u_theory, v_theory, vorticity_error, "Error in computation of vorticity [s^-1]", [0.0, 0.07], lquiver)

plot_field_with_quiver(x, y, u_psi, v_psi, psi_theory, "Stream function [s^-1]", clim, lquiver)


