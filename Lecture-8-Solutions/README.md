# Lecture 8

This directory contains files and resources related to **Lecture 8**. Below is a summary of its contents:

- `taylor_green_poorly_written.py`: A very poorly written code to compute the flow-field of a Taylor-Greeen vortex. This code has multiple syntax errors and numerical errors. One should never write a code like this.
- `taylor_green_no_errors.py`: After resolving all the errors in 'taylor_green_poorly_written.py' script.
- `taylor_green_best.py`: An example script that is written in a much better way following best practice.

The above python scripts demonstrate the common coding mistakes and the best practice coding with example codes that compute the Taylor-Green vortex fields and copares the numerically computed fields with repect to the values obtained from theoretical equations. This document highlights bad coding practices to avoid towards the end.

---

## Taylor-Green Vortex

The Taylor-Green vortex is an unsteady flow of decaying vortex with an exact closed-form analytical solution of the incompressible Navier-Stokes equations in Cartesian coordinates. This makes it an excellent benchmark for validating numerical models.

### Mathematical Formulation of Taylor-Green Vortex

- **Stream Function \($\psi$\)**:
```math
  \psi(x, y, t) = \sin(x) \sin(y) e^{-2\nu t} ,
```
   where, $x$ and $y$ are the Cartesian coordinates in the range $[0, 2\pi]$; $\nu$ is kinematic viscosity, and $t$ indicates time. Based on this stream function, the velocity fields can be derived as,

- **Velocity Fields \($u, v$\)**:
```math
  u = \frac{\partial \psi}{\partial y} = \sin(x) \cos(y) e^{-2\nu t},
  \quad \text{and,} \quad
  v = -\frac{\partial \psi}{\partial x} = -\cos(x) \sin(y) e^{-2\nu t}
```

- **Vorticity ($\omega$) from Velocity Fields**:
```math
  \omega = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y} = 2 \sin(x) \sin(y) e^{-2\nu t} = 2 \psi
```

---

## Coding Steps

1. **Generate 2D Grid**:
   - Generate a uniformly spaced 2D Cartesian grid $(x_{i,j}, y_{i,j})$, having $N_x$ and $N_y$ number of grid points along $x$- and $y$-directions respectively. Consider $N_x \ne N_y$, helps in debugging if there is a mistake with indexing. Use meshgrid.
   - Grid range: $x, y \in [0, 2\pi]$.

2. **Compute Theoretical Fields**:
   - Compute the velocity, stream function, and vorticity fields on the $(x_{i,j}, y_{i,j})$ grid using the theoretical expressions considering kinematic viscosity $\nu = 0.1 \, \text{m}^2/\text{s}$ at time $t = 1.0 \, \text{s}$.

3. **Numerical Differentiation**:
   - Compute the velocity fields ($u^∗$, $v^∗$) from the stream function by performing numerical differentiation:
```math
     u^*_{i,j} = \left. \frac{\partial \psi}{\partial y}\right\rvert_{i,j} = 
     \begin{cases}
      \frac{\psi_{i,j+1}-\psi_{i,j}}{\Delta y} \, , & j=0 \\
      \frac{\psi_{i,j+1}-\psi_{i,j-1}}{2\Delta y} \, , & 0 < j < N_y-1 \\
      \frac{\psi_{i,j}-\psi_{i,j-1}}{\Delta y} \, , & j=N_y-1
     \end{cases}
      \quad \text{and} \quad
     v^*_{i,j} = - \left. \frac{\partial \psi}{\partial x}\right\rvert_{i,j} = 
     \begin{cases}
      -\frac{\psi_{i+1,j}-\psi_{i,j}}{\Delta x} \, , & i=0 \\
      -\frac{\psi_{i+1,j}-\psi_{i-1,j}}{2\Delta x} \, , & 0 < i < N_x-1 \\
      -\frac{\psi_{i,j}-\psi_{i-1,j}}{\Delta x} \, , & i=N_x-1
     \end{cases}
```
   - Compute vorticity from velocity fields:
```math
     \omega^*_{i,j} = \left.\frac{\partial v^*}{\partial x}\right\rvert_{i,j} - \left.\frac{\partial u^*}{\partial y}\right\rvert_{i,j} \; , \text{ where ,}
```
```math
     \left.\frac{\partial v^*}{\partial x}\right\rvert_{i,j} = 
     \begin{cases}
      \frac{v^*_{i+1,j}-v^*_{i,j}}{\Delta x} \, , & i=0 \\
      \frac{v^*_{i+1,j}-v^*_{i-1,j}}{2\Delta x} \, , & 0 < i < N_x-1 \\
      \frac{v^*_{i,j}-v^*_{i-1,j}}{\Delta x} \, , & i=N_x-1
     \end{cases}
      \quad \text{and} \quad
     \left.\frac{\partial u^*}{\partial y}\right\rvert_{i,j} = 
     \begin{cases}
      \frac{u^*_{i,j+1}-u^*_{i,j}}{\Delta y} \, , & j=0 \\
      \frac{u^*_{i,j+1}-u^*_{i,j-1}}{2\Delta y} \, , & 0 < j < N_y-1 \\
      \frac{u^*_{i,j}-u^*_{i,j-1}}{\Delta y} \, , & j=N_y-1
     \end{cases}
```

4. **Error Analysis**:
   - Calculate numerical errors for each field compared to theoretical values. If $\Phi$ is any field variable then,
```math
     \Phi^{\text{error}}_{i,j} = \left| \phi^*_{i,j} - \phi_{i,j} \right|
```
   - Maximum error over the entire domain, 
```math
     \Phi^{\text{error}}_{\text{max}} = \max{\{\Phi^{\text{error}}_{i,j}\}}
```

5. **Plot Results**:
   - Plot velocity and vorticity fields (both theoretical and numerical) along with the velocity quiver.
   - Generate error contours plots as well.

---

## Results

![result](images/results.svg)

---

## BAD Coding Practices to Avoid

### 1. **No Comments or Documentation**
   - Unexplained logic, lack of variable descriptions, and unclear equations make code hard to follow.
   - Lack of input variable descriptions and their units (in case of dimensional physical quantities).
   - Hard for readers to understand the meaning of the variables, functions, complex equations etc.

### 2. **Poor Variable Naming**
   - Avoid cryptic names like `p`, `q`, `v_th` etc, for physical parameters.
   - Use clear, consistent, and descriptive names for parameters and functions.
   - Abbreviated and unclear names make it hard to understand the logic, and complicate debugging.

### 3. **Hardcoded Values**
   - Avoid hard coded value embedding for constants (e.g., grid size, viscosity) directly in the code.
   - Parameterize for flexibility. Any change (e.g., resolution) should not require editing multiple parts of the code
   - Avoid hard coded values of universal constants like $\pi$ and $e$, take valuse from available math libraries.

### 4. **Poor Formatting**
   - Ensure proper indentation, consistent alignment, and readable spacing.
   - Avoid mixed case and unformatted number in print statements. Use proper formatting.

### 5. **Unformatted Plots**
   - Include axis labels, descriptive titles, and uniform scaling in plots.

### 6. **No Functions or Modularity**
   - Modularize repetitive logic into reusable functions.
   - Validate inputs for physical constraints.
   - To make the functions generalized, add checks to validate for incorrect inputs. Add validation for physical constraints (e.g., viscosity must be positive)

### 7. **Inefficiency**
   - Explicit loops are error-prone and computationally expensive for large grids. Avoid explicit loops for operations.
   - Code is harder to debug and maintain.
   - Use vectorized operations for performance.

### 8. **No Error Handling**
   - The code assumes everything will work without exceptions. On should put checks and handle exceptions gracefully for edge cases or invalid inputs.

---

Adopting the best practices (avoiding the common pitfalls listed above) ensures code readability, flexibility, and efficiency. Use this as a guideline for developing robust and maintainable Python scripts.
