# Project 1: Simulating Standing Waves on a 2D Membrane (Chladni Figures)

## Project Overview

This project numerically simulates the oscillations of a 2D membrane to visualize standing waves, also known as **Chladni figures**. The goal is to solve the wave equation using the finite-difference method and to investigate how different boundary conditions and excitation frequencies lead to the formation of resonant patterns.

---

## Theoretical Background

### 1. The Governing Equation (2D Membrane)

The model is based on the two-dimensional linear wave equation:

$$
\frac{\partial^2 U}{\partial t^2} = c^2\left( \frac{\partial^2 U}{\partial x^2} + \frac{\partial^2 U}{\partial y^2} \right)
$$

where $U(x,y,t)$ is the displacement of the membrane and $c$ is the wave speed.

---

### 2. Discretization of Space and Time

We use a uniform grid:

$$
x_i = i\,\Delta x,\quad y_j = j\,\Delta y,\quad t^k = k\,\Delta t
$$

and we denote the displacement at a grid point as:

$$
U_{i,j}^k \approx U(x_i, y_j, t^k)
$$

---

### 3. Finite Difference Approximations

**Second derivative in time (central difference):**

$$
\frac{\partial^2 U}{\partial t^2}\Big|_{i,j}^k \approx \frac{U_{i,j}^{k+1} - 2U_{i,j}^k + U_{i,j}^{k-1}}{\Delta t^2}
$$

**Second derivative in x:**

$$
\frac{\partial^2 U}{\partial x^2}\Big|_{i,j}^k \approx \frac{U_{i-1,j}^k - 2U_{i,j}^k + U_{i+1,j}^k}{\Delta x^2}
$$

**Second derivative in y:**

$$
\frac{\partial^2 U}{\partial y^2}\Big|_{i,j}^k \approx \frac{U_{i,j-1}^k - 2U_{i,j}^k + U_{i,j+1}^k}{\Delta y^2}
$$

---

### 4. The Crank-Nicolson Scheme

To enhance stability, we take the average of the Laplacian at time steps $k+1$ and $k-1$:

$$
\frac{U_{i,j}^{k+1} - 2U_{i,j}^{k} + U_{i,j}^{k-1}}{\Delta t^2}= \frac{c^2}{2}\Bigg(\frac{U_{i-1,j}^{k+1} - 2U_{i,j}^{k+1} + U_{i+1,j}^{k+1}}{\Delta x^2}+ \frac{U_{i,j-1}^{k+1} - 2U_{i,j}^{k+1} + U_{i,j+1}^{k+1}}{\Delta y^2}+\frac{U_{i-1,j}^{k-1} - 2U_{i,j}^{k-1} + U_{i+1,j}^{k-1}}{\Delta x^2}+\frac{U_{i,j-1}^{k-1} - 2U_{i,j}^{k-1} + U_{i,j+1}^{k-1}}{\Delta y^2}\Bigg)
$$

---

### 5. Algebraic Equation for a Node (i,j)

After rearranging the terms for the $U^{k+1}$ time step, we get:

$$
b\,U_{i,j}^{k+1} - a_x\,U_{i-1,j}^{k+1} - c_x\,U_{i+1,j}^{k+1} - a_y\,U_{i,j-1}^{k+1} - c_y\,U_{i,j+1}^{k+1} = -\,d_{rhs,\,i,j}
$$

where the coefficients are:

$$
a_x = c_x = \frac{c^2}{2\,\Delta x^2}, \qquad a_y = c_y = \frac{c^2}{2\,\Delta y^2}, \qquad b = \frac{1}{\Delta t^2} + 2a_x + 2a_y
$$

and the right-hand side (RHS) is defined as:

$$
d_{rhs,\,i,j} = \frac{-2U_{i,j}^{k} + U_{i,j}^{k-1}}{\Delta t^2} - \frac{c^2}{2}\left( \frac{U_{i-1,j}^{k-1} - 2U_{i,j}^{k-1} + U_{i+1,j}^{k-1}}{\Delta x^2} + \frac{U_{i,j-1}^{k-1} - 2U_{i,j}^{k-1} + U_{i,j+1}^{k-1}}{\Delta y^2} \right)
$$

---

### 6. Solving the System with the Relaxation Method

We are left with a system of linear equations. The update formula for a node $(i,j)$ using the Successive Over-Relaxation (SOR) method is:

$$
U_{i,j}^{(new)} = (1-\omega)\,U_{i,j}^{(old)} + \omega \cdot \frac{ a_x(U_{i-1,j} + U_{i+1,j}) + a_y(U_{i,j-1} + U_{i,j+1}) - d_{rhs,\,i,j} }{ b }
$$

where $0 < \omega < 2$ is the relaxation parameter (`alfa` in the code).

---

## How to Run

### Requirements
- Python 3.10+
- Packages: `numpy`, `numba`, `matplotlib`, `tqdm`
- `ffmpeg` installed and available on your system's PATH (for saving animations).

### Installation and Execution

1.  **Clone & Navigate**
    Navigate to the project directory.
    ```bash
    cd lab1-standing-waves
    ```

2.  **Create and Activate a Virtual Environment** (Recommended)
    * **Windows (PowerShell):**
        ```powershell
        python -m venv .venv
        ./.venv/Scripts/Activate.ps1
        ```
    * **macOS / Linux:**
        ```bash
        python3 -m venv .venv
        source .venv/bin/activate
        ```

3.  **Install Dependencies**
    ```bash
    python -m pip install --upgrade pip
    pip install numpy numba matplotlib tqdm
    ```

4.  **Run a Simulation**
    Choose one of the experiment types to run:
    ```bash
    python lab1-standing-waves/main.py --experiment_type FREE_CENTER
    ```
    Other options: `CLAMPED_CENTER`, `FREE_TWO_GENERATORS`.

5.  **Find Outputs**
    The resulting animation files (.mp4) will be saved in the `lab1-standing-waves/results/` directory.

### Notes
- To reduce runtime for testing, you can decrease the `max_steps` variable in `lab1-standing-waves/src/config.py`.
- If you encounter errors related to `ffmpeg`, ensure it is correctly installed and its location is added to your system's PATH environment variable.
