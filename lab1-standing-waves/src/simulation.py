import numpy as np
from tqdm import tqdm
from .solver import relax

# --- Helper functions for boundary conditions ---

def apply_neumann_coeff_mods(b, ax, ay, cx, cy, nx, ny):
    """Coefficient modifications for Neumann boundaries (free edges)."""
    b[:, 1] -= ay[:, 1]
    ay[:, 1] = 0.0
    b[:, ny - 2] -= cy[:, ny - 2]
    cy[:, ny - 2] = 0.0
    b[1, :] -= ax[1, :]
    ax[1, :] = 0.0
    b[nx - 2, :] -= cx[nx - 2, :]
    cx[nx - 2, :] = 0.0

def copy_neumann_boundaries(U, nx, ny):
    """Neumann copying (mirror extrapolation + corners)."""
    U[:, 0] = U[:, 1]
    U[:, ny - 1] = U[:, ny - 2]
    U[0, :] = U[1, :]
    U[nx - 1, :] = U[nx - 2, :]
    U[0, 0] = 0.5 * (U[1, 0] + U[0, 1])
    U[nx - 1, ny - 1] = 0.5 * (U[nx - 2, ny - 1] + U[nx - 1, ny - 2])
    U[0, ny - 1] = 0.5 * (U[0, ny - 2] + U[1, ny - 1])
    U[nx - 1, 0] = 0.5 * (U[nx - 1, 1] + U[nx - 2, 0])

def set_sources(k, ax, ay, cx, cy, b, d_rhs, config):
    """Set sources according to the experiment_type."""
    nx, ny = config['nx'], config['ny']
    AA, KK, dt = config['amplitude'], config['frequency_k'], config['dt']
    experiment_type = config['experiment_type']

    source_value = AA * np.sin(np.pi * KK * k * dt)

    if experiment_type in ('FREE_CENTER', 'CLAMPED_CENTER'):
        I, J = nx // 2, ny // 2
        ax[I, J], ay[I, J], cx[I, J], cy[I, J] = 0.0, 0.0, 0.0, 0.0
        b[I, J] = 1.0
        d_rhs[I, J] = source_value
    elif experiment_type == 'FREE_TWO_GENERATORS':
        I_left, I_right, J_mid = 1, nx - 2, ny // 2
        # Left source
        ax[I_left, J_mid], ay[I_left, J_mid], cx[I_left, J_mid], cy[I_left, J_mid] = 0.0, 0.0, 0.0, 0.0
        b[I_left, J_mid] = 1.0
        d_rhs[I_left, J_mid] = source_value
        # Right source
        ax[I_right, J_mid], ay[I_right, J_mid], cx[I_right, J_mid], cy[I_right, J_mid] = 0.0, 0.0, 0.0, 0.0
        b[I_right, J_mid] = 1.0
        d_rhs[I_right, J_mid] = source_value

# --- Main simulation function ---

def run_simulation(config: dict):
    """
    Run the full membrane simulation loop.
    """
    nx, ny = config['nx'], config['ny']
    Lx, Ly = config['Lx'], config['Ly']
    dt, max_steps = config['dt'], config['max_steps']

    hx = Lx / (nx - 1)
    hy = Ly / (ny - 1)

    x = np.linspace(0.0, Lx, nx)
    y = np.linspace(0.0, Ly, ny)

    # Initialize time layers
    U = np.zeros((nx, ny))
    U0 = U.copy()
    U1 = U.copy()
    dU = U.copy()
    U1[1:nx - 1, 1:ny - 1] = U0[1:nx - 1, 1:ny - 1] + dt * dU[1:nx - 1, 1:ny - 1]

    # Base coefficients
    ax0 = 0.5 / (hx ** 2) * np.ones((nx, ny))
    cx0 = ax0.copy()
    ay0 = 0.5 / (hy ** 2) * np.ones((nx, ny))
    cy0 = ay0.copy()
    b0 = ax0 + cx0 + ay0 + cy0 + 1.0 / (dt ** 2)

    # Storage for results
    U_stack = np.zeros((nx, ny, max_steps), dtype=float)

    print(f"Starting simulation for '{config['experiment_type']}'...")
    for k in tqdm(range(1, max_steps + 1)):
        ax, ay, cx, cy, b = ax0.copy(), ay0.copy(), cx0.copy(), cy0.copy(), b0.copy()

        d_rhs = np.zeros((nx, ny), dtype=float)
        d_rhs[1:nx - 1, 1:ny - 1] = (
            (-2.0 * U1[1:nx - 1, 1:ny - 1] + U0[1:nx - 1, 1:ny - 1]) / (dt ** 2)
            - 0.5 * (
                (U0[0:nx - 2, 1:ny - 1] - 2.0 * U0[1:nx - 1, 1:ny - 1] + U0[2:nx, 1:ny - 1]) / (hx ** 2)
                + (U0[1:nx - 1, 0:ny - 2] - 2.0 * U0[1:nx - 1, 1:ny - 1] + U0[1:nx - 1, 2:ny]) / (hy ** 2)
            )
        )

        set_sources(k, ax, ay, cx, cy, b, d_rhs, config)

        # Boundary conditions
        if config['boundary_condition'] == 'neumann':
            apply_neumann_coeff_mods(b, ax, ay, cx, cy, nx, ny)
        
        # Solve the linear system
        U = relax(ax, ay, cx, cy, b, d_rhs, nx, ny, U1, 
                  alfa=config['relax_alfa'], 
                  eps=config['relax_eps'], 
                  max_iter=config['relax_max_iter'])

        # Boundary correction
        if config['boundary_condition'] == 'neumann':
            copy_neumann_boundaries(U, nx, ny)
        elif config['boundary_condition'] == 'dirichlet':
            U[0, :], U[-1, :], U[:, 0], U[:, -1] = 0.0, 0.0, 0.0, 0.0
        
        U_stack[:, :, k - 1] = U
        U0, U1 = U1.copy(), U.copy()

    return {"U_stack": U_stack, "x": x, "y": y, "config": config}
