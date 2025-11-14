import numpy as np
from numba import jit

@jit(nopython=True, fastmath=True, cache=True)
def relax2d(ax, ay, cx, cy, b, d_rhs, n, m, U_initial, omega=1.8, eps=1e-13, max_iter=10_000):
    """
    Successive Over-Relaxation (SOR) solver for a 2D discrete linear system.

    The function solves a PDE-like linear equation on a structured 2D grid:
        a_x * U[i-1,j] + c_x * U[i+1,j] +
        a_y * U[i,j-1] + c_y * U[i,j+1] -
        b     * U[i,j] = d_rhs

    Parameters
    ----------
    ax, ay, cx, cy : 2D arrays
        Discretization coefficients for neighboring nodes.
    b : 2D array
        Diagonal coefficients of the linear system.
    d_rhs : 2D array
        Right-hand side of the equation.
    n, m : int
        Grid dimensions.
    U_initial : 2D array
        Initial guess for the solution.
    omega : float
        Over-relaxation factor (1 < omega < 2 for standard SOR).
    eps : float
        Convergence tolerance based on max-norm of updates.
    max_iter : int
        Maximum allowed iterations.
    min_iter : int
        Minimum number of iterations (useful for smoothing in multigrid).

    Returns
    -------
    U : 2D array
        Relaxed solution field.
    """
    U = U_initial.copy()  
    err = 1.0
    it = 0

    inv_b = np.zeros_like(b)
    inv_b[1:-1, 1:-1] = 1.0 / b[1:-1, 1:-1]

    U_old = U.copy()
    
    while err > eps and it < max_iter or it < 100:
        it += 1
        
        for i in range(n):
            for j in range(m):
                U_old[i, j] = U[i, j]

        # Interior nodes (i=1..n-2, j=1..m-2)
        for i in range(1, n - 1):
            for j in range(1, m - 1):
                t = (ax[i, j] * U[i - 1, j] + cx[i, j] * U[i + 1, j] +
                     ay[i, j] * U[i, j - 1] + cy[i, j] * U[i, j + 1] -
                     d_rhs[i, j]) * inv_b[i, j]
                

                U[i, j] = t * omega + U[i, j] * (1 - omega)

        err = 0.0
        for i in range(1, n-1):
            for j in range(1, m-1):
                diff = abs(U[i, j] - U_old[i, j])
                if diff > err:
                    err = diff
        if it % 500 == 0:
            print("iter=", it, "err=", err)

    if it >= max_iter:
        print("Warning: Relaxation method did not converge in", max_iter, "iterations. Final error:", err)
    else:
        print("Converged in", it, "iterations. Final error:", err)


    return U


@jit(nopython=True, fastmath=True, cache=True)
def relax3d(ax, ay, az, cx, cy, cz, b, d_rhs, nx, ny, nz, U_initial, omega=1.8, eps=1e-13, max_iter=10000, min_iter=100):
    """
    Successive Over-Relaxation (SOR) solver for a 3D discrete linear system.

    The function solves a PDE-like linear equation on a structured 3D grid:
        a_x * U[i-1,j,k] + c_x * U[i+1,j,k] +
        a_y * U[i,j-1,k] + c_y * U[i,j+1,k] +
        a_z * U[i,j,k-1] + c_z * U[i,j,k+1] -
        b     * U[i,j,k] = d_rhs

    Parameters
    ----------
    ax, ay, az : 3D arrays
        Coefficients for negative x, y, z neighbors.
    cx, cy, cz : 3D arrays
        Coefficients for positive x, y, z neighbors.
    b : 3D array
        Diagonal coefficients of the system.
    d_rhs : 3D array
        Right-hand side of the equation.
    nx, ny, nz : int
        Grid dimensions.
    U_initial : 3D array
        Initial guess for the solution.
    omega : float
        Over-relaxation factor.
    eps : float
        Convergence tolerance (max update).
    max_iter : int
        Hard iteration limit.
    min_iter : int
        Minimum number of smoothing iterations.

    Returns
    -------
    U : 3D array
        Relaxed solution field.
    """
    U = U_initial.copy()
    U_old = U.copy()

    inv_b = np.zeros_like(b)
    inv_b[1:-1, 1:-1, 1:-1] = 1.0 / b[1:-1, 1:-1, 1:-1]

    err = 1.0
    it = 0

    while ((err > eps and it < max_iter) or it < min_iter):
        it += 1

        for i in range(nx):
            for j in range(ny):
                for k in range(nz):
                    U_old[i, j, k] = U[i, j, k]

        for i in range(1, nx - 1):
            for j in range(1, ny - 1):
                for k in range(1, nz - 1):
                    t = ( ax[i,j,k]*U[i-1,j,k] + cx[i,j,k]*U[i+1,j,k]
                        + ay[i,j,k]*U[i,j-1,k] + cy[i,j,k]*U[i,j+1,k]
                        + az[i,j,k]*U[i,j,k-1] + cz[i,j,k]*U[i,j,k+1]
                        - d_rhs[i,j,k]) * inv_b[i,j,k]
                    U[i,j,k] = omega*t + (1.0 - omega)*U[i,j,k]

        err = 0.0
        for i in range(1, nx-1):
            for j in range(1, ny-1):
                for k in range(1, nz-1):
                    diff = abs(U[i,j,k] - U_old[i,j,k])
                    if diff > err:
                        err = diff

        if it % 500 == 0:
            print("iter=", it, "err=", err)

    if it >= max_iter and err > eps:
        print("Warning: SOR did not converge in", max_iter, "iterations. Final error:", err)
    else:
        print("Converged in", it, "iterations. Final error:", err)

    return U


@jit(nopython=True, fastmath=True, cache=True)
def thomas(a, b, c, d):
    """
    Thomas algorithm for solving a tridiagonal linear system A x = d.

    Parameters
    ----------
    a : 1D array
        Sub-diagonal coefficients (a[0] is unused or should be 0).
    b : 1D array
        Main diagonal coefficients.
    c : 1D array
        Super-diagonal coefficients (c[-1] is unused or should be 0).
    d : 1D array
        Right-hand side vector.

    Returns
    -------
    y : 1D array
        Solution vector.
    """
    n = len(b)

    # Forward sweep
    for i in range(1, n):
        w = a[i] / b[i - 1]
        b[i] = b[i] - w * c[i - 1]
        d[i] = d[i] - w * d[i - 1]

    # Back substitution
    y = np.zeros(n)
    y[-1] = d[-1] / b[-1]
    for i in range(n - 2, -1, -1):
        y[i] = (d[i] - c[i] * y[i + 1]) / b[i]

    return y