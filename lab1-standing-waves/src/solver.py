import numpy as np
from numba import jit

@jit(nopython=True, fastmath=True, cache=True)
def relax(ax, ay, cx, cy, b, d_rhs, n, m, U_initial, alfa=1.8, eps=1e-13, max_iter=10_000):
    """
    Solve the linear system using the Successive Over-Relaxation (SOR) method.
    """
    U = U_initial.copy()  # Work with a copy to keep the initial guess intact
    err = 1.0
    it = 0
    inv_b = 1.0 / b  # 25% faster

    while err > eps and it < max_iter:
        it += 1
        err = 0.0
        # Interior nodes (i=1..n-2, j=1..m-2)
        for i in range(1, n - 1):
            for j in range(1, m - 1):
                t = (ax[i, j] * U[i - 1, j] + cx[i, j] * U[i + 1, j] +
                     ay[i, j] * U[i, j - 1] + cy[i, j] * U[i, j + 1] -
                     d_rhs[i, j]) * inv_b[i, j]
                
                deviation = abs(U[i, j] - t)
                err = max(err, deviation)
                U[i, j] = t * alfa + U[i, j] * (1 - alfa)
    
    if it >= max_iter:
        print(f"Warning: Relaxation method did not converge in {max_iter} iterations. Final error: {err}")

    return U
