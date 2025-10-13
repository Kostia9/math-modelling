import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def create_2d_animation(simulation_results: dict, output_filename: str, light_k: int = 2):
    """
    Create a 2D animation of the membrane oscillations and save it to file.
    """
    U_stack = simulation_results['U_stack']
    config = simulation_results['config']
    Lx, Ly = config['Lx'], config['Ly']
    
    # Frame thinning for smoother playback and smaller file
    U_stack_light = U_stack[..., ::light_k]

    frames = np.transpose(U_stack_light, (2, 1, 0))
    vmin, vmax = frames.min(), frames.max()

    fig, ax = plt.subplots()
    im = ax.imshow(frames[0], origin='lower', cmap='gray', vmin=vmin, vmax=vmax,
                   interpolation='nearest', animated=True, extent=[0, Lx, 0, Ly])
    txt = ax.text(0.02, 0.95, '', transform=ax.transAxes, color='w')
    fig.colorbar(im, ax=ax)
    ax.set_title(f"2D Animation: {config['experiment_type']}")
    ax.set_xlabel("x"); ax.set_ylabel("y")

    def init():
        im.set_data(frames[0])
        txt.set_text('t-step: 0')
        return im, txt

    def update(i):
        im.set_data(frames[i])
        im.set_clim(vmin, vmax)
        txt.set_text(f"t-step: {i * light_k}")
        return im, txt

    interval = 1000 * light_k / 30  # Adjust interval based on light_k for ~30 fps
    ani = FuncAnimation(fig, update, init_func=init,
                        frames=frames.shape[0], interval=interval, blit=True)

    print(f"Saving 2D animation to {output_filename}...")
    ani.save(output_filename, writer='ffmpeg', fps=30)
    plt.close(fig)
    print("2D animation saved.")

def create_3d_animation(simulation_results: dict, output_filename: str, light_k: int = 2, ds: int = 4):
    """
    Create a 3D animation of the membrane oscillations and save it to file.
    """
    U_stack = simulation_results['U_stack']
    x, y = simulation_results['x'], simulation_results['y']
    config = simulation_results['config']
        
    U_stack_light = U_stack[::ds, ::ds, ::light_k]
    X, Y = np.meshgrid(y[::ds], x[::ds])

    vmin, vmax = U_stack_light.min(), U_stack_light.max()

    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection='3d')
    
    surf = [ax.plot_surface(X, Y, U_stack_light[:, :, 0], cmap='viridis',
                            linewidth=0, antialiased=False, vmin=vmin, vmax=vmax)]
    ax.set_zlim(vmin, vmax)
    ax.set_xlabel('x'); ax.set_ylabel('y'); ax.set_zlabel('U')
    fig.colorbar(surf[0], shrink=0.6, pad=0.1)

    def update_3d(i):
        surf[0].remove()
        Z = U_stack_light[:, :, i]
        surf[0] = ax.plot_surface(X, Y, Z, cmap='viridis', linewidth=0,
                                  antialiased=False, vmin=vmin, vmax=vmax)
        ax.set_title(f"3D Animation: {config['experiment_type']} (Step {i*light_k})")
        return surf

    interval = 1000 * light_k / 30  # Adjust interval based on light_k for ~30 fps

    ani = FuncAnimation(fig, update_3d, frames=U_stack_light.shape[2],
                        interval=interval, blit=False)
    
    print(f"Saving 3D animation to {output_filename}...")
    ani.save(output_filename, writer='ffmpeg', fps=30, dpi=100)
    plt.close(fig)
    print("3D animation saved.")
