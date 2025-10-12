import argparse
import os
from src.config import get_config
from src.simulation import run_simulation
from src.visualization import create_2d_animation, create_3d_animation

def main():
    parser = argparse.ArgumentParser(description="Run 2D membrane simulation.")
    parser.add_argument(
        '--experiment_type',
        type=str,
        required=True,
        choices=['FREE_CENTER', 'CLAMPED_CENTER', 'FREE_TWO_GENERATORS'],
        help="Type of experiment to run."
    )
    parser.add_argument(
        '--output_dir',
        type=str,
        default='results',
        help="Directory to save animation files."
    )
    args = parser.parse_args()

    # Create the output directory if it does not exist
    os.makedirs(args.output_dir, exist_ok=True)

    # 1) Get configuration
    config = get_config(args.experiment_type)

    # 2) Run simulation
    results = run_simulation(config)
    
    # 3) Create and save animations
    base_filename = f"membrane_{args.experiment_type}"
    
    output_2d = os.path.join(args.output_dir, f"{base_filename}_2d.mp4")
    create_2d_animation(results, output_2d, light_k=config['light_k'])

    output_3d = os.path.join(args.output_dir, f"{base_filename}_3d.mp4")
    create_3d_animation(results, output_3d, light_k=config['light_k'], ds=config['ds'])

    print("\nAll tasks completed successfully!")

if __name__ == "__main__":
    main()
