def get_config(experiment_type: str) -> dict:
    """
    Return a configuration dictionary for the given experiment type.
    """
    # --- Base parameters ---
    config = {
        "nx": 101,
        "ny": 101,
        "Lx": 1.0,
        "Ly": 1.0,
        "amplitude": 2e-4,   # Source amplitude (A)
        "frequency_k": 4,    # Source temporal frequency multiplier (KK)
        "dt": 1e-2,
        "max_steps": 500,
        "relax_alfa": 1.8,
        "relax_eps": 1e-13,
        "relax_max_iter": 5000,
        "experiment_type": experiment_type,
        "light_k": 2,
        "ds": 4
    }

    # --- Experiment-specific parameters ---
    if experiment_type in ('FREE_CENTER', 'FREE_TWO_GENERATORS'):
        config['boundary_condition'] = 'neumann'
    elif experiment_type == 'CLAMPED_CENTER':
        config['boundary_condition'] = 'dirichlet'
    else:
        raise ValueError(f"Unknown experiment_type: {experiment_type}")

    return config
