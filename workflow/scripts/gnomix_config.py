"""Write a gnomix config with the run's seed, threads and cleanup set.

gnomix takes every setting from one YAML file. The rule fixes the seed so training
repeats, caps the cores at the job's threads instead of the whole node, and deletes
the simulated training data once the model is trained.
"""

import argparse

import yaml


def build_parser():
    """Build the gnomix_config.py command-line parser.

    Returns
    -------
    argparse.ArgumentParser
        Parser for the gnomix_config.py flags.
    """
    p = argparse.ArgumentParser(description="Write a gnomix config for one job.")
    p.add_argument("--base", required=True, help="Config YAML to start from.")
    p.add_argument("--seed", required=True, type=int, help="Training seed.")
    p.add_argument("--threads", required=True, type=int, help="Cores gnomix may use.")
    p.add_argument("--output", required=True, help="Config YAML to write.")
    return p


def main(argv=None):
    """Read the base config, set the job's values and write it.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments. Defaults to `sys.argv[1:]`.
    """
    args = build_parser().parse_args(argv)
    with open(args.base) as f:
        config = yaml.safe_load(f)
    config["seed"] = args.seed
    config["model"]["n_cores"] = args.threads
    config["simulation"]["rm_data"] = True
    with open(args.output, "w") as f:
        yaml.safe_dump(config, f, sort_keys=False)


if __name__ == "__main__":
    main()
