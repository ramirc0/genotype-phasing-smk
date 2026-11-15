# Config, path constants, and shared helpers.


configfile: "config/config.yaml"


RUN_ID = config.get("run_id", "default")
OUTDIR = f"{config['outdir']}/{RUN_ID}"
REFDIR = f"{config['outdir']}/refs"
LOGDIR = f"logs/{RUN_ID}"
BENCHDIR = f"benchmarks/{RUN_ID}"

# The workspace path is relative to the rule file. Every rule file sits in rules/.
SOFTWARE_ENV = pixi(workspace="../envs", env=config["pixi_env"], locked=True)
