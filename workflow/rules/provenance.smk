import json


# The resolved config, --config overrides included, so a run records what it ran with.
rule config_snapshot:
    output:
        f"{OUTDIR}/config.snapshot.json",
    params:
        config_json=lambda w: json.dumps(config, indent=2, default=str),
    log:
        f"{LOGDIR}/config_snapshot.txt",
    benchmark:
        f"{BENCHDIR}/config_snapshot.tsv",
    software:
        SOFTWARE_ENV
    shell:
        r"""
        exec &> >(tee {log:q})

        printf '%s\n' {params.config_json:q} > {output:q}

        sync {output:q}
        """
