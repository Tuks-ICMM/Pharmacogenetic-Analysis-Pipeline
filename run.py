from os.path import join
from os import open
from json import load

with open(join("..", "config", "config.json"), "r") as file_content:
    config = load(file_content)

if config.environment:
    if config.environment.email:
        email = "-M "
        email_conditions = "-k "
        for condition in config.environment.email.conditions:
            email_conditions = email_conditions + condition
        if config.environment.email.name:
            email = email + config.environment.email.name


file = [
    "#!/usr/bin/env bash\n",
    "#PBS -q long\n",
    "#PBS -l walltime=900:00:00\n",
    "#PBS -l nodes=1:ppn=1\n",
    "#PBS -k {}\n".format("".join(config.environment.email.conditions)),
    "#PBS -M {}".format(email),
    "#PBS -N Snakemake\n",
    "module load python-3.8.2",
    "cd {};".format(config["environment"]["working-directory"]),
    "snakemake --cluster-config config/cluster.json --profile config/PBS-Torque-Profile",
]

with open("./.run.sh", "w") as run:
    run.writelines(file)

open("qsub ./.run.sh")
