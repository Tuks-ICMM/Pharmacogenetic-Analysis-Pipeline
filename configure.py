from json import load
from os.path import abspath, dirname, join

# from subprocess import run

with open(join("config", "config.json"), "r") as file_content:
    config = load(file_content)

PBS_Headers = [
    "#!/usr/bin/env bash",
    "#PBS -q long",
    "#PBS -l walltime=900:00:00",
    "#PBS -l nodes=1:ppn=1",
    "#PBS -N Snakemake",
    " ",
]

PBS_Body = [
    "module load python-3.8.2",
    "cd {};".format(dirname(abspath(__file__))),
    "snakemake --profile config/PBS-Torque-Profile --use-envmodules",
]

if "environment" in config:
    if "email" in config["environment"]:
        if "condition" in config["environment"]["email"]:
            email_conditions = "#PBS -k "
            for condition in config["environment"]["email"]["conditions"]:
                email_conditions += condition
            PBS_Headers.append(email_conditions)
        if "address" in config["environment"]["email"]:
            PBS_Headers.append("#PBS -M " + config["environment"]["email"]["address"])


with open("run.sh", "w") as file:
    file.writelines("\n".join(PBS_Headers + PBS_Body))

# run(["qsub", ".run.sh"], shell=True)
