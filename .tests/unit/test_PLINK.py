import os
import sys

import subprocess as sp
from tempfile import TemporaryDirectory
import shutil
from pathlib import Path, PurePosixPath

sys.path.insert(0, os.path.dirname(__file__))

import common


def test_PLINK():

    with TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir) / "workdir"
        data_path = PurePosixPath(".tests/unit/PLINK/data")
        expected_path = PurePosixPath(".tests/unit/PLINK/expected")

        # Copy data to the temporary workdir.
        shutil.copytree(data_path, workdir)

        # dbg
        print("results/FINAL/SUPER/ALL_CYP2A6.EUR.acount results/FINAL/SUPER/ALL_CYP2A6.EUR.hardy results/FINAL/SUPER/ALL_CYP2A6.EUR.smiss results/FINAL/SUPER/ALL_CYP2A6.EUR.vmiss results/FINAL/SUPER/ALL_CYP2A6.AFR.acount results/FINAL/SUPER/ALL_CYP2A6.AFR.hardy results/FINAL/SUPER/ALL_CYP2A6.AFR.smiss results/FINAL/SUPER/ALL_CYP2A6.AFR.vmiss results/FINAL/SUPER/ALL_CYP2A6.SAS.acount results/FINAL/SUPER/ALL_CYP2A6.SAS.hardy results/FINAL/SUPER/ALL_CYP2A6.SAS.smiss results/FINAL/SUPER/ALL_CYP2A6.SAS.vmiss results/FINAL/SUB/ALL_CYP2A6.GBR.acount results/FINAL/SUB/ALL_CYP2A6.GBR.hardy results/FINAL/SUB/ALL_CYP2A6.GBR.smiss results/FINAL/SUB/ALL_CYP2A6.GBR.vmiss results/FINAL/SUB/ALL_CYP2A6.GWD.acount results/FINAL/SUB/ALL_CYP2A6.GWD.hardy results/FINAL/SUB/ALL_CYP2A6.GWD.smiss results/FINAL/SUB/ALL_CYP2A6.GWD.vmiss results/FINAL/SUB/ALL_CYP2A6.GIH.acount results/FINAL/SUB/ALL_CYP2A6.GIH.hardy results/FINAL/SUB/ALL_CYP2A6.GIH.smiss results/FINAL/SUB/ALL_CYP2A6.GIH.vmiss", file=sys.stderr)

        # Run the test job.
        sp.check_output([
            "python",
            "-m",
            "snakemake", 
            "results/FINAL/SUPER/ALL_CYP2A6.EUR.acount results/FINAL/SUPER/ALL_CYP2A6.EUR.hardy results/FINAL/SUPER/ALL_CYP2A6.EUR.smiss results/FINAL/SUPER/ALL_CYP2A6.EUR.vmiss results/FINAL/SUPER/ALL_CYP2A6.AFR.acount results/FINAL/SUPER/ALL_CYP2A6.AFR.hardy results/FINAL/SUPER/ALL_CYP2A6.AFR.smiss results/FINAL/SUPER/ALL_CYP2A6.AFR.vmiss results/FINAL/SUPER/ALL_CYP2A6.SAS.acount results/FINAL/SUPER/ALL_CYP2A6.SAS.hardy results/FINAL/SUPER/ALL_CYP2A6.SAS.smiss results/FINAL/SUPER/ALL_CYP2A6.SAS.vmiss results/FINAL/SUB/ALL_CYP2A6.GBR.acount results/FINAL/SUB/ALL_CYP2A6.GBR.hardy results/FINAL/SUB/ALL_CYP2A6.GBR.smiss results/FINAL/SUB/ALL_CYP2A6.GBR.vmiss results/FINAL/SUB/ALL_CYP2A6.GWD.acount results/FINAL/SUB/ALL_CYP2A6.GWD.hardy results/FINAL/SUB/ALL_CYP2A6.GWD.smiss results/FINAL/SUB/ALL_CYP2A6.GWD.vmiss results/FINAL/SUB/ALL_CYP2A6.GIH.acount results/FINAL/SUB/ALL_CYP2A6.GIH.hardy results/FINAL/SUB/ALL_CYP2A6.GIH.smiss results/FINAL/SUB/ALL_CYP2A6.GIH.vmiss",
            "-F", 
            "-j1",
            "--keep-target-files",
    
            "--directory",
            workdir,
        ])

        # Check the output byte by byte using cmp.
        # To modify this behavior, you can inherit from common.OutputChecker in here
        # and overwrite the method `compare_files(generated_file, expected_file), 
        # also see common.py.
        common.OutputChecker(data_path, expected_path, workdir).check()
