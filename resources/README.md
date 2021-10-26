# Resources

Welcome to the resources folder. This folder houses any accompanying resources required for the running of this pipeline.

By default, this includes a copy of some legacy scripts in our case (`liftOverPlink.py` and `rmBadLifts.py`) as well as a copy of the original data used to calibrate the `CONDEL` score plugin used in the e! Ensemble PEARL [Variant Effect Prediction plugin](https://m.ensembl.org/info/docs/tools/vep/script/vep_plugins.html) ([Repository available here](https://github.com/Ensembl/VEP_plugins/blob/release/104/Condel.pm)).
> We store a copy of this data as we have implemented a Python reverse-engineered version of their calculation which derives constants from these reference datasets.

We also store a compatable copy of EigenStrad for internal Admixture analysis use.