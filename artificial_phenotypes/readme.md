This folder contains the notebooks to an artificial phenotype (= target) from real human gut metagenomes. 

**./data/**
The two main tables used for analysis:
- tax_meta.qs: the taxonomic profiles from metagenomes
- taxa_table.qs: the taxonomy from the GTDB 

**./main_replicate/**
Notebooks for the main replicate presented in the paper. The notebooks for each model trained are named with [feature selection method]_[classifier].ipynb. Models further interpreted with endoR or SHAP also have an indication in their names (quite self-explaining). 

endoR intermediary object are included to save time and not have to run everything:
- Pre_bootstraps.qs: intermediary object from endoR from the gRRF RF model = discretized data and decisions + data row numbers for each bootstrap
- All_bootstraps.qs: intermediary object from endoR from the gRRF RF model = endoR bootstraps

**./repetitions/** 
Notebooks for repeated replicates. 
- `./repetitions/comparison_methods/`: notebooks to compare endoR to state-of-the-art methods. Since certain methods were ran separately (gLASSO, SHAP for random forests, and sparCC), their output need to be loaded from the `./repetitions/tmp` folder
- `./evaluation_parameters/`: notebooks to evaluate the effect of data and endoR parameters. A template to generate replicates is provided (generated data available on demand). The `Replicate[number]_B10-B100.ipynb` notebooks correspond to the examples of repeating endoR on a same replicate using 10 or 100 bootstraps. The procedure was performed on 6 replicates. Otherwise, notebooks' names are self-explaining.