This folder contains the notebooks to predict the presence/absent of *Methanobacteriaceae* in human guts using taxonomic and metabolic profiles obtained from metagenomes and host variables.  

**./data/**
Data used for analysis (raw data available on demand). 

**./models** 
Notebooks with results from CV for models:
- without feature selection and random forest: none_ranger.ipynb
- feature selection performed using the taxa-aware gRRF algorithm, random forest fitted with ntrees = 250 or 500 : ta-gRRF_nt250.ipynb and ta-gRRF_nt200.ipynb
- feature selection performed using the taxa-aware gRRF algorithm and gradient boosted model: ta-gRRF_xgboost.ipynb 

Best predictive performances were reached with the ta-gRRF FS and random forest. The final model fitted on all taxa and its analysis with endoR are contained in the Final_model.ipynb notebook. 

Order_count_vs_rdm.ipynb: notebook to assess whether orders are significantly more represented in the set of endoR selected features compared to what would be expected by random.

**./marker_genes/**
Raw data blast are available on demand. The *_n_copies.ipynb notebooks contain the formatting of data: representative genomes were grouped into species and the number of copies of each gene were averaged across genomes. Sub-sequently, for genera and families, the weighted average of the number of copies were calculated across species, with weights corresponding to the relative abundances of species across samples. 

**./tmp/**
Intermediary objects such as CV models, the final model, endoR bootstraps, etc. These can be directly loaded instead of running the whole notebooks. 
The following files could not be included due to size but are available on demand (and can be generated using the notebooks in any case):
- ta_ranger_nt500.qs
- Final_model-permutated.qs
- endoR_boots.qs
- ta_ranger_nt250.qs
- Sub_ta-gRRF_range.qs

