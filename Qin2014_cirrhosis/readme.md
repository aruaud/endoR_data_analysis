This folder contains the data generated by Qin et al., 2014 and notebooks to analyse them using random forests and endoR.

**./data/**
Original data downloaded from the MLrepo (see: https://knights-lab.github.io/MLRepo/docs/qin_healthy_cirrhosis.html) 
\+ formatted data generated with the Data_formatting.ipynb notebook

**./models** 
Notebooks with results from CV for models:
- without feature selection: noFS.ipynb
- feature selection performed using Boruta: Boruta.ipynb
- feature selection performed using the taxa-aware gRRF algorithm: ta-gRRF.ipynb

Best predictive performances were reached with the ta-gRRF FS. The final model fitted on all taxa and its analysis with endoR are contained in the same ta-gRRF.ipynb notebook. 

Qin_vs_endor.ipynb: notebook to make plots comparing results from the endoR analysis versus results from Qin et al., 2014.

**./tmp/**
Intermediary objects such as CV models, the final model, endoR bootstraps, etc. These can be directly loaded instead of running the whole notebooks. 