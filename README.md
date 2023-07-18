# RBM_Optimization
The optimization framework for the coupled mizuRoute RBM model for the WRR paper.

Files are organized by directory. The statistical model (aka Reression_Model) uses the statistial relationship between air temperature and water temperature. The rest use the River Basin Model (RBM) to calculate river temperature through heat flux equations. The Control folder contains control files to run both mizuRoute and RBM. Then the workflow goes Optimization, Validation, Error_Quantification, Production_Runs.
