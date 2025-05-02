## Git repository for companion code and data for 'Modelling uncertainty around free-list cultural salience scores' paper

This repository contains four data files and one script file.

Script 'FreeList_Uncertainty_CompanionCode.R' provides R code to replicate the analyses in the main paper (in addition to some extensions), as well as demonstrating how to apply some of the new AnthroTools functions to automate some of these processes.

The four data files are:

 - 'FL10_Virtues.txt' - Tyvan free-list responses for what is means to be a "good Tyvan person"
 - 'demo.txt' - Tyvan demographic data (here, we're only interested in the variable 'sex')
 - 'tyva_domains.txt' - Tyvan free-list data regarding what Buddha (big/moralising god), spirits (local god) and the police dislike. This data is used to demonstrate how to incorporate uncertainty in Jaccard's similarity
 - 'Cross-cultural_ERM1.csv' - Cross-cultural Evolution of Religion and Morality (ERM) data for eight societies. Used to demonstrate how to incorporate uncertainty into Cultural FST estimates, using whether participants listed 'Morality' as something big/moralising gods dislike.

Note that all of these datasets were previously openly-available, or adapted from openly-available datasets (see https://github.com/bgpurzycki/Tyvan-Values and https://github.com/bgpurzycki/free-list_QASS).