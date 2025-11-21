This is the code corresponding to the paper: 

https://arxiv.org/abs/2505.02587

The folder stores chiefly the following information: 
1_Code: Code from data download to model fitting - This is your main file for understanding the fitting procedure
2_Data: Training data stored
   -> Data is stored in: https://www.dropbox.com/scl/fo/a5do9x3rik8wpzzb9685n/h?rlkey=fle6gv7gypr41etuivjqxg34a&dl=0 
3_Results: Results from the COVID-19 Data application and Simulation
5_Supplementary Material: Supplementary Material

-> Caveat: 
Both the simulation and the COVID-19 application in 1. Code has its own homemade functions. Make sure that your machine is wiped and you have the correct functions loaded into your R-session. 
They actually perform the same tasks but have some difference to account for the difference in the tasks and data structured performed in the scripts.  


Fantastic figures and where to find them (in the Scripts- 1. Code): 
Figure_1 -> 2_1_DataVis.R
Figure_2 -> 4_1_Simulation_paper_Unadj.R
Figure_3 -> 4_Simulation_Paper.R
Figure_4 -> 4_Simulation_Paper.R
Figure_5 -> 4_3_Simulation_NegativeBinomial.R
Figure_6 -> 3_Model_Covid.R
Figure_7 -> 3_Model_Covid.R 
Figure_8 -> 3_Model_Covid.R 
Figure_9 -> 4_Simulation_Paper.R
Figure_10 -> 4_Simulation_Paper.R
Figure_11- Figure_14 -> 4_3_Simulation_NegativeBinomial.R


Tables: 
Table 1 -> 4_Simulation_Paper.R
Table 2 -> 4_3_Simulation_NegativeBinomial.R
Table 3 -> 3_Model_Covid.R 

Supplementary Material:
5_Supplementary_Material.R


Let's have a look at the content on 1_Code:


Set up: 
0_Libraries.R  ---> Libraries for modelling 

0_3_Home_Functions_Covid.R # ---> Homemade functions for modelling COVID-Data 3_Model_Covid.R also 5_Supplementary_Material.R
0_4_Home_Functions_Sim.R # ---> Homemade functions for modelling Simulated data 4_0_Simulation_Paper.R
0_41_Home_Functions_Sim_Unadj.R # ---> Homemade functions for modelling simulated data (in an unadjusted way and then adjust) 4_1_Simulation_Paper_Unadj.R
0_43_Home_Functions_Sim_NegativeBinomial.R # ---> Homemade functions for modelling misspecified simulated data 4_3_Simulation_NegativeBinomial.R
0_3_Home_Functions_Sim_NegativeBinomial.R # ---> 5_Supplementary_Material.R

Data download and wrangling: (You can jump these two data sets if you plan to use train_data.R directly) [Running time very fast]
1_Data_Download.R
2_Data_Wrangling.R

Data visualisation: [Running time very fast]
2_1_DataVis.R

Modelling the COVIDData with lage 30:  
3_Model_Covid.R [Running time around 3h on MAC OS M2 and 20h on Thinkpad T490s]

Modelling the COVIDData with lage 40:
5_Supplementary_Material.R [Running time around 3.5h on MAC OS M2 and 20h on Thinkpad T490s]


Modelling simulated data:
4_0_Simulation_Paper.R (Simulation) [Running time around 1h on MAC OS M2 and 7h on Thinkpad T490s]
4_1_Simulation_Paper_Unadj.R (Bias adjustment check) [Running time around 1.5h on MAC OS M2 and 7h on Thinkpad T490s]
4_3_Simulation_NegativeBinomial.R  (Misspecified)  [Running time around 5h on MAC OS M2 and 30h on Thinkpad T490s]

Flow chart:


                --> 1_Data_Download.R --> 2_Data_Wrangling.R  --> 2_1_DataVis.R
              /                                           \_
0_Libraries.R   ------- train_data.R -------------------->  \-->  3_Model_Covid.R <------------- 0_3_Home_Functions_Covid.R 
              \                                              \--> 5_Supplementary_Material.R <-/
               \
                \  0_4_Home_Functions_Sim.R ---------------------> 4_0_Simulation_Paper.R
                 \ 0_41_Home_Functions_Sim_Unadj.R --------------> 4_1_Simulation_Paper_Unadj.R
                   0_43_Home_Functions_Sim_NegativeBinomial.R ---> 4_3_Simulation_NegativeBinomial.R




(Lazyweave and RSQL do not exist in 4.5.2 anymore- Lazyweave is not important and can be replaced by kable BUT MYSQL is important and can be replaced by RMariaDB) 


The author's current R-version:

R version 4.5.1 (2025-06-13 ucrt) -- "Great Square Root"
Copyright (C) 2025 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

[Workspace loaded from ~/.RData]

