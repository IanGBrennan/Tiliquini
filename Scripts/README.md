# Scripts

Scripts used to analyze morphological data and generate figures for this manuscript. 

All analyses can be recreated using just the scripts provided here. Most scripts are sourced within the others that are run chronologically (00_...---11_...). 

+ *00_Data_Preparation.R*: Process raw data and generate log shape ratios.

+ *01_Process_BayesTraits_Output.R*: Read in and process the output from BayesTraits runs of the VariableRates model.

+ *02_Model_Fitting.R*: Fit a series of comparative models to the univariate trait data and plot support.

+ *02_VarRates_Trees.R*: Plot the output of the Variable Rates model as in the left column of Figure 4. 

+ *03_Extract_Ancestral_Traits.R*: Estimate ancestral trait values at nodes and along branches following the rate-transformed trees from BayesTraits.

+ *04_Disparity_Through_Time_Traits.R*: Plot the trends in accumulation of disparity through for traits as in Figure 2.

+ *05_Rate_Through_Time_Traits.R*: Plot the trends in evolutionary rates through time for each trait. 

+ *06a_Disparity_Rate_Through_Time_Modules.R*: Plot the disparity through time and rates through time as in Figure 2 and Fig.S31 for modules for all Tiliquini.

+ *06b_Disparity_Rate_Through_Time_Modules_Oz.R*: Plot the disparity through time and rates through time as in Figure 2 and Fig.S31 for modules for just Australian Tiliquini.

+ *07_Rate_Trajectory_Traitgram.R*: Plot the trait trajectories through time for a given trait, and colored by rate, as in the right column of Figure 4. 

+ *08_Rate_to_Node.R*: Plot the rate trajectories through time for shifted nodes as in the central column of Figure 4. 

+ *09_Sim_to_Node.R*: Simulate traits and compare them against the evolutionary path of an observed trait as in the left column of Figure 6. 

+ *10_Distance_Btwn_Nodes.R*: Plot the difference in observed - expected trait change between nodes as in the central column of Figure 6. 

+ *11_Innovate_Elaborate.R*: Analyse the trait data, identify primary PC axes, and plot the trends along individual branches as elaborative or innovative. 

---

Individual scripts are often called within the chronological series, but are listed and explained below.

+ *Calculate_AICs.R*: helper functions for calculating AIC values from log-likelihoods.

+ *morphotrajectory.R*: functions for the sim-to-node and distance-btwn-node plots presented in Figure 6, and the morphotrajectory plots in Figure 5.

+ *PLOT_innovate_elaborate.R*: functions for the elaborate/innovate plot of Figure 2. 

+ *PLOT_rate.trajectory.R*: script to plot the rate-trajectory plots (right column rate/trait-grams) of Figure 4. 

+ *plotting_BayesTraits.R*: functions for processing BayesTraits output files and visualizing rate heterogeneity along on the focal tree (left column trees) of Figure 4. 

+ *process_BTVarRates.R*: function for processing BayesTraits log file to assess convergence of estimated parameters.

+ *rate.trajectory.R*: functions used in the *PLOT_rate.trajectory.R* script, as well as for generating the rate-to-node plots of Figure 4 (center column).

+ *RUN_BayesTraits_VarRates_Likelihoods.R*: calculate AIC values for the VR model and compare against other comparative models. 

+ *RUN_pulsR.R*: script for running the common and Levy models across all morphological traits.

+ *SaveTraitsToFiles_BayesTraits.R*: helper function for writing individual morphological traits to BayesTraits input formatted files.

+ *shifts.to.simmap.l1ou.R*: functions for visualizing l1ou results.

+ *trait.at.time.R*: many functions for extracting trait, rate, and disparity metrics through time for traits and modules

+ *visualizeBayesTraitsLog.R*: used in conjunction with `process_BTVarRates.R` to visually assess convergence of parameter estimates


___

