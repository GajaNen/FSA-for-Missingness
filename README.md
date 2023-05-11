# Thesis

This repository contains the code for my thesis titled *Identifying Predictors of Missingness with Feature Selection Algorithms: A Simulation Study*. 

*Design:*

In each study, each scenario was simulated with 500 repetitions. In each study, there was one incomplete variable (*Y*) with its non-response indicator (*R*).

In Simulation 1, the effects of the predictor correlation, percentage of relevant variables, percentage of missing values and missingness mechanism on feature selection in artificial data were explored. Variables are denoted as *relevant*, if they have some dependency with *Y*, and under MAR and MNAR, with *R*. The factors were fully crossed, except for the percentage of relevant variables, which was set to 5% in MCAR scenarios. 

-   All information about varying conditions is in the [setup file](https://github.com/GajaNen/Thesis/blob/main/Simulation%201/setup.R).
-   [Dependencies for the whole simulation](https://github.com/GajaNen/Thesis/blob/main/Simulation%201/depend.R).
-   [Code with functions of data generating and fitting algorithms](https://github.com/GajaNen/Thesis/tree/main/Simulation%201/getOutput).
-   [Code to run the simulation from start to end](https://github.com/GajaNen/Thesis/blob/main/Simulation%201/runSim.R).

In Simulation 2, the effects of the percentage of missing values, explained variance in the latent variable for missingness and missingness mechanism on feature selection in real predictors & simulated outcome were explored. In each repetition, a resample with simple random sampling with replacement (*N*=299) of [the heart failure data](https://archive.ics.uci.edu/ml/datasets/Heart+failure+clinical+records) was taken for features and missingness was simulated the same way as in Simulation 1. Factors were fully crossed, except for explained variance, which was set to 0 in MCAR scenarios.

-   All information about varying conditions is in the [setup file](https://github.com/GajaNen/Thesis/blob/main/Simulation%202/setup.R).
-   [Heart failure data](https://github.com/GajaNen/Thesis/blob/main/Simulation%202/heart_failure.csv).
-   [Dependencies for the whole simulation](https://github.com/GajaNen/Thesis/blob/main/Simulation%202/depend.R).
-   [Code with functions of data generating and fitting algorithms](https://github.com/GajaNen/Thesis/tree/main/Simulation%202/getOutput).
-   [Code to run the simulation from start to end](https://github.com/GajaNen/Thesis/blob/main/Simulation%202/runSim.R).

