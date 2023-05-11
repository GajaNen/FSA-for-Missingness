# Thesis

This repository contains the code for my thesis titled *Identifying Predictors of Missingness with Feature Selection Algorithms: A Simulation Study*. 

*Abstract:*

In biomedical and social sciences, researchers are often faced with missing data. Most of the techniques for dealing with such data assume that missingness depends only on the observed parts of the data, which need to be accounted for in the analysis. Thus, feature selection for missingness has to be performed. Current approaches for predictor identification may run into problems in the presence of class imbalance, multicolinearity, and non-linearity, and other feature selection algorithms (FSAs) can be used instead. This study explores the effects of predictor correlation, strength of the relationship between the predictors and missingness, missingness mechanism, percentage of true predictors and percentage of missing values on performance of 9 filter, wrapper and embedded FSAs and logistic regression. Two simulation studies were conducted, one on an artificially created data and one resampling study of a real-world data with simulated missingness. Subset correctness and completeness, stability, power to select all true predictors and similarity between FSAs were evaluated. The results show that percentage of missing values and percentage of relevant variables affect FSAs the most and, generally, random forest-based methods achieve the highest performance in terms of all three measures under a variety of scenarios, especially when the percentage of relevant variables is low. Alternatively, fast correlation-based filter, logistic regression and regularised SVM are a good choice under certain conditions.

*Design:*

In each study, there were 500 repetitions per scenario. 

In Simulation 1, the effects of predictor correlation, percentage of relevant variables, percentage of missing values and missingness mechanism on feature selection in artificial data were explored. The factors were fully crossed, except for the percentage of relevant variables, which was set to 5% in MCAR scenarios. 

-   All information abot varying condtions is in the [setup file](https://github.com/GajaNen/Thesis/blob/main/Simulation%201/setup.R).
-   [Code with functions of data generating and fitting algorithms](https://github.com/GajaNen/Thesis/tree/main/Simulation%201/getOutput).
-   [Code to run sim from start to end](https://github.com/GajaNen/Thesis/blob/main/Simulation%201/runSim.R).

In Simulation 2, the effects of missingness percentage, explained variance in the latent variable for missingness and missingness mechanism on feature selection in real predictors & simulated outcome were explored. In each repetition, a resample with simple random sampling with replacement (N=299) of [the heart failure data](https://archive.ics.uci.edu/ml/datasets/Heart+failure+clinical+records) ws used for features and missingness was simulated the same way as in Simulation 1.

-   All information abot varying condtions is in the [setup file](https://github.com/GajaNen/Thesis/blob/main/Simulation%202/setup.R).
-   [Code with functions of data generating and fitting algorithms](https://github.com/GajaNen/Thesis/tree/main/Simulation%202/getOutput).
-   [Code to run sim from start to end](https://github.com/GajaNen/Thesis/blob/main/Simulation%202/runSim.R).

