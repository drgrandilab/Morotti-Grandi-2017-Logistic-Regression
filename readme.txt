Logistic regression analysis of populations of electrophysiological models to assess proarrythmic risk.

                                      Matlab code

Contents:
__________________________________________________________________________________________________________

                                      MAIN FILES
__________________________________________________________________________________________________________
1) SA_generate_parameters.m
File for generation of random perturbations of model parameters.
__________________________________________________________________________________________________________
2) SA_analyze_EAD.m
File for cyclic execution of EAD protocol and analysis of model outputs. It uses the matrix created with
SA_generate_parameters.m.
__________________________________________________________________________________________________________
3) SA_logistic.m
File for logistic regression analysis and plot of the results. This file uses the matrices created with
SA_generate_parameters.m and SA_analyze_EAD.m.
__________________________________________________________________________________________________________

                                      SUPPORTING FILES
__________________________________________________________________________________________________________
4) morotti_et_al_ham_ina_ran_model_SA.m
ODE model of excitation-contraction coupling in the human atrial cell used in SA_analyze_EAD.m. This file
was modified from the one published in Morotti et al. J Mol Cell Cardiol. 2016 Jul;96:63-71. for allowing
modulation of model parameters).
__________________________________________________________________________________________________________
5) EAD_occurrence.m
Function called by SA_analyze_EAD.m for the assessment of EAD occurrence.
__________________________________________________________________________________________________________
6) rotateXLabels.m
Function called by SA_logistic.m for rotating the labels on the X axis in some figures.
__________________________________________________________________________________________________________

                                      MAT FILES
__________________________________________________________________________________________________________
7) yf_ham_ina_ran_ACh0p1_1Hz.mat
Initial conditions (obtained with the baseline model at 1-Hz pacing in presence of 0.1 uM acetylcholine)
used for simulations in SA_analyze_EAD.m.
__________________________________________________________________________________________________________
8) SA_par_matrix_1000_s0p1.mat
Matrix of the scaling factors used in the example shown in this paper (this file was created and saved
with the SA_generate_parameters.m, and then used in SA_analyze_EAD.m and SA_logistic.m).
__________________________________________________________________________________________________________
9) SA_outputs_matrix_1000_s0p1.mat
Matrix of the model outputs obtained in the example shown in the paper (this file was created and saved 
with SA_analyze_EAD.m, and then used in SA_logistic.m).
__________________________________________________________________________________________________________


Reference:
Morotti S & Grandi E.
Logistic regression analysis of populations of electrophysiological models to assess proarrythmic risk.
MethodsX. 2016 Dec 23; 4:25-34. doi: http://dx.doi.org/10.1016/j.mex.2016.12.002

Please cite the above paper when using this code.
