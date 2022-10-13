# OConnell2015
Matlab code from O'Connell and Reeves, 2015.
Michael O'Connell 
08/02/2016

Selected Contents:
decon_log_nest_JPattern.m
deconvolutionInterphase2Dcopy.m
deconvolutionResults.mat
genef_log_nest_avg.m
kanodia_log_nest_JPattern.m
dir: Other Files
script_isres_decon.m
sub_decon.sh


decon_log_nest_JPattern.m
This file contains the NESTed-functions version of the DECONvolution model, in which the
input parameters are assumed to be on a LOG scale. It also uses a Jacobian Pattern matrix, 
which tells the ODE solver the structure of the Jacobian. The published results are found
by running opening deconvolutionResults.mat and using xb(row,:) as input to 
decon_log_nest_JPattern(). 

deconvolutionInterphase2Dcopy.m
The 2D version of the deconvolution model that shows a shuttling-like phenotype. 

deconvolutionResults.mat
Results from different versions of the deconvolution model. XB is the set of outputs used
for publication. Same as PARENTS{1} and PARENTS{3}. PARENTS{5} is a set of parameters for
the Kanodia model (kanodia_log_nest_JPattern.m).


genef_log_nest_avg.m
Nested version of the gene expression function downstream of the
deconvolution model, using the average of the last 4 dl gradient time
points. 

kanodia_log_nest_JPattern.m
This file contains the NESTed-functions version of the KANODIA model, in which the
input parameters are assumed to be on a LOG scale. It also uses a Jacobian Pattern matrix, 
which tells the ODE solver the structure of the Jacobian. Parameters can be found in 
deconvolutionResults.mat > Parents{5}

Other Files: contains some files necessary to run the other .m files, such as gregsdata.mat

script_isres_decon.m
To be used on the HPC, this script runs the ISRES optimization function on the deconvolution
model (or any other model you want to feed it; use it as a template).

sub_decon.sh
This is the file that needs to be submitted as a job on the HPC to run script_isres_decon.m
correctly. The parameters are mutable. As of this writing, the command to run this on
the hpc is bsub < sub_decon.sh
