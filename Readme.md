

This is the Matlab code for the software "Hierachical Bayesian Regression" (HBR) and analytic gradient calculation (GCGP) published in 

    [1] Aderhold, A., Husmeier, D., & Grzegorczyk, M. (2014). Statistical inference of regulatory networks for circadian regulation. Statistical applications in genetics and molecular biology, 13(3), 227-273.

The software contains two components: The HBR code and the gradient calculation code (GCGP) that uses the Gaussian Process package GPstuff (http://research.cs.aalto.fi/pml/software/gpstuff/) by 

    [2] Jarno Vanhatalo, Jaakko Riihimäki, Jouni Hartikainen, Pasi Jylänki, Ville Tolvanen, and Aki Vehtari (2013). GPstuff: Bayesian Modeling with Gaussian Processes. Journal of Machine Learning Research, 14(Apr):1175-1179

The HBR code was originally written by Marco Grzegorczyk, and adapted by Andrej Aderhold for the Biopepa data. It is located in the subfolder 'HBR/', see Section A. for details.
The GCGP code was written by Andrej Aderhold, and uses sub-routines from the GPstuff, version 4.4 software that is under the GNU GPL license. It is located in the subfolder GCGP, see Section B. for details.

This software is also used in 

    [3] Grzegorczyk, M., Aderhold, A., & Husmeier, D. (2015). Inferring bi-directional interactions between circadian clock genes and metabolism with model ensembles. Statistical applications in genetics and molecular biology, 14(2), 143-167.

    [4] Aderhold, A., Grzegorczyk, M., and Husmeier, D. (2016). Approximate Bayesian inference in semi-mechanistic models. Statistics and Computing, 1-38



Usage
----------------------------------------------

Note: We recommend consulting the example scripts for the gradient calculation located in 

    GCGP/Example_process_Biopepa_Data.m

and for the main method HBR:

    HBR/Example_Biopepa_full.m
    HBR/Example_Biopepa_simple.m


Prior to running the HBR software, the gradients need to be generated with GCGP. Change into the directory 'GCGP' and 
examine the script 'main_GCGP.m' that contains the folloging function header:

    function [All_gradients] = main_GCGP(ts_data, ts_data_changepoints, gradient_type) 

This function takes three arguments: The time-series data as a matrix ('ts_data'), a vector that defines the 
change-points in the data in the case the time-series is segmented ('ts_data_changepoints'), and the gradient type 
that should be generated ('gradient_type'). This can be 'RBF' or 'PER' for the analytic RBF and periodic kernel 
using the Gaussian Process, or 'coarse' for a numeric calculation using the difference quotient. 

As an example consider the following code line from the file 'GCGP/Example_process_Biopepa_Data.m':

   gradient = main_GCGP(Data_concentrations(mRNA_response_set,:), change_points, gradient_type);
 
The first argument is a matrix that contains a time-series for specific response variables for which we'd like to 
calculate a gradient. The second argument 'change_points' is a vector with indices that define the boundaries of 
the segments contained in the time-series. This is important if the time-series is an assembly of multiple experiments. 
The last argument 'gradient_type' was previously defined as a string with 'RBF', 'PER', or 'coarse'. 
Finally the return variable 'gradient' contains the matrix of gradients for all variables passed in the first argument. 
Note, the gradient matrix is transposed and then saved to a .csv file to be read in by HBR:

    % Save the response gradients to this file
    outfile = sprintf('Processed_Data_CSV/%s_%sGradient_id%i.csv', file_name_out, gradient_type, dataid);

    % Write with header
    dlmwrite(outfile, head_Gradient, '');
    dlmwrite(outfile, gradient, '-append')

This .csv file contain a header defined with 'head_Gradient', see 'GCGP/Example_process_Biopepa_Data.m' on how 
to define this array of strings. Note that the time-series data is it-self pre-processed inside 
'GCGP/Example_process_Biopepa_Data.m' by a z-score transformation of the concentration data. 
One could move this code into the main 'HBR' code if preferred. 

The main program 'HBR' is implemented in the script 'HBR/main_HBR.m'. The header of this function is defined as:

    function [Run] = main_HBR(Time_series, Gradients, Iterations, DAG)

This script runs the main procedure for one data set defined in 'Time_series' with the gradient data located in 'Gradients'. 
Two example scripts are provided that call this function for a more complex and a simple data set:

    HBR/Example_Biopepa_full.m
    HBR/Example_Biopepa_simple.m

Note, that both scripts make use of data pre-processing that involves the deletion of specific observations 
dependent on the response variable. Thus, the time series passed to 'main_HBR()' is a cell-structure containing 
multiple matrices as explained in more detail in Section A.1. This can be simplified by using a single time-series matrix. 
The example script will execute the main method with:  

    Run = main_HBR(Time_series, Gradients, Iterations, DAG);

The returned variable 'Run' contains all information necessary to calculate the edge posterior probabilities:
    
    % Convert to edge probilities
    Edge_matrix = DAGs_to_edge_probabilities(Run.dag);
     
Finally 'Edge_matrix' can be used to calculate a score if a gold standard is provided. 
For instance, after converting Edge_matrix to a vector 'Edge_vec', and  providing an equivalent 
vector with gold standard information 'realnet_vec', the AUROC and AUPREC scores can be calculated with:

    [auroc, auprec] = AUROC_AUPREC_scores(realnet_vec, Edge_vec);

See the example scripts for more details and comments on how to execute HBR. 



A.1 Hierarchical Bayesian Regression - HBR
----------------------------------------------

The HBR software samples parent sets for response variables given a time-series of potential predictor variables 
and a gradient vector for each response variable. The sample of the resulting network of parent-response 
configurations is saved for each iteration of the MCMC chain in the form of a directed acyclic graph (DAG). 
These graph samples can be used to calculate the edge probabilities of the full network. 
The software directory is 'HBR/' and contains the follwing sub-directories and scripts:

     HBR/main_HBR.m                            The main software that samples predictor-response networks. 
        /Example_Biopepa_simple.m	       A single example network calculation from [1] and [3]. This scripts calls main_HBR(). 
					       It also converts the DAG to edge probabilities and calculates AUROC and AUPREC scores 
					       using DAGs_to_edge_probabilities() and AUROC_AUPREC_scores().
	/Example_Biopepa_full.m		       The full Biopepa data calculation for all networks and gradient settings from [1] and [3].
	/Evaluate_Biopepa_full.m 	       Evaluation of the results generated by Example_Biopepa_full.m.
	/Results/			       Directory containing the results of the example scripts.
	/Data_Biopepa/			       Directory containing the Biopepa data.
	/Data_Biopepa/Goldstandards	       The true network, i.e. the gold standard, of the Biopepa example data.
	/Code_HBR/			       Directory containing the main code of HBR.
	         /MCMC_HBR.m		       The main MCMC loop.
		 /COMPUTE_LOCAL_LOG_SCORE.m    Computes the log score for a single reponse variable.
		 /COMPUTE_LOG_SCORE.m	       Computes the log score for all response variables.
		 /UPDATE_LAMBDA_VEC_SIGMA.m    Computes the lambda SNR vector and the sigma variable.
		 /DAGs_to_edge_probabilities.m Converts multiple DAG samples returned by the main routine into edge probabilities.
		 /AUROC_AUPREC_scores.m	       Computes the area under the receiver operating characteristic (AUROC) and precison/recall 
                                               curves (AUPREC) given a gold standard.		     

Run the example scripts and see the comments inside the scripts for practical application of the main script main_HBR. 
The main method is defined with 

    function [Run] = main_HBR(Time_series, Gradients, Iterations, DAG) 

located in main_HBR.m. This function takes the following four arguments:

    - Time_series  : A cell structure that contains a time-series in each element, or a single time-series matrix. 
                     In the former case, each element in the cell structure is a matrix that contains the time-series 
                     data for each response variable. The predictor time-series for a response 'i' should then be 
                     located in 'Time_series{i}. This might be necessary in the case the observations differ between 
                     response variables, e.g. because of disabled observations (see 'HBR/Example_Biopepa_simple.' as an example). 
                     In the latter case, 'Time_series' is a single time-series matrix used for all response variables. 
                     A time-series matrix always has a dimension with variables/features in the rows, and observations/samples 
                     in the columns, i.e. each row corresponds to the time-series of one predictor variable, and each column 
                     is one observation or data time-point. Note: the user does not have to specify which kind of variable 
                     (cell structure or single matrix) was used - The function main_HBR() will examine 'Time_series' and detect what was used. 

    - Gradients   :  A cell structure that contains a gradient vector for each response as one element. The gradient for a response 'i' should be placed 
                    into Gradients{i} = gradient_vector. The vector has a size of the number of samples in the data. 

    - Iterations  : Integer value with the number of MCMC iterations. 

    - DAG         : A matrix that contains the inital network stucture as a directed acyclic graph (DAG). Each column
                    defines the parent nodes of a response variable, and the rows are the predictor variables. 
                    A network can be pre-defined by setting edges of interest as always on and fixed (value 2), or to
                    be ignored (value -1). Edges with the value 0 (edge is off) or 1 (edge exists) are updated by HBR 
                    in each iteration of the MCMC. Note: Leave the matrix empty with '[]' in order to generate a generic matrix.  

 
 Returns: 
  
   + Run  :  Contains samples of the MCMC simulation with the following members: 

      - Run.dag         : a cell structure that contains the DAG of each iteration.
      - Run.Log_Scores  : a vector with the log scores in each iteration.
      - Run.steps       : integer value with last MCMC iteration.
      - Run.lambda_snr_vec : a cell structure with lambda SNR vectors.

  Note: Samples start after the burn-in phase that finishes at Iterations/2 iterations. 
       The DAG samples in Run.dag can be converted to edge probabilities with the function DAGs_to_edge_probabilities(Run.dag)


The input argument 'Gradients' can be taken from the output of the GPCG method defined in Section B. The gradient can also come from another source.


A.2. License
-------------------------------

HBR is available for use under the following license, commonly known
as the 3-clause (or "modified") BSD license:

Copyright (c) 2013-2016 Marco Gzegorczyk and Andrej Aderhold

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. The name of the author may not be used to endorse or promote products
   derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



B.1 Gradient Calculation - GCGP
--------------------------------

The GCGP software is located in the directory GCGP/. The main script is 'main_GCGP.m'. 
The contained function main_GCGP() takes a time-series as an input that can contain any 
type of concentration measurements over time. In [1] and [3], these measurements were synthetically 
generated mRNA samples from an artificial circadian clock model of the plant Arabidopsis thaliana. 
GCGP creates analytic solution for the gradients of each observation in the time-series that is derived 
from a Gaussian Process with RBF kernel. GCGP can also calculate a numeric solution using a difference quotient. 
The output of GCGP is a .csv file that contains the gradients in a table format that can be imported into the HBR software. 
The files contained GCGP are

     GCGP/main_GCGP.m                            The main software that calculates the analytic solution of a time-series gradient using a Gaussian Process.
	 /Example_process_Biopepa_Data.m         Processes the Biopepa raw data into data for HBR and calls main_GCGP(). 
	 /arabidopsis_clock.biopepa              The Biopepa circadian clock model script. 
	 /Biopepa_Protocol_data_generation.txt   A description of the Biopepa data. 
	 /Raw_Biopepa_Data/                      The input raw time-series data for the Biopepa example data is located here.
	 /Processed_Data_CSV/                    The output of the example data generated in Example_process_Biopepa_Data.m is written to this directory.
	 /GPstuff-4.4				   Version 4.4 of the Gaussian Process code of [3]. 

Run Example_process_Biopepa_Data.m to see the gradient calculation in action. Note that this script contains additional 
pre-processing steps that are specific for the Biopepa example data. 

The function main_GCGP(ts_data, ts_data_changepoints, gradient_type) takes the following arguments:
 
    - ts_data : a matrix with the variables in the rows and samples/observations in the columns, 
                i.e. the time-series of multiple variables is in this matrix. 
   
    - ts_data_changepoints : a vector with changepoints that segment ts_data into segments. 
                             A gradient calculation is executed for each segment. Important 
                             for multi-segment data, e.g. time-series that contain different experiments. 
                             Leave this vector empty if you do not desire any segmentation.
    
    - gradient_type : 'RBF' for the analytic gradient, or 'coarse' for the numeric gradient using a difference quotient. 

  Returns: 

   A matrix that contains the gradient data for each variable in the rows, i.e. the gradient vector 
   for the whole time-series (all segments) of one variable is stored in one row. 
   The columns correspond to the observations. Note that ts_data_changepoints should be an empty vector, if you do not wish to segment the data. 
   This means that the Gaussian Process is applied to the whole time-series at ones. 


The sub-function main_GCGP_core(y, x, kernel, initial_length) is called by main_GCGP() for each variable and segment in 'ts_data'. 
It has the following arguments:

    - y : input data vector of one variable
    
    - x : input x (time point) data vector. In most cases this can be just (1,2,..,nr_samples) where nr_samples is the number of samples in the data 'y'. 
    
    - kernel : 'RBF' or 'PER' for the RBF or periodic kernel. If set to 'coarse' it will use the difference quotient to approximate a gradient.
    
    - inital_length : inital value for the length hyper-parameter. If empty, it is set to 1. 

  Returns: 

   A gradient vector with the size of the input vector y. 

The script main_GCGP.m also contains the derivatives of the RBF and periodic kernel function with the names 

    rbf_kernel_deriv(x_star, x, length, sigma2)
    per_kernel_deriv(x_star, x, length, sigma2, period)

Both return the covariance matrix of the function derivative.  

Note that the fminscg() optimization routine in GPstuff was slightly modified to avoid infinite loops caused by invalid optimization values. 
GPstuff will now throw an error after 500 invalid iterations, which is caught by GPCG. 
GPCG will then adapt the initial hyper-parameter value of the Gaussian process to avoid GPstuff becoming stuck again.


B.2. License
---------------------------------------------

This software is distributed under the GNU General Public License (version 3 or later); please refer to the file License.txt, 
included with the software, for details.
