%  HBR is available for use under the following license, commonly known
%  as the 3-clause (or "modified") BSD license: 
%
%  Copyright (c) 2013-2016, Marco Gzegorczyk and Andrej Aderhold
%
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without 
%  modification, are permitted provided that the following conditions 
%  are met:
%
%  1. Redistributions of source code must retain the above copyright 
%     notice, this list of conditions and the following disclaimer.
%  2. Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the 
%     documentation and/or other materials provided with the distribution.
%  3. Neither the name of the copyright holder nor the names of its contributors 
%     may be used to endorse or promote products derived from this software 
%     without specific prior written permission.
% 
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
%  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
%  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
%  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
%  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
%  OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
%  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%



% Function: main_HBR()
%
% Main method of HBR that takes the following arguments: 
%
%  Time_series  : A cell structure, each containing a time-series matrix or 
%   a single time-series matrix. In both cases, the matrix has 'nr_predictors-by-nr_samples' 
%   dimensions. I.e. each row corresponds to the time-series of one predictor variable, and
%   each column is one data time-point.  
% 
%   The cell structure is used when each response variable shall have a specific 
%   time-series predictor matrix. The predictor time-series for a response 'i' should 
%   be placed into Time_series{i} = time_series_matrix.
%
%   In the case that a single predictor matrix exists for all response
%   variables, the matrix should be placed directly into 'Time_series' with
%   Time_series = time_series_matrix. main_HBR() detects what kind of setup
%   is used. 
%
%  Gradients  :  A cell structure that contains the gradient vector for each
%   response as one element. The gradient for a response 'i' should be placed 
%   into Gradients{i} = gradient_vector. The vector has '1-by-nr_samples'
%   dimensions. 
%
%  Iterations  : Integer value with the number of MCMC iterations. 
%
%  DAG         : A 'nr_predictors-by-nr_responses' matrix that contains the
%   inital network stucture as a directed acyclic graph (DAG). Each column
%   corresponds to the incoming edges of one response variable. Defining this 
%   matrix allows to define certain edges as always on and fixed (value 2), 
%   or to be ignored (value -1). Edges with value 0 (edge is off) or 1 (edge exists) 
%   are updated by HBR in each iteration of the MCMC. 
%   Note: Leave the matrix empty with '[]' in order to generate a generic matrix.  
%
% 
% Returns: 
%  
%  Run  :  Contains samples of the MCMC simulation with the following
%  members: 
%
%  Run.dag         : a cell structure that contains the DAG of each iteration.
%  Run.Log_Scores  : a vector with the log scores in each iteration.
%  Run.steps       : integer value with last MCMC iteration.
%  Run.lambda_snr_vec : a cell structure with lambda SNR vectors.
%
% Note: Samples start after the burn-in phase that finishes at Iterations/2 iterations. 
%       The DAG samples in Run.dag can be converted to edge probabilities
%       with the function DAGs_to_edge_probabilities(Run.dag)
%

function [Run] = main_HBR(Time_series, Gradients, Iterations, DAG)

    % Setup random number seed and add script path
    rng('shuffle')
        
    addpath('Code_HBR')
   
    % This structure is passed around the sub-routines
    DATA_ALL.Time_series = Time_series;
    DATA_ALL.Gradients = Gradients;
    
    
    %
    % Define HBR model parameters
    %
    Model.mcmc_iterations = Iterations; 
    Model.burnin = Model.mcmc_iterations / 2;
    Model.p_transition = 0.2;
    Model.k_transition = 1;
    Model.max_fanin = 4;
    Model.alpha_snr = 2;
    Model.beta_snr = 0.2;
    Model.lambda_snr = 1;
    Model.nue_var = 0.01;
    Model.k_max = 10;
    
    % Check if the input time-series is composes of several matrices, one
    % for each response. If 'DATA_ALL.Time_series' is a cell structure, and
    % not a matrix, then this is true. Otherwise 'DATA_ALL.Time_series'
    % will be a single matrix with the time-series data.
    Model.multi_timeseries = 0;
    if iscell(DATA_ALL.Time_series)
        Model.multi_timeseries = 1;
        
        % get number of predictor variables from first time-series in the
        % cell structure
        Model.n_pred_nodes = size(DATA_ALL.Time_series{1},1);
    else
        % get number of predictor variables from the matrix itself.
        Model.n_pred_nodes = size(DATA_ALL.Time_series,1);
    end
                
    % set number of response and predictor nodes
    Model.n_resp_nodes = length(DATA_ALL.Gradients);
    

    % Define the prior probability for the number of parents. 
    % If everything is set to zero, we have an equal fan-in preference.
    Model.fanin_prior = zeros(Model.n_pred_nodes + 1, 1);   
   
  
    % If no pre-defined network graph (DAG) was passed, i.e. the matrix 
    % DAG was empty ([]), then create a generic DAG based solely on the 
    % number of response and predictor variables.
    %
    % The DAG matrix has 'nr_predictors-by-nr_responses' dimensions, i.e. 
    %  each column corresponds to the edges of one response. 
    %  If an edge is added, removed, or flipped then the DAG is updated 
    %  accordingly. 
    %
    % The elements in the DAG can have following values
    % 
    %  -1 : ignore this parent - it is never selected or unselected
    %   0 : potential parent that can become a selected parent ('1')
    %   1 : selected parent that can become a potential parent ('0')
    %   2 : selected parent that is fixed and is never changed
    %
    % * Added edges become '1'.
    % * Removed edges become '0'.
    % * Flipped edges change from '1'->'0' or '0'->'1' 
    %
    if isempty(DAG)
        
        % Initialize the generic DAG with zeros, i.e. all parents (edges)
        % can become potential predictors. 
        DAG = zeros(Model.n_pred_nodes, Model.n_resp_nodes);
        
        % If the number of responses and predictor are the same, which is
        % often the case when all nodes in a network are potential
        % responses as well as predictors, then set the diagonal element to
        % the flag '2'. This marks the response it self as a fixed parent.
        %
        % Set it to '1' to mark it as an edge that can be removed.
        %
        % Set it to '-1' to disable the self-loop. 
        %
        % Disable these lines if the self-loop should be treated as all
        % other potential parent. 
        %
        if Model.n_pred_nodes == Model.n_resp_nodes
            DAG(1:(Model.n_resp_nodes+1):end) = 2;
        end
        
        % Note, if the predictor set differs from the response set you need
        % to customize the DAG if you like certain parents to be disabled
        % (use '-1') or fixed predictors (use '2'). 
    end
       
    fprintf('\n######################################################################################\n')
    fprintf('The new HBR MCMC simulation has been started, Total number of iterations: %i \n', Model.mcmc_iterations);
    fprintf('######################################################################################\n\n')
 
    %
    %
    % Start of the MCMC simulation
    %
    %  This returns the structure 'Run' that has following members:
    %
    %  Run.dag            = {};
    %  Run.Log_Scores     = [];
    %  Run.steps          = 0;
    %  Run.lambda_snr_vec = {};
    % 
    % The relevant data is Run.dag - it contains the DAG for each
    % iteration after the burn-in defined by Model.burnin.
    % Use these DAGs to calculate the parent posterior probabilities by
    % averaging over all DAGs found in this structure. 
    %
    % See the script Evaluate_Biopepa_full.m or the end of script 
    %                Example_Biopepa_simple.m for an example.
    %
    Run = MCMC_HBR(DATA_ALL, Model, DAG);

     


return;



