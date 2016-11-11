%
%
% This will run HBR for a single data set of the Biopepa data, and 
%  evaluate the result. Part of the evaluation is
%    
%   * calculate edge probabilities
%   * calculate AUROC and AUPREC scores given a gold standard
% 
%
%

function Example_Biopepa_simple()

    Iterations = 1000;
   
    % gradient can be 'RBF', 'PER', or 'coarse'
    gradient = 'RBF';
    
    % directory for the input data
    data_dir = 'Data_Biopepa';
       
    % The prediction type can be 'mRNA' or 'protein'.
    % This controls which variables can become potential parents, and
    % the structure of the DAG. 
    predictor_type = 'mRNA';
    
    % The indices for the response variables.
    % The Biopepa data set is a mix of mRNA and protein data points
    % The following indices correspond to the mRNA data points.
    response_vec = [2,4,8,10,12,14,17];
    
    % The indices for the predictor variables, which correspond to the mRNA 
    % data points. 
    % In addition, element 7 is included corresponding to variable 'Light*P'. 
    % Light is crucial for the circadian clock and should not be missed!
    % The 'P' is a light modulating protein that causes a light peak in the
    % morning, but decreases light intensity afterwards. 
    predictor_vec = [2,4,7,8,10,12,14,17]; 

    % In case the prediction types are proteins, change the predictor
    % indices.
    if strcmp(predictor_type, 'protein') 
        % 
        % The protein indices are actually [3,5,6,7,9,11,13,15,16,18]
        % However, we include the mRNA self-loop concentration for each 
        % response as fixed parent, hence we have all the indices from 1 to
        % 18.
        % The DAG matrix defines, which indices are actually potential
        % parents, and which will be ignored or are fixed.

        predictor_vec = 1:18; 

    end
    

    % read the time-series concentration data
    conc_file = sprintf('%s/biopepa_wildtype_ts-concentrations_id1.csv', data_dir);
    tmp = importdata(conc_file);
    TS = tmp.data';
            
    % read the gradient data
    grad_file = sprintf('%s/biopepa_wildtype_%sGradient_id1.csv', data_dir, gradient);
    tmp = importdata(grad_file);
    GRAD = tmp.data';
           
            
    fprintf('*** Network: wildtype, data instance: 1, gradient: %s ', gradient);
    
    %For each response variable there is one matrix in
    % 'DATA_ALL.Time_series' and 'DATA_ALL.Gradients'. 
    %
    % DATA_ALL.Time_series can also be a single time-series matrix.
    %
    % Having a cell structure with one matrix per response, allows
    % us to change the number of observations for each response,
    % which we need for Biopepa. 
    %
    % Note: main_HBR() detects if there is a single time-series or
    % multiple time-series in 'DATA_ALL.Time_series'.
    Time_series = {};
    Gradients = {};

    % Biopepa data specific pre-processing:
    %
    %    
    % For some of the response variables, specific observations
    % are excluded. These are observations, which are known to 
    % be knocked out and should not have an impact on network
    % inference. However, we showed previously that network 
    % inference deteriorates whenever these specific observations
    % are kept in the time-series.
    %
    % Here we loop over each response and take out some
    % observations. In addition, each response is assigned it's own design
    % matrix with 'DATA_ALL.Time_series{response}' since the dimensions can change. 
    %
    for ridx = 1:length(response_vec) 

        response = response_vec(ridx);

        exclude_obs = [];

        if(response==2)  exclude_obs = 1:13;      end
        if(response==4)  exclude_obs = 14:26;     end
        if(response==14) exclude_obs = 27:39;     end
        if(response==10) exclude_obs = 40:52;     end
        if(response==12) exclude_obs = 40:52;     end

        nrobs = size(TS, 2);

        take = setdiff(1:nrobs, exclude_obs);

        % extract relevant concentration data
        Time_series{end+1} = TS(predictor_vec, take);

        % extract gradient
        Gradients{end+1} = GRAD(ridx, take); 

    end

    
    
    %
    % Initialize the DAG and start MCMC
    %
    
    % Setup the network graph as DAG -- predictor-by-response, i.e. 
    %  each column corresponds to the edges for one response
    %
    % The elements in the DAG can have following values
    % 
    %  -1 : ignore this parent - it is never selected or unselected
    %   0 : potential parent that can become a selected parent ('1')
    %   1 : selected parent that can become a potential parent ('0')
    %   2 : selected parent that is fixed and is never changed
    %
    n_resp_nodes = length(Gradients);
    n_pred_nodes = size(Time_series{1},1);

    if strcmp(predictor_type, 'mRNA')
        %
        % The DAG for mRNA contains potential parents ('0')  
        % and some fixed selected parents that act as self-loops ('2')
        %
        DAG = zeros(n_pred_nodes, n_resp_nodes);

        % Define the self-loop indices
        % Note, here index 3 is missing. This is because 3 is associated to 
        % the 'Light' variable in the Biopepa data and it is a potential 
        % parent ('0')
        SELF_LOOPS = [1,2,4,5,6,7,8];
        
         % set the self-loop for each response
        for i_node = 1:n_resp_nodes
            DAG(SELF_LOOPS(i_node), i_node) = 2;
        end
    else
        %
        % The DAG for the protein predictors is more complicated!
        %
        %  Each response node also has its mRNA as a fixed self-loop (marked
        %  with '2') in the DAG. This is the same as with the "mRNA"
        %  predictor type.
        %
        %  In addition, the proteins that act as potential parents can
        %  differ for different responses. In fact, the protein associated 
        %  to each response is excluded from the list of potential parents.
        %  Some responses even have two proteins -- the basic one and a 
        %  modified version.
        %
        %  The cell structure 'POTENTIAL_PARENTS_VEC' lists the potential 
        %  proteins for each response (marked with '0'), and 'SELF_LOOPS' 
        %  lists the mRNA self-loop (marked with '2'). 
        %
        %  The DAG is initially fully defined with '-1' values. Elements 
        %  with this value are always ignored in the upcoming edge moves.
        %
        
        DAG = zeros(n_pred_nodes, n_resp_nodes) - 1;

        POTENTIAL_PARENTS_VEC{1}  = [  5,6,7,9,11,13,15,16,18];
        POTENTIAL_PARENTS_VEC{2}  = [3,    7,9,11,13,15,16,18]; 
        POTENTIAL_PARENTS_VEC{3}  = [3,5,6,7,  11,13,15,16,18]; 
        POTENTIAL_PARENTS_VEC{4}  = [3,5,6,7,9,   13,15,16,18];
        POTENTIAL_PARENTS_VEC{5}  = [3,5,6,7,9,11   ,15,16,18];
        POTENTIAL_PARENTS_VEC{6}  = [3,5,6,7,9,11,13,      18]; 
        POTENTIAL_PARENTS_VEC{7}  = [3,5,6,7,9,11,13,15,16   ];
                
        SELF_LOOPS = [2,4,8,10,12,14,17];
        
        % set the potential parents and the self-loop for each response
        for i_node = 1:n_resp_nodes
            DAG(POTENTIAL_PARENTS_VEC{i_node}, i_node) = 0;
            DAG(SELF_LOOPS(i_node), i_node) = 2;
        end
        
    end
    
  

    % Executes the main procedure of HBR.
    %
    % Note: If you set the last argument to an empty vector, such as this:
    %
    %   Run = main_HBR(Time_series, Gradients, Iterations, []);
    %
    % then main_HBR() creates a generic DAG with all zero entries. If the
    % number of responses and predictor match, the diagonal (the
    % self-loops) are set to '2', which corresponds to a fixed edge.
    %
    Run = main_HBR(Time_series, Gradients, Iterations, DAG);
    
    % Convert to edge probilities
    Edge_matrix = DAGs_to_edge_probabilities(Run.dag);
    
 
    %
    % Evaluate the results, i.e. the Edge_matrix
    %
    %  1. Read in the gold standard matrix and convert to a vector
    %  2. Convert the edge probabilities to a vector
    %  3. Execute the AUROC_AUPREC_scores.m script that calculates these
    %  scores based on the 1. and 2. 
    %
    
    %
    % 1. Read gold standard for this network and convert to vector
    %
    goldfile = sprintf('Data_Biopepa/Goldstandards/true_network_wildtype_%s.txt', predictor_type);
    fprintf('read %s\n', goldfile);
    realnet = importdata(goldfile);  % reads in response-by-predictor matrix
    
    % transpose this matrix so it is in the same format as 'Edge_matrix', i.e. predictor-by-response 
    realnet = realnet';  
    
    % convert to vector and remove entries that should be ignored (flagged
    % with '-1')
    realnet_vec = realnet(:);
    realnet_vec = realnet_vec(realnet_vec ~= -1);
    
    
    %
    % 2. Convert to edge probabilities to vector and remove invalid
    % elements. Invalid elements are marked with NaN, see the script 
    % DAGs_to_edge_probabilities.m for the marking rule.
    %
    Edge_vec = Edge_matrix(:);
    Edge_vec = Edge_vec(~isnan(Edge_vec));
    
    %
    % 3. Run the score calculation
    %
    
    [auroc, auprec] = AUROC_AUPREC_scores(realnet_vec, Edge_vec);
            
    fprintf(' *   network: wildtype, data-id: 1, AUROC: %g, AUPREC: %g\n', auroc, auprec);

           
    
end
  
