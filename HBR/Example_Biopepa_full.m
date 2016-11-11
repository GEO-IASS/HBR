

%
% Batch run for different settings of the Biopepa data from the SAGMB
% publication
%
%
function Example_Biopepa_full()

    %
    % Run HBR for the Biopepa data of the SAGMB 2014 paper
    % with different gradients and predictor types (mRNA and protein)
    %
    HBR_Biopepa_full_proc('coarse', 'mRNA')
    HBR_Biopepa_full_proc('coarse', 'protein')
    %HBR_Biopepa_full_proc('RBF', 'mRNA')
    %HBR_Biopepa_full_proc('RBF', 'protein')

end


%
% This function will run HBR for the Biopepa data for two specific settings
% - the gradient and the predictor type:
%
% The gradient can be 'coarse' or 'RBF', and the predictor
% type can be 'mRNA' or 'protein'
%
function HBR_Biopepa_full_proc(gradient, predictor_type)


    % 20000 Iterations were used in the SAGMB2014 paper, but half that number
    % gave also very similar results.
    Iterations = 10000;
   
    % directory for the input data
    data_dir = 'Data_Biopepa';
       
    % Define data related parameters
    networks = {'wildtype'}; %, 'mutant_PRR7_PRR9', 'mutant_PRR7_PRR9_TOC1', ...
               %'mutant_PRR5_PRR7_PRR9', 'mutant_TOC1', 'mutant_PRR5_PRR7_PRR9_TOC1'};

    % number of data instances for each network
    data_instances = 1:5;
      
    % all results go into here
    RESULTS = [];
   
    % The indices for the response variables
    response_vec = [2,4,8,10,12,14,17];
    
    % The indices for the mRNA data (same indices as response), in addition
    % element 7 is included corresponding to variable 'Light*P'. Light is a
    % crucial for the circadian clock and should not be missed!
    % The 'P' is a light modulating protein that causes a light peak in the
    % morning, but decreases light intensity afterwards. 
    predictor_vec = [2,4,7,8,10,12,14,17]; 

    
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
    
 
   
    % Loop over the different networks
    for network_idx = 1:length(networks)

        network = networks{network_idx};

        % Loop over different data instances of the network
        for i=data_instances

            % read the time-series concentration data
            conc_file = sprintf('%s/biopepa_%s_ts-concentrations_id%i.csv', data_dir, network, i);
            tmp = importdata(conc_file);
            TS = tmp.data';
            
            % read the gradient data
            grad_file = sprintf('%s/biopepa_%s_%sGradient_id%i.csv', data_dir, network, gradient, i);
            tmp = importdata(grad_file);
            GRAD = tmp.data';
           
            
            fprintf('*** Network: %s, data instance: %i, gradient: %s', network, i, gradient);
    
            % For each response variable there is one matrix in
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
            DATA_ALL.Time_series = {};
            DATA_ALL.Gradients = {};
            
            % Biopepa data specific pre-processing:
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
                DATA_ALL.Time_series{end+1} = TS(predictor_vec, take);
                
                % extract gradient
                DATA_ALL.Gradients{end+1} = GRAD(ridx, take); 
                
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
            n_resp_nodes = length(DATA_ALL.Gradients);
            n_pred_nodes = size(DATA_ALL.Time_series{1},1);

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
            % Note: If you set the third argument to an empty set, such as this:
            %
            %   Run = main_HBR(DATA_ALL, Iterations, []);
            %
            % then main_HBR() creates a generic DAG with all zero entries. If the
            % number of responses and predictor match, the diagonal (the
            % self-loops) are set to '2' corresponding to a fixed edge.
            %
            Run = main_HBR(DATA_ALL, Iterations, DAG);
            
            RESULTS{network_idx}{i} = Run;

            savefile = sprintf('Results/OUT_HBR_%s_%sGradient.mat', predictor_type, gradient);    
            save(savefile,'RESULTS'); 
        end        

    end


end
  






