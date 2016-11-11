%
% Main MCMC code for HBR.
%
%


function [Run] = MCMC_HBR(DATA_ALL, Model, DAG)

    % Initialize result data structure. Values are saved after Model.burnin
    % iterations are passed.
    % MCMC run
    Run.dag            = {};   % all DAGs 
    Run.Log_Scores     = [];   % Log scores
    Run.steps          = 0;    % last iteration step 
    Run.lambda_snr_vec = {};   % vector of lambda SNR values
    
    % Calculate the lambda SNR vector the first time, this vector is updated 
    % in each MCMC iterations
    lambda_snr_vec  = Model.lambda_snr  * ones(Model.n_resp_nodes,1);
   
    
    %
    %
    %  Main HBR MCMC Loop 
    %
    %
    for t = 1:Model.mcmc_iterations

        
        % Loop over each response (child) node
        for child_node = 1:Model.n_resp_nodes


            % If chance lower 0.5, attempt addition/deletion move
            if (rand(1) < 0.5) 

                % extract current parents of this child, including
                % self-loop, ignore parents marked with '-1'
                existing_parents = find(DAG(:,child_node) > 0);
               
                % check if we can add another parent, the constraint is max_fanin
                if (length(existing_parents) < Model.max_fanin);

                    unselected_parents = find(DAG(:,child_node) == 0);
                    parent_change_index = unselected_parents(randsample(length(unselected_parents), 1));
                    
                    % number of potential parents, excluding fixed parents
                    % such as self-loop, goes into the Hastings-factor
                    neighbours_old = sum(DAG(:,child_node) == 0) + sum(DAG(:,child_node) == 1); 
                    
                else % max max_fanin was reached, attempt to delete a parent

                    % extract parent, but not the self-loop (these have flag '2')
                    selected_parents = find(DAG(:,child_node) == 1);
                    
                    % sample the parents to delete
                    parent_change_index = selected_parents(randsample(length(selected_parents), 1));
                                      
                    % set old neighbours to the max_fanin
                    neighbours_old = Model.max_fanin; 
                    
                end
 
               
                % update the DAG 
                DAG_NEW = DAG;
                DAG_NEW(parent_change_index,child_node) = 1 - DAG_NEW(parent_change_index,child_node);

                % set neighbour number for the new DAG, that goes into the
                % Hastings-factor of the acceptance probability
                neighbours_new = Model.max_fanin;   % set to max_fanin 
                
                % But if number of parents in the new DAG is smaller than the
                % fan-in,
                if sum(DAG_NEW(:,child_node) > 0) < Model.max_fanin
                
                    % .. set to number of potential parents
                    neighbours_new = sum(DAG_NEW(:,child_node) == 0) + sum(DAG_NEW(:,child_node) == 1); 
                end
                
                  
                % compute local log score for new and old parent
                % configuration. "Local" means the score only for a single 
                % response node
                local_score_new = COMPUTE_LOCAL_LOG_SCORE(DATA_ALL, Model, DAG_NEW, child_node, lambda_snr_vec);
                local_score_old = COMPUTE_LOCAL_LOG_SCORE(DATA_ALL, Model, DAG,     child_node, lambda_snr_vec);

                % compute acceptance probability
                A = exp(local_score_new - local_score_old + (log(neighbours_old) - log(neighbours_new) ) );

                % if acceptance 'A' is higher than random uniform sample in [0,1]
                % accept the new DAG
                if (rand(1) < A) 
                    DAG = DAG_NEW; 
                end

            else % Perform an exchange move
                
                % extract existing parents, excluding self-loop, which has
                % a '2' flag in the DAG
                existing_parents = find(DAG(:,child_node) == 1);
                
                % Check if there is at least one parent
                if (length(existing_parents) >= 1)
                % perform an exchange move

                    DAG_NEW = DAG;
                    
                    % Sample one parent from the existing ones
                    %
                    % Note: Use the more complicated form of randsample() below,
                    %       instead of just calling randsample with
                    %
                    %    parent_old_index = randsample(existing_parents, 1);
                    %
                    % In case 'existing_parents' has a length of 1, it will
                    % not return a sample from this 1-element vector but 
                    % a sample from 1:(value of element). 
                    %
                    parent_old_index = existing_parents(randsample(length(existing_parents), 1));
            
                    
                    % sample new parent from potential parents
                    unselected_parents = find(DAG(:,child_node) == 0);
                    parent_new_index = unselected_parents(randsample(length(unselected_parents), 1));
                    
                    % setup the new DAG: disable the previously selected
                    % and existing parent, and enable the new non-existing
                    % parent
                    DAG_NEW(parent_old_index,child_node) = 0;
                    DAG_NEW(parent_new_index,child_node) = 1;

                    % compute local score
                    local_score_new = COMPUTE_LOCAL_LOG_SCORE(DATA_ALL, Model, DAG_NEW, child_node, lambda_snr_vec);
                    local_score_old = COMPUTE_LOCAL_LOG_SCORE(DATA_ALL, Model, DAG,     child_node, lambda_snr_vec);

                    log_hastings = 0;

                    A = exp(local_score_new - local_score_old + log_hastings);

                    % if acceptance 'A' is higher than random uniform sample in [0,1]
                    % accept the new DAG
                    if (rand(1) < A) 
                        DAG = DAG_NEW; 
                    end

                end    
            end  % End of if clausel, addition/deletion move      
        end  % End loop over nodes


        % Update lambda_snr_vec and the log score
        lambda_snr_vec  = UPDATE_LAMBDA_VEC_SIGMA(DATA_ALL, Model, DAG, lambda_snr_vec);
        log_score = COMPUTE_LOG_SCORE(DATA_ALL, Model, DAG, lambda_snr_vec);
        
        
        % Save some parameters after a burnin of half simulation length
        if (t >= Model.burnin)
            
            Run.dag{end+1}              = DAG;
            Run.Log_Scores(end+1)       = log_score;
            Run.steps                   = t;
            Run.lambda_snr_vec{end+1}   = lambda_snr_vec;
            
        end

     
    end % End of main loop

    
return
