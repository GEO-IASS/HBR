%
% Compute log score for single response node
%
function [log_score] = COMPUTE_LOCAL_LOG_SCORE(DATA, Model, DAG, child_node, lambda_snr_vec)

    %
    % get time-series for this child/response
    %
    if Model.multi_timeseries 
        pred_data = DATA.Time_series{child_node};
    else
        % for the single time-series, get the matrix directly
        pred_data = DATA.Time_series;
    end
   
    n_obs = size(pred_data, 2);
    
    % Extract log prior for the fanin
    log_prob_graph =  Model.fanin_prior(length(find(DAG(:,child_node) > 0)) + 1);

    %
    % Compute the local score for the response child_node:
    %
    parents = find(DAG(:,child_node) > 0);

    lambda  = lambda_snr_vec(child_node,1);

    sum_log_det_Sigma_tilde = 0;
    sum_Delta2              = 0;

    
    % extract data for log likelihood calculation
    X = [ones(1,n_obs); pred_data(parents,:)];  % #pred x #obs, first '1' row is the intercept 
    y = DATA.Gradients{child_node}';                % #obsx1

    % The mean prior vector. In HBR it is all zero. One element for each
    % parent and one for the intercept
    mue_prior = zeros(length(parents) + 1, 1);
    m_tilde     = X'*mue_prior;   

    % #pred x #obs   
    Identity = eye(length(parents) + 1);
    Sigma_tilde = eye(n_obs) + lambda*X' * Identity * X;

    % #pred x #obs
    inv_Sigma_tilde = eye(n_obs) - X'*inv(inv(lambda * Identity) +X*X')*X;

    % (1x#obs) * (#obsx#obs) * (#obsx1) 
    sum_Delta2 = sum_Delta2 + (y-m_tilde)'*inv_Sigma_tilde*(y-m_tilde); 

    sum_log_det_Sigma_tilde = sum_log_det_Sigma_tilde + log(det(Sigma_tilde));                                   

    sum_1 = gammaln((n_obs + Model.nue_var)/2) - gammaln(Model.nue_var / 2);

    sum_2 = (Model.nue_var / 2)*log(Model.nue_var) - (n_obs / 2) * log(pi) - 0.5 * sum_log_det_Sigma_tilde;

    sum_3 = -(n_obs + Model.nue_var) / 2 * log(Model.nue_var + sum_Delta2);

    log_score_i = sum_1 + sum_2 + sum_3;

    log_score = log_prob_graph + log_score_i;

return
