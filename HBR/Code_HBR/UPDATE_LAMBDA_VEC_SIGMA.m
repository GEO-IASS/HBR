
function [lambda_snr_vec] = UPDATE_LAMBDA_VEC_SIGMA(DATA, Model, DAG, lambda_snr_vec)


    % Calculate the lambda SNR for each response node
    % these values are saved back into lambda_snr_vec
    for child_node = 1:Model.n_resp_nodes
    
        %
        % For this operation we need the hyper-parameters alpha and beta
        %  which depend on the variance and data. 
        %  First, the variance is calculated
        %  Second, alpha_snr and Third beta_snr.
        %
        % The last step is to sample a new lambda_snr from a Gamma
        % distribution.
        
        lambda_snr  = lambda_snr_vec(child_node,1);
    
        % extract parents
        parents  = find(DAG(:,child_node) > 0);
    
        % hyper-parameters
        alpha_sigma = Model.nue_var/2;
        beta_sigma  = Model.nue_var/2;

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

        
        % extract data, X : design matrix, y : response vector
        X = [ones(1,n_obs); pred_data(parents,:)];  % #pred x #obs, first '1' row is the intercept 
        y = DATA.Gradients{child_node}';                % #obsx1

        % The mean prior vector. In HBR it is all zero. One element for each
        % parent and one for the intercept
        mue_prior = zeros(length(parents) + 1, 1);
        m_tilde = X'*mue_prior; 
                         
        % #obs x #obs                  
        Identity = eye(length(parents) + 1);
        inv_Sigma_tilde = eye(n_obs) - X'*inv(inv(lambda_snr * Identity) + X*X')*X;
                 
        Delta2 = (y - m_tilde)'*inv_Sigma_tilde*(y - m_tilde); 
                
        alpha_sigma = alpha_sigma + n_obs/2;
        beta_sigma  = beta_sigma  + Delta2/2;
   
        % sample inverse of variance with Gamma 
        inv_var_all = gamrnd(alpha_sigma,(1/beta_sigma));
        var_all     = 1/inv_var_all;
        
        % for each component, we would calculate the lambda SNR
        K_snr        = 1;  % <-- this is 1 because we only have a single component in the homogeneous HBR
        card_parents = length(parents)+1;

        % get alpha and beta SNR base values
        alpha_snr  = Model.alpha_snr + (K_snr*card_parents)/2;
                        
        % #pred #pred              
        Sigma_inv = inv(lambda_snr * Identity ) + X*X';                        
        mue_apost = inv(Sigma_inv)*(inv(lambda_snr * Identity )*mue_prior+X*y);

        % multivariate normal sample with mean posterior and new covariance
        W_i = mvnrnd(mue_apost,var_all*inv(Sigma_inv));
        W_i = W_i';

        
        % Get beta_snr
        beta_snr = Model.beta_snr + 0.5 * inv_var_all * (W_i - mue_prior)' * inv(Identity)*(W_i - mue_prior); 
       
        % Finally:
        % Sample inverse lambda SNR from Gamma
        inv_lambda_snr           = gamrnd(alpha_snr,(1/beta_snr));
        lambda_snr_vec(child_node,1) = (1/inv_lambda_snr);

    end  % loop over nodes

return
         
