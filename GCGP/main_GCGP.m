%
%  Copyright (C) 2013-2016 Andrej Aderhold
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <http://www.gnu.org/licenses/>.


%
%
% This file contains the following functions:
%
%   main_GCGP()       : main GCGP method that takes a time-series with
%   multiple variables and segments in the time-series as input.
%
%   main_GCGP_core()  : called by main_GCGP(), takes a single variable
%   time-series vector as input and calculated the gradient. 
% 
%   rbf_kernel_deriv() : derivative of the RBF kernel, used by main_GCGP_core()
%   per_kernel_deriv() : derivative of the Periodic kernel, used by main_GCGP_core() 
%



%
% Function main_GCGP() 
%
% Arguments:
%
%  ts_data : a 'nr_variables-by-nr_samples' matrix, i.e. the time-series of 
%            multiple variables is in this matrix. Each row corresponds to a 
%            variable, and each column to a data sample.
%
%  ts_data_changepoints : a vector with changepoints that segment ts_data 
%           into segments. A gradient calculation is executed for each segment. 
%           Important for multi-segment data, e.g. time-series that contain 
%           different experiments. 
%
%  gradient_type : 'RBF' for the analytic gradient, 'PER' for the periodic gradient, or 'coarse' for the 
%           numeric gradient using a difference quotient. 
%
% Returns: 
%
%   Data :  A matrix 'nr_variables-by-nr_samples' of gradients, i.e. one 
%   gradient vector for the whole time-series (all segments) in each row.
%
% Note:
%
%   Output gradient is z-score transformed at the very end. Disable this if
%   you do not wish transformation.
%

function [All_gradients] = main_GCGP(ts_data, ts_data_changepoints, gradient_type)

    
    
    % load the GPstuff package
    currpath = pwd();
    cd GPstuff-4.4;
    startup;
    cd(currpath);
    
    % The gradients for each variable go into here
    All_gradients = [];
    
    
    %
    % Calculate analytic gradient with GPstuff 
    %   This calls main_GCGP_core() for the time-series of each variable.
    %  
    if(strcmp(gradient_type, 'RBF') || strcmp(gradient_type, 'PER'))
            
        last_cp = 1;
        
        if isempty(ts_data_changepoints)
            ts_data_changepoints = [size(ts_data, 2)];
        end
        
        for next_cp = ts_data_changepoints 
        
            % extract the next segment from the time-series 
            profile = ts_data(:,last_cp:next_cp); 

            % number of time-points
            nr_ts = size(profile, 2);
            
            %
            % These are the time-points, or 'x' variables for the Gaussian
            % Process. For Biopepa these time-points have a 2h gap. For all 
            % other time-series you could just replace it by
            %
            % ts = (1:size(profiles, 2))';
            %
            % Note, the gap should be invariant to the outcome of the GP because the
            % length parameter will be optimized anyway. However, the
            % initial length value can change the optimization path in
            % dependence of the gap. Here we have a gap of 2, and the
            % inital length will be 2, which leads to the best results in 
            % most cases. 
            
            ts = (0:2:(nr_ts*2-2));

            % Execute the gradient calculation for each variable
            tmpmat = [];
            for rowi = 1:size(profile,1)
                
                fprintf(' ** Gradient calculation for variable: %i, time-points: %i-%i\n', rowi, last_cp, next_cp);

                
                %
                % Execute the core gradient calculation script.
                % Note: This is a while loop because we try to execute it as long
                %       as we have warnings. This can happen with bad initial 
                %       parameters. See the comments in the catch clause below.
                %
                while 1
            
                    % set initial length, if the GP optimizer fails, adapt this
                    initial_GP_kernel_length = 2;
    
                    % run the GP
                    try

                        % attempt to run it -- should in most cases be successfull
                        result =  main_GCGP_core(profile(rowi,:), ts, gradient_type, initial_GP_kernel_length); 
                        

                        break; % if we are here, everything went fine.

                    catch %exception
                        %
                        % getGPgradient_main(), in particular the function fminscg can run
                        % into difficulties with optimization. It can hang with an endless
                        % loop with the warning message:
                        %
                        %      'Function value at xnew not finite or a number'
                        %
                        % This happens when there is a bad match of the length
                        % hyperparameter and the data. I inserted code that throws an error 
                        % whenever this message comes up 500 times, i.e. the warning loop 
                        % repeats 500 times.
                        %
                        % This error is caught here. By changing the initial length
                        % parameter of the GP, the warning messages can be avoided.
                        % Note: happens very rarely with Biopepa.
                        %
                        % Change the initial length -- this often helps with fminscg
                        initial_GP_kernel_length = initial_GP_kernel_length + 1;

                        % If we tried too many times, throw error
                        if initial_GP_kernel_length > 20
                           error('Could not adapt initial_length so that gpstuff optimization does not hang with warning: Function value at xnew not finite or a number'); 
                        end
                    end
                end            
                               
                tmpmat = [tmpmat; result];

            end

            % append gradient matrix for this segment to the main matrix
            All_gradients = [All_gradients, tmpmat];
            
            % remember last change-point
            last_cp = next_cp+1;
        end

    end  % end of gradient type - GP kernel
    
    
    
    %
    % Calculate numeric gradient with a difference quotient
    %
    if(strcmp(gradient_type, 'coarse'))   % i.e. coarse gradient

        % calculate the gradient, can be either coarse approximated (4 hours) or fine 
        All_gradients = zeros(size(ts_data));

        % first column
        All_gradients(:,1) = ts_data(:,2) - ts_data(:,1);

        % 2nd to second last column
        All_gradients(:,2:(end-1)) = (ts_data(:,3:end) - ts_data(:,1:(end-2))) / 2; 

        % last column
        All_gradients(:,end) = ts_data(:,end) - ts_data(:,end-1);

    end

 
    %
    % z-score transformation of gradients.
    %    ! Each row (variable) is transformed at ones, i.e. there is no
    %    separate tranformation for different time segments (as happened with the 
    %    gradient above), and there is no joint transformation for all
    %    variables.
    fprintf('z-score transforming gradient vectors!\n');
    for i_grad =1:size(All_gradients,1)       

        std_of_i  = std(All_gradients(i_grad,:));

        % avoid NaN's
        if(std_of_i == 0)
            std_of_i = 1;
        end

        % do the transformation
        All_gradients(i_grad,:) = (All_gradients(i_grad,:) - mean(All_gradients(i_grad,:))) / std_of_i;

    end
 
end




% Function: main_GCGP_core()
%
% Core method that does the gradient calculation for an input vector 'y'
% and 'x'. 'x' contains the time line points and 'y' the data. 
%
% Arguments:
%
%  y : input data vector of one variable
%  x : input x (time point) data vector. In most cases this can be just (1,2,..,nr_samples)
%  kernel : 'RBF' or 'PER'
%  inital_length : inital value for the length hyper-parameter.
%
% Returns:
%
%  pred_grad : A gradient vector with the size of the input vector y. 
%
function [pred_grad] = main_GCGP_core(y, x, kernel, initial_length)
  
    % The initial length scale is a hyper-parameter of the GP covariance function 
    % If it is not passed as argument, just set it to 1.
    if isempty(initial_length)
        initial_length = 1;
    end

    if initial_length < 0 
       error('initial_length in main_GCGP_core() too small, select larger or equal to 0.') 
    end

    % transpose input x vector
    x = x';
    
    % Set initial GP model hyper-parameters, the optimizer will adapt them 
    % to the data.
    initial_period_length = 6; % Only used for PER (periodic kernel)
    initial_magnSigma2 = 1;
    initial_noiseSigma2 = 0.1;
 
    % shift to a mean of zero so the GP can handle it
    y = (y - mean(y)); 
    y = y';

    % Hyperparameter function priors
    % Currently not used but could be used for tweaking.
    pl = prior_t();
    pm = prior_sqrtt();

    
    %
    % This is also the default optimization function in GPstuff. Might be 
    % overwritten with another method in case it is necessary for a
    % specific kernel. The Matern kernel can have difficulties with fminscg
    % for some data. 
    %
    optim_function = @fminscg;

    
    if strcmp(kernel, 'RBF')    

        % set 'RBF' kernel data structure
        gpcf = gpcf_sexp('lengthScale', initial_length, 'magnSigma2', initial_magnSigma2); % , ...
                    % 'lengthScale_prior', pl, 'magnSigma2_prior', pm);
                    % using the priors above will result in slightly improved biopepa results
   
    elseif strcmp(kernel, 'PER')  % periodic kernel

        % set 'PER' kernel data structure
        gpcf = gpcf_periodic('lengthScale', initial_length, 'period', initial_period_length, ...
            'magnSigma2', initial_magnSigma2); 

        % breaks it..
        % gpcf = gpcf_periodic(gpcf, 'lengthScale_prior', pl, 'magnSigma2_prior', pm);


    elseif strcmp(kernel, 'MAT32')  % Matern nu = 3/2 kernel

        initial_magnSigma2 = 1;
        initial_length = initial_length - 1;  % make it smaller than the input
        
        if initial_length <= 0
            initial_length = initial_length + 1;
        end
        
        gpcf = gpcf_matern32('lengthScale', initial_length, 'magnSigma2', initial_magnSigma2); 

      
    elseif strcmp(kernel, 'MAT52')  % Matern nu = 5/2 kernel

        initial_magnSigma2 = 1;
        initial_length = 1;
        gpcf = gpcf_matern52('lengthScale', initial_length, 'magnSigma2', initial_magnSigma2); 

    end


    % Define Gaussian likelihood with noise hyper-parameter 'sigma2' 
    lik = lik_gaussian('sigma2', initial_noiseSigma2);

    % Define the GP with likelihood and kernel function
    gp = gp_set('lik', lik, 'cf', gpcf, 'jitterSigma2', 1e-9);

    % set the numeric precision of the optimization 
    opt = optimset('TolFun',1e-3,'TolX',1e-3); 

    % Optimize with the scaled conjugate gradient method
    gp = gp_optim(gp, x, y,'opt', opt, 'optimf', optim_function);
    
    % get covariance matrix 'C', this is:   C = K + sigma2_noise * I
    [~, C] = gp_trcov(gp, x);

    % invert it
    inv_C = inv(C);

    % Extract GP hyperparameters
    length = gp.cf{1}.lengthScale;
    sigma2 = gp.cf{1}.magnSigma2;

    % if it is the periodic kernel, extract period hyper-parameter
    if strcmp(kernel, 'PER')
        period = gp.cf{1}.period;
    end

    % the gradients go into here
    pred_grad = [];

    %
    % Calculate the gradient for each point in 'y'.
    % Note, this could be done in vector notation, but I keep it simple to
    % make it clearer what happens. 
    %
    for j = 1:size(y,1)

        % The mean could be calculated with
        % pred_mean = [pred_mean, K(j,:) * inv_C * y];

            
        if strcmp(kernel, 'RBF')
            K_star = rbf_kernel_deriv(x(j), x, length, sigma2);
        elseif strcmp(kernel, 'PER')
            K_star = per_kernel_deriv(x(j), x, length, sigma2, period);
        %elseif strcmp(kernel, 'MAT32')
        %    K_star = rbf_kernel_deriv(x(j), x, length, sigma2);
        % elseif strcmp(kernel, 'MAT52')
        %    K_star = rbf_kernel_deriv(x(j), x, length, sigma2);
        end

        % calculates the gradient analytically every 2h
        pred_grad = [pred_grad, K_star * inv_C * y];

    end

   
end
    
    
% Function rbf_kernel_deriv()
%
% Derivative of the RBF kernel
%
function [K_star] = rbf_kernel_deriv(x_star, x, length, sigma2)

    K_star = zeros(1,size(x,1));

    for ind_i = 1:size(x,1)

        x2 = x(ind_i);
       
        % squared euclidian distance
        n2 = (x_star - x2)^2;   
        
        l2 = 1./length;
        l2 = l2^2;
        
        C = 1/2 * l2 * n2;
      
        % the derivative
        deriv_tmp = - (sigma2 * (2*x_star - 2*x2)) / (2*length^2 * exp(C));
        
        K_star(ind_i) = deriv_tmp; 
    
    end
end


% Function per_kernel_deriv()
%
% Derivative of the periodic kernel
%
function [K_star] = per_kernel_deriv(x_star, x, length, sigma2, period)

    K_star = zeros(1,size(x,1));

    for ind_i = 1:size(x,1)

        x2 = x(ind_i);
       
        % squared euclidian distance
        n2 = (x_star - x2)^2;   
        
        l2 = 1./length;
        l2 = l2^2;
        
        C = 1/2 * l2 * n2;
        
        % the derivative
        deriv_tmp = -(4*pi*sigma2*exp(-(2*sin((pi*(x_star - x2))/period)^2)/length^2)*cos((pi*(x_star - x2)) / period) * sin((pi*(x_star - x2)) / period))/(length^2 * period);
        
        K_star(ind_i) = deriv_tmp;
   
    end
end

   
    
    
