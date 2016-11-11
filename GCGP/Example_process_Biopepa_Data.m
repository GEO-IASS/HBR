
% This is the function in which the file names from the Biopepa txt files 
% are specified and the folder they are located in. 
% Everything else is done in function Protocol() which generates a single 
% data structure containing all the covariable matrices with appended
% response row and each of the indendent data replicas (see description of
% the return type of function Protocol() below)


function Example_process_Biopepa_Data()


   
    
    % Index of mRNAs in the Biopepa data time-series
    mRNA_response_set = [2,4,8,10,12,14,17]; 


    % Specify the gradient response to generate
    % Can be  
    %    'RBF'    :  analytic gradient using a Gaussian Process with RBF 
    %                kernel, the GPstuff package is used for this.
    %         
    % or 'coarse' :  numeric gradient using difference quotient 
    %
    gradient_type = 'RBF';
  
    % From where to read the raw data. This data comes from Biopepa runs 
    % using the Biopepa Eclipse Plug-in and the Biopepa script
    % 'arabidopsis_clock.biopepa'. 
    folder_raw_data = 'Raw_Biopepa_Data';
      
    % This is the header for the gradient .csv file
    head_Gradient = {'"GI_mRNA"','"LHY_mRNA"','"PRR5_NI_mRNA"','"PRR7_mRNA"','"PRR9_mRNA"','"TOC1_mRNA"','"Y_mRNA"'};
    head_Gradient = sprintf('%s,',head_Gradient{:});
    head_Gradient(end)='';
    
    % This is the header for the time-series concentration .csv file
    head_allVars = {'"GI:ZTL"', '"GI_mRNA"', '"GI_prot"', '"LHY_mRNA"', '"LHY_prot"', '"LHY_prot_modif"', '"Light*P"', '"PRR5_mRNA"', '"PRR5_prot"', '"PRR7_mRNA"', '"PRR7_prot"', '"PRR9_mRNA"', '"PRR9_prot"', '"TOC1_mRNA"', '"TOC1_prot"', '"TOC1_prot_mod"', '"Y_mRNA"', '"Y_prot"'};
    head_allVars = sprintf('%s,',head_allVars{:});
    head_allVars(end)='';
    
    % There are 20 variables in the Biopepa data set but at the end of this
    % script we only save the concentrations of the first 18, corresonding
    % to the 18 variables in 'head_allVars'
    relevant_predictor_idx = 1:18;

    % We have six network variants, one is the wildtype that corresponds
    % to the P2010 network, and five pruned networks that contain different knock-outs
    %
    % The coding is:
    % type = 0 -> no knock out (wildtype)
    %        1 -> PRR7/PRR9 knock out
    %        2 -> PRR7/PRR9/TOC1 knock-out
    %        3 -> TOC1 knock-out
    %        4 -> PRR5/PRR7/PRR9 knock-out
    %        5 -> PRR5/PRR7/PRR9/TOC1 knock-out

    for knockout_type=0:5


        file_list = {};
        file_folder = '';

        if knockout_type == 0   % without any constant knockout

            file_list = {'biopepa_knockout_GI' 
                         'biopepa_knockout_LHY'
                         'biopepa_knockout_TOC1'
                         'biopepa_knockout_PRR7_PRR9'
                         'biopepa_wt_DD'
                         'biopepa_wt_LD'
                         'biopepa_wt_LL'
                         'biopepa_wt_photoperiod_4h'
                         'biopepa_wt_photoperiod_6h'
                         'biopepa_wt_photoperiod_8h'
                         'biopepa_wt_photoperiod_18h'
                     };

            % set the folder where the files reside
            file_folder = 'wildtype/';
            file_name_out = 'biopepa_wildtype';

        elseif knockout_type == 1   % The constant knock out of PRR7,PRR9 proteins


            file_list = {'biopepa_knockout_GI_PRR7_PRR9'
                         'biopepa_knockout_LHY_PRR7_PRR9'   
                         'biopepa_knockout_TOC1_PRR7_PRR9'
                         'biopepa_knockout_PRR7_PRR9'
                         'biopepa_wt_DD_knockout_PRR7-PRR9'
                         'biopepa_wt_LD_knockout_PRR7_PRR9'
                         'biopepa_wt_LL_knockout_PRR7-PRR9'
                         'biopepa_wt_photoperiod_4h_knockout_PRR7-PRR9'
                         'biopepa_wt_photoperiod_6h_knockout_PRR7-PRR9'
                         'biopepa_wt_photoperiod_8h_knockout_PRR7-PRR9'
                         'biopepa_wt_photoperiod_18h_knockout_PRR7-PRR9'
                     };

            file_folder = 'mutant_PRR7_PRR9/';       
            file_name_out = 'biopepa_mutant_PRR7_PRR9';

        elseif knockout_type == 2  % PRR7, PRR9, TOC1

            file_list = {'biopepa_knockout_GI_PRR7_PRR9_TOC1'
                         'biopepa_knockout_LHY_PRR7_PRR9_TOC1'
                         'biopepa_knockout_PRR7_PRR9_TOC1'
                         'biopepa_knockout_PRR7_PRR9_TOC1_dupl1'
                         'biopepa_DD_knockout_PRR7_PRR9_TOC1'
                         'biopepa_LD_knockout_PRR7_PRR9_TOC1'
                         'biopepa_LL_knockout_PRR7_PRR9_TOC1'
                         'biopepa_4h_knockout_PRR7_PRR9_TOC1'
                         'biopepa_6h_knockout_PRR7_PRR9_TOC1'
                         'biopepa_8h_knockout_PRR7_PRR9_TOC1'
                         'biopepa_18h_knockout_PRR7_PRR9_TOC1'
                     };

            file_folder = 'mutant_PRR7_PRR9_TOC1/';     
            file_name_out = 'biopepa_mutant_PRR7_PRR9_TOC1';

        elseif knockout_type == 3   % TOC1

            file_list = {'biopepa_knockout_GI_TOC1'
                         'biopepa_knockout_LHY_TOC1'
                         'biopepa_knockout_TOC1'
                         'biopepa_knockout_PRR7_PRR9_TOC1'
                         'biopepa_DD_knockout_TOC1'
                         'biopepa_LD_knockout_TOC1'
                         'biopepa_LL_knockout_TOC1'
                         'biopepa_4h_knockout_TOC1'
                         'biopepa_6h_knockout_TOC1'
                         'biopepa_8h_knockout_TOC1'
                         'biopepa_18h_knockout_TOC1'
                     };

            file_folder = 'mutant_TOC1/';     
            file_name_out = 'biopepa_mutant_TOC1';

        elseif knockout_type == 4   % PRR5, PRR7, PRR9

            file_list = {'biopepa_knockout_GI_PRR5_PRR7_PRR9'
                         'biopepa_knockout_LHY_PRR5_PRR7_PRR9'
                         'biopepa_knockout_TOC1_PRR5_PRR7_PRR9'
                         'biopepa_knockout_PRR5_PRR7_PRR9'
                         'biopepa_DD_knockout_PRR5-7-9'
                         'biopepa_LD_knockout_PRR5-7-9'
                         'biopepa_LL_knockout_PRR5-7-9'
                         'biopepa_4h_knockout_PRR5-7-9'
                         'biopepa_6h_knockout_PRR5-7-9'
                         'biopepa_8h_knockout_PRR5-7-9'
                         'biopepa_18h_knockout_PRR5-7-9'
                     };

            file_folder = 'mutant_PRR5_PRR7_PRR9/';     
            file_name_out = 'biopepa_mutant_PRR5_PRR7_PRR9';

        elseif knockout_type == 5

            file_list = {'biopepa_knockout_GI_PRR5_PRR7_PRR9_TOC1'
                         'biopepa_knockout_LHY_PRR5_PRR7_PRR9_TOC1'
                         'biopepa_knockout_PRR5_PRR7_PRR9_TOC1'
                         'biopepa_knockout_PRR5_PRR7_PRR9_TOC1_dupl1'
                         'biopepa_DD_knockout_PRR5-7-9_TOC1'
                         'biopepa_LD_knockout_PRR5-7-9_TOC1'
                         'biopepa_LL_knockout_PRR5-7-9_TOC1'
                         'biopepa_4h_knockout_PRR5-7-9_TOC1'
                         'biopepa_6h_knockout_PRR5-7-9_TOC1'
                         'biopepa_8h_knockout_PRR5-7-9_TOC1'
                         'biopepa_18h_knockout_PRR5-7-9_TOC1'
                     };

            file_folder = 'mutant_PRR5_PRR7_PRR9_TOC1/';     
            file_name_out = 'biopepa_mutant_PRR5_PRR7_PRR9_TOC1';

        end

      
        % Repeat for the five independent data instances of each network 
        for dataid = 1:5

            % Run the main routine for the extraction of the relevant
            % observations, and the calculation of the gradient
                       
         
            fprintf('\n ** gradient_type: %s, knockout_type: %i, dataid: %i \n\n', gradient_type, knockout_type, dataid);

            % Read the raw data into this structure
            Raw_data = {};

            % Read in all the files -- each file contains a specific experiment. 
            % This can be a light experiment or a knockout
            for n_data = 1:length(file_list)

                filetoread = sprintf('%s/%s/%s_id%i.txt', folder_raw_data, file_folder, file_list{n_data}, dataid );
                filecontent = importdata(filetoread); 

                % read the data, but only from observation 1200 to 1440. Everything
                % before 1200 is the entrainment of the clock and ignored. This is 
                % Biopepa specific.
                data_tmp   = filecontent.data(1200:1440,2:end);
                data_tmp   = data_tmp';

                % Set 'Light*P' variable. This variable is important and remains 
                % in the output at index 7. The original element 7 is 'P' and
                % element 20 is 'Light'
                data_new(7,:) = data_tmp(7,:) .* data_tmp(20,:);
                
                % Put each time-series in a single cell, they are later combined to the full
                % time-series, but before that, they are thinned out and the
                % gradient is calculated. 
                Raw_data{n_data}  = data_tmp;

            end

            %
            % The extracted samples from the raw time-series are saved in
            % here
            %
            Data_concentrations = [];
            change_points = [];
            
            %
            % Loop over each component/experiment that goes into the final
            % time-series. For each component we calculate the gradient. It is not
            % advisable to calculate the gradient over the whole time-series if it
            % is segmented into different small time-series.
            %
            for index2 = 1:length(Raw_data) 

                % get the time-series for the current component/experiment
                ts_data = Raw_data{index2};

                % Thinning: Only keep every 20th observation.
                points_to_take = 1:20:size(ts_data, 2);

                % Extract these observations      
                ts_data  = ts_data(:,points_to_take); 


                % The following is Biopepa data specific handling, but we have seen
                % that it is beneficial to replace Zero concentration data from knock-outs
                % with random data.
                %
                %  If there is a knock-out of a specific gene, replace both, the 
                %  mRNA and protein concentrations with a random number because they are 0
                %  
                if index2 == 1         % -> GI knockout experiment

                    ts_data(2,:) = randn(1, size(ts_data,2));  % GI_mRNA
                    ts_data(3,:) = randn(1, size(ts_data,2));  % GI_prot

                elseif index2 == 2     % -> LHY knockout experiment

                    ts_data(4,:) = randn(1, size(ts_data,2));  % LHY_mRNA
                    ts_data(5,:) = randn(1, size(ts_data,2));  % LHY_prot

                elseif index2 == 3     % -> knock out of TOC1 

                    ts_data(14,:) = randn(1, size(ts_data,2));  % TOC1_mRNA
                    ts_data(15,:) = randn(1, size(ts_data,2));  % TOC1_prot
                    ts_data(16,:) = randn(1, size(ts_data,2));  % TOC1_modif

                elseif index2==4      % -> knock out of PRR7/PRR9 

                    ts_data(10,:) = randn(1, size(ts_data,2));  % PRR7_mRNA
                    ts_data(12,:) = randn(1, size(ts_data,2));  % PRR9_mRNA

                    ts_data(11,:) = randn(1, size(ts_data,2));  % PRR7_prot
                    ts_data(13,:) = randn(1, size(ts_data,2));  % PRR9_prot

                end

                if(knockout_type == 1)  % constant knock out of PRR7 and PRR9

                    ts_data(11,:) = randn(1, size(ts_data,2));  % PRR7
                    ts_data(13,:) = randn(1, size(ts_data,2));  % PRR9

                elseif(knockout_type == 2)  % constant knock out of PRR7, PRR9 and TOC1

                    ts_data(11,:) = randn(1, size(ts_data,2));  % PRR7
                    ts_data(13,:) = randn(1, size(ts_data,2));  % PRR9
                    ts_data(15,:) = randn(1, size(ts_data,2));  % TOC1
                    ts_data(16,:) = randn(1, size(ts_data,2));  % TOC1_modif

                elseif (knockout_type == 3)  % TOC1

                    ts_data(15,:) = randn(1, size(ts_data,2));  % TOC1
                    ts_data(16,:) = randn(1, size(ts_data,2));  % TOC1_modif

                elseif (knockout_type == 4)  % constant knock out of PRR5-7-9

                    ts_data(9,:) = randn(1, size(ts_data,2));  % PRR5
                    ts_data(11,:) = randn(1, size(ts_data,2));  % PRR7
                    ts_data(13,:) = randn(1, size(ts_data,2));  % PRR9

                elseif (knockout_type == 5)  % constant knock out of PRR5-7-9 and TOC1

                    ts_data(9,:) = randn(1, size(ts_data,2));  % PRR5
                    ts_data(11,:) = randn(1, size(ts_data,2));  % PRR7
                    ts_data(13,:) = randn(1, size(ts_data,2));  % PRR9
                    ts_data(15,:) = randn(1, size(ts_data,2));  % TOC1
                    ts_data(16,:) = randn(1, size(ts_data,2));  % TOC1_modif

                end

                % append data to covariable table
                Data_concentrations = [Data_concentrations, ts_data];

                % Save the total size of observations after each segment is
                % added. This change_point vector is used to split the
                % complete time-series into segments for which a gradient
                % is calculated separately. 
                change_points = [change_points, size(Data_concentrations, 2)];
                
                
            end  % end loop over different segment (experiments) in the time-series

            
            %
            % Calculate the gradients for all variables in the time-series
            %
            gradient = main_GCGP(Data_concentrations(mRNA_response_set,:), change_points, gradient_type);

            gradient = gradient';
            
            % Save the response gradients to this file
            outfile = sprintf('Processed_Data_CSV/%s_%sGradient_id%i.csv', file_name_out, gradient_type, dataid);

            % Write with header
            dlmwrite(outfile, head_Gradient, '');
            dlmwrite(outfile, gradient, '-append')

            
            %
            % z-score transformation of time-series concentration data
            %
            for i_node=(1:size(Data_concentrations,1)-1)  % exclude last one (light)

                std_of_i  = std(Data_concentrations(i_node,:));
                mean_of_i = mean(Data_concentrations(i_node,:));

                % if stdev is 0, replace by 1 to avoid division by zero.
                if(std_of_i == 0)
                    std_of_i = 1;
                end

                Data_concentrations(i_node,:) = (Data_concentrations(i_node,:) - mean_of_i) / std_of_i;			

                % set lowest value to 0
                Data_concentrations(i_node,:) = Data_concentrations(i_node,:) - min(Data_concentrations(i_node,:));		
            end
            
            % extract and save only the relevant elements
            Data_concentrations = Data_concentrations(relevant_predictor_idx,:)';
            
            % Save the time-series concentrations to this file
            outfile = sprintf('Processed_Data_CSV/%s_ts-concentrations_id%i.csv', file_name_out, dataid);

            % Write with header
            dlmwrite(outfile, head_allVars, '');
            dlmwrite(outfile, Data_concentrations, '-append')

        end
    end

end




%
