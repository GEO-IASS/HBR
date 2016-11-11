%
% Code by Marco  (2013)
%  Some minor changes by Andrej Aderhold (2014)
%
% Input: 
%      TrueVector   - vector of true network, i.e. gold standard
%      Aposterioris - vector of inferred network from method
%
 
function [auroc, auprec] = AUROC_AUPREC_scores(TrueVector, Aposterioris)


    n_non_edges   = length(find(TrueVector == 0));
    n_edges       = length(find(TrueVector == 1));
   
    % sort and delete duplicates
    APO_values = unique(Aposterioris);
    
    % reverse
    APO_values = APO_values(end:-1:1);
    
    ROC_x = 0;
    ROC_y = 0;
    
    TPS = 0;
    FPS = 0;

    
    for i = 1:length(APO_values)

        % find elements in Aposterioris that is above theshold
        pos_indicis = find(Aposterioris >= APO_values(i));
        % and these are the ones below
        neg_indicis = setdiff(1:length(Aposterioris), pos_indicis);

        % extract number of edges that are also true and false 
        TP = sum(TrueVector(pos_indicis) == 1);
        FP = sum(TrueVector(pos_indicis) == 0);
        TN = sum(TrueVector(neg_indicis) == 0);
        
        % original code:
        %    TP = length(find(MATRIX==1 & True_Matrix==1));
        %    TN = length(find(MATRIX==0 & True_Matrix==0)); % minus diagonal elements
        %    FP = length(find(MATRIX==1 & True_Matrix==0));
        
        
        % save this for later AUPREC calculation 
        TPS = [TPS,TP];
        FPS = [FPS,FP];

        % AUROC related
        Sensitivity = TP/n_edges;
        inv_Specif  = 1 - (TN/n_non_edges);

        ROC_y = [ROC_y, Sensitivity];
        ROC_x = [ROC_x, inv_Specif];

    end

    ROC_x(end+1) = 1;
    ROC_y(end+1) = 1;

    auroc = trapz(ROC_x,ROC_y);

    %plot(ROC_x, ROC_y);

    %
    % start AUPREC calc.
    %

    


    for i=2:length(TPS)
        if (TPS(i)-TPS(i-1))>1

            NEW_TPS = [];
            NEW_FPS = [];

            for x = 1:(TPS(i)-TPS(i-1)-1)
                skew    = (FPS(i)-FPS(i-1))/(TPS(i)-TPS(i-1));
                NEW_TPS = [NEW_TPS,TPS(i-1)+x];
                NEW_FPS = [NEW_FPS,FPS(i-1)+ skew*x];
            end

            TPS = [TPS(1:i-1),NEW_TPS,TPS(i:end)];
            FPS = [FPS(1:i-1),NEW_FPS,FPS(i:end)];

        end

    end


    PRECISION = TPS(2:end)./(TPS(2:end)+FPS(2:end));
    RECALL    = TPS(2:end)/n_edges;

    auprec = trapz(RECALL,PRECISION);



return;

