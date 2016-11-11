%
% This script adds up and normalizes the DAGs from an MCMC Run.
% 
% It return a matrix with the edge probabilities. Entries that should be 
%  ignored are marked with NaN. These are fixed edges (flag '2') and
%  invalid/ignore edges (flag '-1'). 
%
% These NaN entries can later be removed from further score calculations.
%

function [Edge_matrix] = DAGs_to_edge_probabilities(DAGs)

    % get second DAG, this is the first one after the burn-in
    all_edges = DAGs{1};
            
    % Convert edges marked with '-1' and '2' to NaN -- these are
    % the "ignore-edges" and self-loops and they will be removed 
    % from the score calculation.

    all_edges(all_edges == 2) = NaN;
    all_edges(all_edges == -1) = NaN;

    % add up all the rest
    for sample_idx=2:length(DAGs)

        tmp = DAGs{sample_idx};

        % again mark self-loops and ignore-edges with NaN so they 
        % can be removed further below
        tmp(tmp == 2) = NaN;
        tmp(tmp == -1) = NaN;

        all_edges = all_edges + tmp;
    end

    % normalize
    Edge_matrix = all_edges / (length(DAGs) - 1);
            
end