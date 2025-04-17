function [expr,gene_used] = mapExpressionToReactionsKnownGPR(model,expressionData,GPR)
% Determines the expression data associated to each reaction present in
% the model from known GPR rules. Useful when you need to map a lot of data
% to the same model.
%
% USAGE:
%    expr = mapExpressionToReactionsKnownGPR(model,expressionData,GPR) 
%
% INPUTS:
%	model                   model structure
%	expressionData          mRNA expression data structure
%       .gene               	cell array containing GeneIDs in the same
%                               format as model.genes
%       .value                  Vector containing corresponding expression
%                               value (FPKM/RPKM)
%   GPR                     parsed GPR rules as in the output of
%                           mapExpressionToReactions
% 
% OUTPUTS:
%   expr:         reaction expression, corresponding to model.rxns.
%
% .. Author:
%        Bastien Nihant, using code extracted from Anne Richelle's
%        function: mapExpressionToReactions (Anne Richelle, May 2017 -
%        integration of new extraction methods)
 
% Find wich genes in expression data are used in the model
[gene_id, gene_expr] = findUsedGenesLevels(model,expressionData);
% Link the gene to the model reactions
[expr, gene_used] = selectGeneFromGPR(model, gene_id, gene_expr, GPR, 0);
end

