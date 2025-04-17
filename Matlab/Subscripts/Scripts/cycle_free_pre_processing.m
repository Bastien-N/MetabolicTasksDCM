function model = cycle_free_pre_processing(model0,allowCycles)
% A modification of the preprocessing function developed by Braunstein et
% al for the paper:
% 
%   "Braunstein, A. et al., An analytic approximation of the feasible space
%       of metabolic networks. 
%   Nat Commun 8, 14915 (2017). https://doi.org/10.1038/ncomms14915"
%
% The original code for that study from above paper which can also be accessed at:
% https://github.com/anna-pa-m/Metabolic-EP
% 
% The modification is to use cycle-free FVA to find the feasible reaction
% bounds
% USAGE:
%   model = pre_processing(model0);
%
% INPUT
%   model0:   a genome-scale metabolic model that is loaded into MATLAB
%   allowCycles: Whether loops are allowed in solution or which method to block loops.
% 
%                       * 1 (or true) : loops allowed (default)
%                       * 0 (or false): loops not allowed. Use  to find loopless solutions
%                       * 'original'  : original loopless FVA
%                       * 'fastSNP'   : loopless FVA with with Fast-SNP preprocessing of nullspace
%                       * 'LLC-NS'    : localized loopless FVA using information from nullsapce
%                       * 'LLC-EFM'   : localized loopless FVA using information from EFMs.
%                       Require CalculateFluxModes.m from EFMtool to calculate EFMs.
%
% OUTPUT
%   model:    preprocessed model
%
% ..Author: 
%   Bastien Nihant, modifying the work of:
%   Alfredo Braunstein, Andrea Pagnani and % Anna Paola Muntoni (April 2017)
%   Chaitra Sarathy, 31 Aug 2020 (added documentation)

if ~exist("allowCycles","var")
    allowCycles = 0;
end
 %#function fmincon
    model = [];
%    Nfluxes = length(model0.lb);
    fprintf('Preprocessing..\n');

    if sum(model0.c ~= 0) >0
        model0.c(model0.c ~= 0) = 0;
        disp('Removed the objective')
    end
    [lb,ub] = fluxVariability(model0,0,[],[],[],allowCycles);
    fixed = (lb == ub);
    fprintf('%d variables are fixed\n', nnz(fixed));
    notfixed = (lb < ub);
    model.S = model0.S(:,notfixed);
    model.c = model0.c(notfixed);
    usedMets = sum(model.S ~= 0,2) > 0;
    model.S = model.S(usedMets,:);
    model.b = model0.b(usedMets) - model0.S(usedMets,fixed) * lb(fixed);
    model.lb = lb(notfixed);
    model.ub = ub(notfixed);
    model.rxns = model0.rxns(notfixed);
    model.rxnNames = model0.rxnNames(notfixed);
    model.mets = model0.mets(usedMets);
    model.metNames = model0.metNames(usedMets);
    model.metFormulas = model0.metFormulas(usedMets);
    model.genes = model0.genes;
    model.rules = model0.rules(notfixed);
    str = sprintf('%s pre-processed', model0.description);
    model.description = str; 
    model = removeUnusedGenes(model);
end