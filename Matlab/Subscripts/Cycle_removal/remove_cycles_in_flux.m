function [cfFlux] = remove_cycles_in_flux(model,sol,er,restrictObj,par)
%An implementation of CycleFreeFlux that is closer to the original R
%implementation
%   INPUT
%               model: a model
%               sol: either a single vector of flux values, e.g. from an
%                     FBA, or a matrix of flux values  (rxn x N)
%               er: optional, a previously determined logical vector of
%                    exchanges reactions, lowers computationnal demand
%               restrictObj: optional, useless if er is provided, logical
%               indicating if the objective function should be included in
%               er.
%               par: optional (defaults to false) logical indicating if parallel computing is allowed.
%               if a parpool already exists it will be used and left as is.
%               Otherwise a pool will be created if N > 50 and destroyed
%               afterwards.
%
%   OUTPUT
%               cfFlux: a matrix of the same dimensions as sol, with
%               cycling fluxes removed
if ~exist('restrictObj','var')
    restrictObj = false;
end
if ~exist('er','var') || isempty(er)
     % third parameter does not exist, so default it to something
     if restrictObj
      [er,~] = findExcRxns(model,1); 
     else
       [er,~] = findExcRxns(model,0);  
     end
     
end

if ~exist('sol','var')
    sol = optimizeCbModel(model);
    sol = sol.v;
end
if nargin < 5
    par = false;
end
if size(sol,1) ~= 1 && size(sol,2) ~= 1 
    cfFlux = zeros(size(sol));
    sol0 = sol;
    if par && size(gcp('nocreate'),2) ~=0
        nestedWork = @Work;
        env = getEnvironment;
        parfor i = 1:size(sol,2)
            restoreEnvironment(env);
            sol = sol0(:,i);
            cfFlux(:,i) = nestedWork(sol,model,er);
        end
    elseif par && size(gcp('nocreate'),2) ==0 && size(sol0,2) > 50
        p = parpool;
        env = getEnvironment;
        nestedWork = @Work;
        parfor i = 1:size(sol,2)
            restoreEnvironment(env);
            sol = sol0(:,i);
            cfFlux(:,i) = nestedWork(sol,model,er);
        end
        delete(p);

    else
        for i = 1:size(sol,2)
            sol = sol0(:,i);
            cfFlux(:,i) = Work(sol,model,er);
        end
    end
else
    cfFlux = Work(sol,model,er);
end

%function used 
function restrictedFlux = Work(soltmp,modelTmp,erTmp)
[resrModel,solSigns] = restrict_to_sol(modelTmp,soltmp,erTmp);
resrModel.c = solSigns;
minOpt = optimizeCbModel(resrModel,'min');
padding = 0;
maxPadding = 1;
if minOpt.stat == 0 
    while minOpt.stat == 0 && padding <= maxPadding
        padding = padding + 0.1;
        [resrModel,solSigns] = restrict_to_sol(modelTmp,soltmp,erTmp,padding);
        resrModel.c = solSigns;
        minOpt = optimizeCbModel(resrModel,'min');
    end
end
restrictedFlux = minOpt.v;
end
end