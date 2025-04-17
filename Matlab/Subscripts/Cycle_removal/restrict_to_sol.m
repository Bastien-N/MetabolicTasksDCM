function [restrictedModel,signs] = restrict_to_sol(model,sol,er,relax)
%RESTRICT_TO_SOL Summary of this function goes here
%   Detailed explanation goes here
restrictedModel = model;
if ~exist('relax','var')
    relax = 0;
end
for ii = 1:size(sol,1)            
        if sol(ii) >= 0
            restrictedModel.ub(ii) = sol(ii)+relax;
            restrictedModel.lb(ii) = max(0,restrictedModel.lb(ii)-relax);
        elseif sol(ii) < 0
            restrictedModel.lb(ii) = sol(ii)-relax;
            restrictedModel.ub(ii) = min(0,restrictedModel.ub(ii)+relax);
        end
end
restrictedModel.ub(er) = sol(er);
restrictedModel.lb(er) = sol(er);

signs = sign(sol);
end
