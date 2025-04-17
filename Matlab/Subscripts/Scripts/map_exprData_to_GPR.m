function [mappedData] = map_exprData_to_GPR(model,exprData)
%MAP_EXPRDATA_TO_GPR maps expression data to gpr for a whole expression
%dataset
%   USAGE:   [mappedData] = map_exprData_to_GPR(model,exprData)
%   INPUTS   model, the (Cobra) model
%                  exprData, the gene expression data as a matlab table
%                  (type of struct) in which the first column is called
%                  'gene' (!Important) and contains the genes IDs in the
%                  same format as in model.genes
%   OUTPUTS mappedDat, the mapped expression data as a table.
gprRules = GPRparser(model);
exprData.Properties.VariableNames{1} = 'gene';

mappedData = array2table(zeros(length(gprRules),size(exprData,2)));
mappedData.Var1 = model.rxns;
mappedData.Properties.VariableNames = exprData.Properties.VariableNames;
mappedData.Properties.VariableNames{1} = 'Reaction';

for i = 1:(size(exprData,2)-1)
    temp = exprData(:,[1 i+1]);
    temp.Properties.VariableNames = {'gene','value'};
    [mappedTemp,~] = mapExpressionToReactionsKnownGPR(model,temp,gprRules);
    mappedData(:,i+1) = num2cell(mappedTemp);
end
end

