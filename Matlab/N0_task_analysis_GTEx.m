%----------------------------------------------------------------%
%   Running Cellfie on all the clean GTEx data and  saving    %
%  the resulting scores in the appropriate folder -                 %
%----------------------------------------------------------------%

initCobraToolbox(false);
p = parpool();
%% Loading 
model = load("../Models/Human-GEM_Calcium.mat");model = model.model;
model.lb(model.lb >0) = 0;
model.c(model.c ~= 0) = 0;


taskReport = readcell('../Results/0_Clean_Data/task_report.txt');
isEssential = readmatrix('../Results/0_Clean_Data/essential_logi.txt');
isUsed = readmatrix('../Results/0_Clean_Data/used_logi.txt');


taskStruct = generateTaskStructure_2022('..\Data\MetabolicTasks_CARDIO.xlsx');


%% Cellfie with the complete dataset
data.gene = model.genes;
vals = readtable("../Results/0_Clean_Data/geTMM_GTEx_noBC.csv");
vals = vals(ismember(vals.gene,model.genes),:);
%vals{:,2:end} = log10(vals{:,2:end} + 1);
data.value = nan(size(model.genes,1),size(vals,2)-1);
for i = 1:length(model.genes)
    if sum(strcmp(vals.gene,model.genes{i})) > 0
        data.value(i,:) = vals{strcmp(vals.gene,model.genes{i}),2:end};
    end
end

taskData.taskStructure = taskStruct;
taskData.essentialRxns = cell(size(taskStruct,1),1);
for i = 1:size(taskStruct,1)
    taskData.essentialRxns{i} = model.rxns(isEssential(:,i) == 1) ;
end

param.doParallel = 4;
param.getDetails = 0;

[score, score_binary ,~, ~]=CellFie_2022(data,taskData,model,param);
writematrix(score,'../Results/1_GTEx_clustering/Cellfie_scores/cellfie_score.txt');
writematrix(score_binary,'../Results/1_GTEx_clustering/Cellfie_scores/ cellfie_score_binary.txt');

