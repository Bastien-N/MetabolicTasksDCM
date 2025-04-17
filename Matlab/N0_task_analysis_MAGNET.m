
%-------------------------------------------------------------%
%   Running Cellfie on all the clean MAGNET data and      %
% saving the resulting scores in the appropriate folder -   %
%-------------------------------------------------------------%
%Initiating the cobra toolbox
initCobraToolbox(false);
%Starting the parallel processing toolbox
p = parpool(4);
%% Loading 
model = load("../Models/Human-GEM_Calcium.mat");model = model.model; %Loading the model
model.lb(model.lb >0) = 0;%This makes sure that the model does not have an objective nor forced flux (lower bound > 0)
model.c(model.c ~= 0) = 0;

taskReport = readcell('../Results/0_Clean_Data/task_report.txt');
isEssential = readmatrix('../Results/0_Clean_Data/essential_logi.txt'); %Reading the essential reactions as a logical array
isUsed = readmatrix('../Results/0_Clean_Data/used_logi.txt');

taskStruct = generateTaskStructure_2022('..\Data\MetabolicTasks_CARDIO.xlsx'); %Reading the task list
taskStruct([277 258]) = [];
isEssential(:,[277 258]) = [];
%% Cellfie with the complete dataset
%Creating the data structure
data.gene = model.genes;
vals = readtable("../Results/0_Clean_Data/geTMM_MAGNET_BC.csv");
vals = vals(ismember(vals.gene,model.genes),:);%We only want the genes that are in the model
%vals{:,2:end} = log10(vals{:,2:end} + 1);
data.value = nan(size(model.genes,1),size(vals,2)-1);
for i = 1:length(model.genes)
    if sum(strcmp(vals.gene,model.genes{i})) > 0
        data.value(i,:) = vals{strcmp(vals.gene,model.genes{i}),2:end};
    end
end
%data.value(data.value < 1) = 1;
taskData.taskStructure = taskStruct;
taskData.essentialRxns = cell(size(taskStruct,1),1);
for i = 1:size(taskStruct,1)
    taskData.essentialRxns{i} = model.rxns(isEssential(:,i) == 1) ;%This loops is where we get rid of genes not in the model
end

param.doParallel = 4;%Not necessary if you skipped parpool()
param.getDetails = 0;% Absolutely Required!
% [mappedData,genes_used] = map_exprData_to_GPR(model,vals);
% writetable(mappedData,'../Results/3_MAGNET_NFvDCM/Cellfie_scores/mapped_expressionData.txt')
% writecell(genes_used,'../Results/3_MAGNET_NFvDCM/Cellfie_scores/genes_used.txt');

[score, score_binary ,~, ~,expression]=CellFie_2022(data,taskData,model,param);
writematrix(expression.Rxns,'../Results/3_MAGNET_NFvDCM/Cellfie_scores/reaction_activity_score.txt');
for i = 1:size(expression.gene_used,1)
    for ii = 1:size(expression.gene_used,2)
        if ~isempty(expression.gene_used{i,ii})
            expression.gene_used{i,ii} = expression.gene_used{i,ii}{1}; 
        end
    end
end
writecell(expression.gene_used,'../Results/3_MAGNET_NFvDCM/Cellfie_scores/genes_used.txt');
writematrix(expression.count,'../Results/3_MAGNET_NFvDCM/Cellfie_scores/gene_used_promiscuity.txt');

writematrix(score,'../Results/3_MAGNET_NFvDCM/Cellfie_scores/cellfie_score.txt');
writematrix(score_binary,'../Results/3_MAGNET_NFvDCM/Cellfie_scores/cellfie_score_binary.txt');

%% Cellfie with the DCM dataset
data.gene = model.genes;
vals = readtable("../Results/0_Clean_Data/geTMM_MAGNET_BC_DCM.csv");
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

[score, score_binary ,~, ~,expression]=CellFie_2022(data,taskData,model,param);
writematrix(expression.Rxns,'../Results/4_MAGNET_DCM_clustering/Cellfie_scores/reaction_activity_score.txt');
for i = 1:size(expression.gene_used,1)
    for ii = 1:size(expression.gene_used,2)
        if ~isempty(expression.gene_used{i,ii})
            expression.gene_used{i,ii} = expression.gene_used{i,ii}{1}; 
        end
    end
end
writecell(expression.gene_used,'../Results/4_MAGNET_DCM_clustering/Cellfie_scores/genes_used.txt');
writematrix(expression.count,'../Results/4_MAGNET_DCM_clustering/Cellfie_scores/gene_used_promiscuity.txt');

writematrix(score,'../Results/4_MAGNET_DCM_clustering/Cellfie_scores/cellfie_score.txt');
writematrix(score_binary,'../Results/4_MAGNET_DCM_clustering/Cellfie_scores/cellfie_score_binary.txt');

%% Cellfie with the NF dataset
data.gene = model.genes;
vals = readtable("../Results/0_Clean_Data/geTMM_MAGNET_BC_NF.csv");
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

[score, score_binary ,~, ~,expression]=CellFie_2022(data,taskData,model,param);
writematrix(expression.Rxns,'../Results/2_MAGNET_NF_MvF/Cellfie_scores/reaction_activity_score.txt');
for i = 1:size(expression.gene_used,1)
    for ii = 1:size(expression.gene_used,2)
        if ~isempty(expression.gene_used{i,ii})
            expression.gene_used{i,ii} = expression.gene_used{i,ii}{1}; 
        end
    end
end
writecell(expression.gene_used,'../Results/2_MAGNET_NF_MvF/Cellfie_scores/genes_used.txt');
writematrix(expression.count,'../Results/2_MAGNET_NF_MvF/Cellfie_scores/gene_used_promiscuity.txt');

writematrix(score,'../Results/2_MAGNET_NF_MvF/Cellfie_scores/cellfie_score.txt');
writematrix(score_binary,'../Results/2_MAGNET_NF_MvF/Cellfie_scores/cellfie_score_binary.txt');

