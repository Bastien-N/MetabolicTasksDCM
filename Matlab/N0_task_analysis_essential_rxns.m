%--------------------------------------------------------------------------------%
%   Getting the essential reactions for the model and the current taskList   %
%-------------------------------------------------------------------------------%
initCobraToolbox(false);
model = load("../Models/Human-GEM_Calcium.mat");model = model.model;
model.lb(model.lb >0) = 0;
model.c(model.c ~= 0) = 0;

taskStruct = generateTaskStructure_2022('..\Data\MetabolicTasks_CARDIO.xlsx');
param.taskStructure = taskStruct;
param.printOutput = 1;
param.getEssential = 1;
param.saveUsed = 1;
[taskReport,essentialRxns,~,usedRxns] =  checkMetabolicTasksHumanGEM_2023_2(model,param);

isEssential = zeros(size(model.lb,1),size(taskStruct,1));
isUsed = isEssential;

for i = 1:length(taskStruct)
    if ~isempty(essentialRxns{i})
        isEssential(:,i) = ismember(model.rxns,essentialRxns{i});
    end
    if ~isempty(usedRxns{i})
        isUsed(:,i) = ismember(model.rxns,usedRxns{i});
    end
end

writecell(taskReport,'..\Results\0_Clean_Data\task_report.txt');
writematrix(isEssential,'..\Results\0_Clean_Data\essential_logi.txt');
writematrix(isUsed,'..\Results\0_Clean_Data\used_logi.txt');

