function [coreMat] = localT2_core(mappedDat,lowerThres,upperThres)
%LOCALT2_CORE Uses localT2 to define the core set of reactions
%   INPUTS:     
%       mappedDat, a matrix of expression values (gene/reaction X Sample). Usually gene
%               expression values or gene expression values mapped onto a models
%               reactions using GPR rules.
%       lowerThres, the lower threshold for the localT2
%           algorithm as a percentile from 0 to 100 
%       upperThres, the upper threshold for the localT2
%           algorithm as a percentile from 0 to 100.
%
%   OUTPUTS:    
%       coreMat, the binary core set matrix (gene/reaction X Sample).
%
%   Reference (Richelle et al. 2019, doi: 10.1371/JOURNAL.PCBI.1007185)


if exist('lowerThres','var') && exist('upperThres','var') && lowerThres < 1 && upperThres < 1
    warning('The thresholds should be percentiles in the 0 to 100 range, you may be using the 0 to 1 range here.')
end
if ~exist('lowerThres','var')
    lowerThres = 25;
end
if ~exist('upperThres','var')
    upperThres = 75;
end

linData = reshape(mappedDat,[],1);
lowerThres = quantile(linData,lowerThres/100);
upperThres = quantile(linData,upperThres/100);

coreMat = false(size(mappedDat));
meanDat = mean(mappedDat,2);

allAbove = meanDat >= upperThres;
allBellow =  meanDat <= lowerThres;

coreMat(allAbove,:) = true;
coreMat(allBellow,:) = false;
for i = 1:size(coreMat,1)
    if allAbove(i) || allBellow(i)
        %Nothing
    else
        coreMat(i,:) = mappedDat(i,:)>= meanDat(i);
    end
end
end

