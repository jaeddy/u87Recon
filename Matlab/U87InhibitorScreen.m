
% hockenbery screen validation
clear
clc

load u87inhibitors
% load HR2_CbModel_Dec2012
load U87_CbModel
model = U87;

sum(ismember(ecNum,model.rxnECNumbers))

%% Filter out inhibitors with poor Z-factors

zCutoff = 'moderate';

switch zCutoff
    case 'good'
        cutoff = 0.5;
    case 'moderate'
        cutoff = 0;
end

keepIdx = zfactor >= cutoff;
ecNum = ecNum(keepIdx);
brendaEnzyme = brendaEnzyme(keepIdx);
percSurvival = percSurvival(keepIdx);
gene = gene(keepIdx);
        
%% Parse inhibitor target ECs, identify matching reactions

ecNumClean = regexprep(ecNum,'\s(\(|\w|\)|\.|=|/)*','');
sum(ismember(ecNumClean,model.rxnECNumbers));

ecMatch = {};
for i = 1:numel(ecNumClean)
    if ~isempty(ecNumClean{i})
        ecMatch{i} = strmatch(ecNumClean{i},model.rxnECNumbers,'exact');
    end
end

matchLength = cellfun('length',ecMatch);
rxnMatch = repmat({''},numel(ecMatch),max(matchLength));
rxnNameMatch = repmat({''},numel(ecMatch),max(matchLength));
for i = 1:numel(ecMatch)
    if matchLength(i) > 0
        rxnMatch(i,1:matchLength(i)) = model.rxns(ecMatch{i});
        rxnNameMatch(i,1:matchLength(i)) = model.rxnNames(ecMatch{i});
    end
end


%% Calculate WT growth yield
modelDMEM = formulateDMEM(U87,3);
FBAsolution = optimizeCbModel(modelDMEM);
biomassWT = FBAsolution.f;


%% Knockout reactions corresponding to inhibitor targets
numKOs = size(rxnMatch,1);

inhibitorResult = zeros(numKOs,1);
for i = 1:numKOs
    if matchLength(i) > 0
        inhibitorTargets = rxnMatch(i,1:matchLength(i));
        modelTmp = changeRxnBounds(modelDMEM,inhibitorTargets,0,'b');
        FBAsolution = optimizeCbModel(modelTmp);
        inhibitorResult(i) = FBAsolution.f/biomassWT;
    else inhibitorResult(i) = -1;
    end
end
numAffected = sum(abs(inhibitorResult) < 1)    

%% Compile results
matchedIdx = inhibitorResult > -1;

matchedExp = percSurvival(matchedIdx);
matchedSim = inhibitorResult(matchedIdx);

growthPerc = 70;

TP = sum(matchedExp >= growthPerc & matchedSim >= growthPerc/100);
FP = sum(matchedExp < growthPerc & matchedSim >= growthPerc/100);
TN = sum(matchedExp < growthPerc & matchedSim < growthPerc/100);
FN = sum(matchedExp >= growthPerc & matchedSim < growthPerc/100);

sensitivity = TP/(TP+FN);
specificity = TN/(TN+FP);

