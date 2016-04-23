
% compile U87 mCADRE results

%% Load mCADRE results for U87
load U87mCADREResults


%% Identify blocked reactions in the pruned model
% inactivePM = checkModelConsistency(PM,[],{},1);

%% Add blocked core reactions back into the pruned model
blockedCore = intersect(model_C,inactiveRxns);
blockedCoreFormulas = printRxnFormula(model,blockedCore,0);

[~,blockedIdx,modelIdx] = intersect(blockedCore,model.rxns);
blockedCoreLB = zeros(size(blockedCore));
blockedCoreLB(blockedIdx) = model.lb(modelIdx);

blockedCoreUB = zeros(size(blockedCore));
blockedCoreUB(blockedIdx) = model.ub(modelIdx);

blockedCoreModel = extractSubNetwork(model,blockedCore);

%%
U87 = PM;
for i = 1:numel(blockedCore)
    U87 = addReaction(U87,blockedCore{i},blockedCoreFormulas{i});
    rxnIdx = strmatch(blockedCore{i},U87.rxns);
    U87.lb(rxnIdx) = blockedCoreLB(i);
    U87.ub(rxnIdx) = blockedCoreUB(i);
end


%% Identify blocked reactions in the final U87 model
% inactiveU87 = checkModelConsistency(U87,[],{},1);


%% Collect and update all metadata for model reactions and metabolites

% Add any missing fields to the U87 model
modelFields = fieldnames(model);
U87Fields = fieldnames(U87);

modelDims = [size(model.S,1); size(model.S,2)];
U87Dims = [size(U87.S,1); size(U87.S,2)];

missingFields = setdiff(modelFields,U87Fields);
for i = 1:numel(missingFields)
    fieldSize = size(getfield(model,missingFields{i}));
    newSize = U87Dims(ismember(modelDims,fieldSize));
    
    fieldClass = class(getfield(model,missingFields{i}));
    switch fieldClass 
        case 'double'
            U87 = setfield(U87,missingFields{i},zeros(newSize,1));
        case 'cell'
            U87 = setfield(U87,missingFields{i},repmat({''},newSize,1));
    end
end
U87 = orderfields(U87,model);       

% Locate and extract relevant info from Recon 2
[rxnsInt,modelIdx,U87Idx] = intersect(model.rxns,U87.rxns);
rxnFields = modelFields(~cellfun('isempty',regexp(modelFields,'rxn(?!s)')));
for i = 1:numel(rxnFields)
    tmpModelField = getfield(model,rxnFields{i});
    tmpU87Field = getfield(U87,rxnFields{i});
    tmpU87Field(U87Idx) = tmpModelField(modelIdx);
    U87 = setfield(U87,rxnFields{i},tmpU87Field);
end
    
[metsInt,modelIdx,U87Idx] = intersect(model.mets,U87.mets);
metFields = modelFields(~cellfun('isempty',regexp(modelFields,'met(?!s)')));
for i = 1:numel(metFields)
    tmpModelField = getfield(model,metFields{i});
    tmpU87Field = getfield(U87,metFields{i});
    tmpU87Field(U87Idx) = tmpModelField(modelIdx);
    U87 = setfield(U87,metFields{i},tmpU87Field);
end

% Update GPRs
for i = 1:numel(blockedCore)
    rxnName = blockedCore{i};
    rxnIdx = find(~cellfun('isempty',regexp(model.rxns,rxnName)));
    grRule = model.grRules{rxnIdx};
    U87 = changeGeneAssociation(U87,rxnName,grRule);
end
    

% Extract any remaining info from the primed input model (added reactions)

%% Compile model statistics

% Number of reactions (w/ and w/o blocked core)
HR2rxns = numel(model.rxns);
GMrxns = numel(GM.rxns);
PMrxns = numel(PM.rxns);
bCrxns = numel(blockedCoreModel.rxns);
U87rxns = numel(U87.rxns);
rxnsRow = cellstr(num2str([HR2rxns,GMrxns,PMrxns,bCrxns,U87rxns]'))';

% - non-exchange/transport reactions?
% - unique reactions?

% Number of metabolites (w/ and w/o blocked core)
HR2mets = numel(model.mets);
GMmets = numel(GM.mets);
PMmets = numel(PM.mets);
bcMets = numel(blockedCoreModel.mets);
U87mets = numel(U87.mets);
metsRow = cellstr(num2str([HR2mets,GMmets,PMmets,bcMets,U87mets]'))';

% - unique metabolites?
HR2mets_u = numel(unique(regexprep(model.mets,'\[\w\]','')));
GMmets_u = numel(unique(regexprep(GM.mets,'\[\w\]','')));
PMmets_u = numel(unique(regexprep(PM.mets,'\[\w\]','')));
bCmets_u = numel(unique(regexprep(blockedCoreModel.mets,'\[\w\]','')));
U87mets_u = numel(unique(regexprep(U87.mets,'\[\w\]','')));
mets_uRow = cellstr(num2str([HR2mets_u,GMmets_u,PMmets_u,bCmets_u,U87mets_u]'))';

% Number of genes (w/ and w/o blocked core)
HR2geneList = model.genes(sum(model.rxnGeneMat)>0);
GMgeneList = GM.genes(sum(GM.rxnGeneMat)>0);
PMgeneList = PM.genes(sum(PM.rxnGeneMat)>0);
bCgeneList = PM.genes(sum(blockedCoreModel.rxnGeneMat)>0);
U87geneList = U87.genes(sum(U87.rxnGeneMat)>0);

HR2xcripts = numel(HR2geneList);
GMxcripts = numel(GMgeneList);
PMxcripts = numel(PMgeneList);
bCxcripts = numel(bCgeneList);
U87xcripts = numel(U87geneList);
xcriptsRow = cellstr(num2str([HR2xcripts,GMxcripts,PMxcripts,bCxcripts,U87xcripts]'))';

HR2genes = numel(unique(regexprep(HR2geneList,'\.[0-9]*','')));
GMgenes = numel(unique(regexprep(GMgeneList,'\.[0-9]*','')));
PMgenes = numel(unique(regexprep(PMgeneList,'\.[0-9]*','')));
bCgenes = numel(unique(regexprep(bCgeneList,'\.[0-9]*','')));
U87genes = numel(unique(regexprep(U87geneList,'\.[0-9]*','')));
genesRow = cellstr(num2str([HR2genes,GMgenes,PMgenes,bCgenes,U87genes]'))';

% Number of blocked reactions (w/ and w/o blocked core)

% Create table
statsTable = [rxnsRow;...
    metsRow;...
    mets_uRow;...
    xcriptsRow;...
    genesRow];

%% Compartment distribution
compartments = {'e','c','m','n','r','x','l','g'};
modelNames = {'model','GM','PM','blockedCoreModel','U87'};

compMetsTable = repmat({''},numel(compartments),numel(modelNames));
for i = 1:numel(compartments)
    for j = 1:numel(modelNames)
        compModel = eval(['extractCompModel_JE(',...
            modelNames{j},',compartments{i},0);']);
        numMets = sum(~cellfun('isempty',regexp(compModel.mets,...
            ['\[',compartments{i},'\]'])));
        compMetsTable{i,j} = num2str(numMets);
    end
end

%% Pathway distribution
subSysVec = zeros(numel(model.rxns),1);
for i = 1:numel(subSystems)
    subSysIdx = ismember(model.subSystems,subSystems{i});
    subSysVec(subSysIdx) = subSystemCodes(i);
end

groups = {[5],[1],[2],[7],[3],[4],[6,9,10],[8]};
groupVec = zeros(numel(groups),1);
for i = 1:numel(groups)
    groupVec(i) = sum(ismember(subSysVec,groups{i}))/numel(subSysVec);
end
%% Test model growth


% Compartment distribution



%% Test model growth