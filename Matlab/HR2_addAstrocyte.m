
% Add astroctye reactions to Human Recon 2
clear
clc
load GLU_norm
% load HR2_CbModel_new
% load 11_07_08_Recon20
load HR2_CbModel_Dec2012
Recon2 = model;

%% Pull out cytosolic and mitochondrial reactions from astrocyte
cA_model = extractCompModel_JE(GLU_model,'cA',0);
mA_model = extractCompModel_JE(GLU_model,'mA',0);

cA_rxns = cA_model.rxns;
mA_rxns = mA_model.rxns;
HR2_rxns = Recon2.rxns;

A_rxns = union(mA_rxns,cA_rxns);

% Fix formatting issue with Recon2 reaction IDs
% for i = 1:138
% HR2_rxns{i}(1) = [];
% end

%% Format reaction IDs to match Recon2
A_rxns_fmt = regexprep(A_rxns,',','_');
A_rxns_fmt = regexprep(A_rxns_fmt,'\-','_');
A_rxns_fmt = regexprep(A_rxns_fmt,'\(','_');
A_rxns_fmt = regexprep(A_rxns_fmt,'\)','_');
A_rxns_fmt = regexprep(A_rxns_fmt,'_Int','');
[Aonly_rxns_fmt,diff_idx] = setdiff(A_rxns_fmt,HR2_rxns);

%%
Aonly_rxns = A_rxns(diff_idx);
[rxnsInt,Aidx,GLUidx] = intersect(Aonly_rxns,GLU_model.rxns);
Aonly_lbs = GLU_model.lb(GLUidx);
Aonly_ubs = GLU_model.ub(GLUidx);
Aonly_c = GLU_model.c(GLUidx);
Aonly_subSys = GLU_model.subSystems(GLUidx);
Aonly_grRules = GLU_model.grRules(GLUidx);
Aonly_annotations = GLU_model.annotations(GLUidx);

%% Format subsystems
for i = 1:length(Aonly_subSys)
    Aonly_subSys{i} = ['',Aonly_subSys{i}];
end

%% Format annotations as reaction notes
for i = 1:length(Aonly_annotations)
    if numel(Aonly_annotations{i})
        Aonly_annotations{i} = ['NOTES: ',Aonly_annotations{i}];
    else Aonly_annotations{i} = '';
    end
end

%% Format gene rules
Aonly_grRules_fmt = repmat({''},size(Aonly_grRules));
for i = 1:numel(Aonly_grRules)
    if numel(Aonly_grRules{i})
        ruleParts = regexp(Aonly_grRules{i},'\s','split');
        ruleParts = regexprep(ruleParts,'\(|\)','');
        boolean = regexp(ruleParts,'and|or|\(|\)');
        openParen = repmat({'('},1,numel(ruleParts));
        closeParen = repmat({'.1)'},1,numel(ruleParts));
        openParen(cellfun(@length,boolean)>0) = {''};
        closeParen(cellfun(@length,boolean)>0) = {''};
        ruleParts = strcat(openParen,ruleParts,closeParen);
        ruleParts = strcat(ruleParts,repmat({'_'},1,numel(ruleParts)));
        rule = strcat(ruleParts{:});
        rule = regexprep(rule,'_',' '); rule(end) = '';
        Aonly_grRules_fmt{i} = rule;
    end
end

%% Print formulas for astrocyte reactions
rxnFormulas = repmat({''},numel(Aonly_rxns),1);
for i = 1:numel(Aonly_rxns)
    rxnFormulas(i) = printRxnFormula(GLU_model,Aonly_rxns{i});
end

%% Format reaction formulas
rfsNew = regexprep(rxnFormulas,'\[cA\]','[c]');
rfsNew = regexprep(rfsNew,'\[mA\]','[m]');
rfsNew = regexprep(rfsNew,'\[I\]','[e]');
rfsNew = regexprep(rfsNew,'\-','_');
rfsNew = regexprep(rfsNew,'_>','->');

%% Check metabolites unique to astrocyte
% A_mets = {};
% for i = 1:numel(rfsNew)
%     A_mets = union(A_mets,parseRxnFormula(rfsNew{i}));
% end
% A_mets = A_mets';
% 
% Aonly_mets = setdiff(A_mets,Recon2.mets);

%% Add reactions to Recon2
added = false(size(Aonly_rxns));
Recon2_plusA = Recon2;
for i = 1:numel(rfsNew)
    [metList,stoichList,revFlag] = parseRxnFormula(rfsNew{i});
    lb = Aonly_lbs(i);
    ub = Aonly_ubs(i);
    objCoeff = Aonly_c(i);
    subSystem = Aonly_subSys{i};
    grRule = Aonly_grRules_fmt{i};
    [Recon2_plusA,exists] = addReaction_JE(Recon2_plusA,...
        Aonly_rxns_fmt{i},metList,stoichList,revFlag,...
        lb,ub,objCoeff,subSystem,grRule);
    
%     [Recon2_plusA,exists] = addReaction(Recon2_plusA,Aonly_rxns_fmt{i},rfsNew{i});
    added(i) = numel(exists)==0;
end

%% Extract info for added reactions
plusA_rxns = Aonly_rxns_fmt(added);
plusA_rxnFormulas = rfsNew(added);
plusA_grRules = Aonly_grRules_fmt(added);
plusA_lb = Aonly_lbs(added);
plusA_ub = Aonly_ubs(added);
plusA_subSys = Aonly_subSys(added);
plusA_annotations = Aonly_annotations(added);

%% Make sure meta info is included for added reactions
[rxnInt,HR2pAidx,plusAidx] = intersect(Recon2_plusA.rxns,plusA_rxns);
Recon2_plusA.rxnNotes(HR2pAidx) = plusA_annotations(plusAidx);
Recon2_plusA.subSystems(HR2pAidx) = plusA_subSys(plusAidx);
Recon2_plusA.confidenceScores(HR2pAidx) = repmat({''},size(plusAidx));
% Recon2_plusA.rxnConfidenceEcoIDA(HR2pAidx) = repmat({''},size(plusAidx));
% Recon2_plusA.rxnsboTerm(HR2pAidx) = repmat({''},size(plusAidx));

%% Make sure molecular formulas are included for added metabolites
[plusA_mets,HR2pAidx] = setdiff(Recon2_plusA.mets,Recon2.mets);
plusA_mets = regexprep(plusA_mets,'\]','A]');
plusA_mets = regexprep(plusA_mets,'_','-');
[metInt,GLUidx,plusAidx] = intersect(GLU_model.mets,plusA_mets);
plusA_formulas = GLU_model.metFormulas(GLUidx);

Recon2_plusA.metFormulas(HR2pAidx) = plusA_formulas;
% Recon2_plusA.metNeutralFormula(HR2pAidx) = {'C10H14N5O6P';'';''};
Recon2_plusA.metCharge(HR2pAidx) = [-2;0;0];
Recon2_plusA.metHMDBID(HR2pAidx) = {'';'';''};
Recon2_plusA.metEHMNAbbrevs(HR2pAidx) = {'';'';''};
Recon2_plusA.metHepatonetAbbrevs(HR2pAidx) = {'';'';''};
Recon2_plusA.metKEGGID(HR2pAidx) = {'';'';''};
Recon2_plusA.metInChIString(HR2pAidx) = {'';'';''};
% Recon2_plusA.metCompartment(HR2pAidx) = {'m';'m';'m'};

Recon2_plusA.metFormulas(5065:5066) = {''};


%% Update description
Recon2_plusA.description = 'Recon2(IT)_plus_Astrocyte(NL)';

%%
save('HR2plusA_info')
save('HR2plusA_CbModel','Recon2_plusA','GLU_model','Recon2')