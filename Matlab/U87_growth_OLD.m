
% U87 growth

%% Load model
clear
load U87CbModel
clc

%% Find inactive reactions / dead ends
% inactiveRxns = checkModelConsistency(U87,[],'glpk');
% deadEndMets = detectDeadEnds_fast(U87);
% deadEndMets = U87.mets(deadEndMets);

%% Add glycerophospholipid metabolism reactions
rxnList = {'CLS_hs','PGPPT','PGPP_hs'};
for i = 1:numel(rxnList)
    indHR1 = strmatch(rxnList{i},model.rxns,'exact');
    rxnForm = printRxnFormula(model,rxnList{i});
    rxnRev = model.rev(indHR1);
    rxnLB = model.lb(indHR1);
    rxnUB = model.ub(indHR1);
    rxnSubSys = model.subSystems(indHR1);
    % ~~ GO BACK AND ADD GENE RULES ~~
    
    U87 = addReaction(U87,rxnList{i},rxnForm{:});
    indU87 = strmatch(rxnList{i},U87.rxns,'exact');
    U87.rev(indU87) = rxnRev;
    U87.lb(indU87) = rxnLB;
    U87.ub(indU87) = rxnUB;
    U87.subSystems(indU87) = rxnSubSys;
end

%% Add sphingolipid metabolism reactions
rxnList = {'SGPL11r','SPH1Ptr'};
for i = 1:numel(rxnList)
    indHR1 = strmatch(rxnList{i},model.rxns,'exact');
    rxnForm = printRxnFormula(model,rxnList{i});
    rxnRev = model.rev(indHR1);
    rxnLB = model.lb(indHR1);
    rxnUB = model.ub(indHR1);
    rxnSubSys = model.subSystems(indHR1);
    % ~~ GO BACK AND ADD GENE RULES ~~
    
    U87 = addReaction(U87,rxnList{i},rxnForm{:});
    indU87 = strmatch(rxnList{i},U87.rxns,'exact');
    U87.rev(indU87) = rxnRev;
    U87.lb(indU87) = rxnLB;
    U87.ub(indU87) = rxnUB;
    U87.subSystems(indU87) = rxnSubSys;
end
% add transport reaction
U87 = addReaction(U87,'ETHAMPtr',{'ethamp[c]','ethamp[r]'},[-1,1],1,...
    -1000,1000,0,'Transport, Endoplasmic Reticular');

%% Add missing exchanges from HR1
excIndHR1 = findExcRxns(model);
excRxnsHR1 = model.rxns(excIndHR1);

excIndU87 = findExcRxns(U87);
excRxnsU87 = U87.rxns(excIndU87);

[missingExcRxns,missingInd] = setdiff(excRxnsHR1,excRxnsU87);

for i = 1:length(missingExcRxns);
    excRxn_i = printRxnFormula(model,missingExcRxns{i});
    U87 = addReaction(U87,missingExcRxns{i},excRxn_i{:});
end
    
model = U87;
clc

%% Define biomass & add glycogenin primer
biomass_spec = read_biomass(model, 'biomass_3.txt');
if any(biomass_spec(match_regexp(model.metNames, 'glycogenin-')))
    % add glycogenin primer uptake (protein always available) to enable 
    % glycogen production
    model = add_reaction(model,...
        'UP_Tyr_DASH_ggn[c]',{'Tyr_DASH_ggn[c]'},-1,false,-1000,0);
end
biomassMets = model.mets(biomass_spec ~= 0);
biomassCoeff = biomass_spec(biomass_spec ~= 0);

%% Add biomass
model = addReaction(model,'biomass',biomassMets,biomassCoeff,0);
model = changeObjective(model,'biomass',1);

%% Add ATP maintenance reaction
model = addReaction(model,'ATPmaint(c)','atp[c] + h2o[c] -> adp[c] + h[c] + pi[c]');
model = changeRxnBounds(model,'ATPmaint(c)',0,'l');

%% Find inactive reactions / dead ends
inactiveRxns = checkModelConsistency(model,[],'glpk');
deadEndMets = detectDeadEnds_fast(model);
deadEndMets = model.mets(deadEndMets);

%% Clear model inputs
excInd = findExcRxns(model);
excRxns = model.rxns(excInd);
excLBs = model.lb(excInd);
excUBs = model.ub(excInd);

model = changeRxnBounds(model,excRxns,0,'l');

%% Define media
% media components:
% - oxygen (o2)
% - sodium (na1)
% - potassium (k)
% - calcium (ca2)
% - iron (fe2)
% - chlorine (cl)
% - phosphate (pi)
% - sulfate (so4)
% - ammonia (nh4)
% - pyruvate (pyr)

mediaList = {'EX_o2(e)','EX_na1(e)','EX_k(e)','EX_ca2(e)','EX_fe2(e)',...
    'EX_cl(e)','EX_pi(e)','EX_so4(e)','EX_nh4(e)','EX_pyr(e)'};
model = changeRxnBounds(model,mediaList,-1000,'l');

%% Add amino acids
% amino acids:
% - L-glutamine (gln)
% - L-isoleucine (ile)
% - L-leucine (leu)
% - L-lysine (lys)
% - L-tyrosine (tyr)
% * L-arginine (arg)
% * L-threonine (thr)
% * L-valine (val)

aaList = {'EX_gln_L(e)'};
% aaList = {'EX_gln_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
%     'EX_tyr_L(e)'};
% aaList = {'EX_gln_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
%     'EX_tyr_L(e)','EX_arg_L(e)','EX_thr_L(e)','EX_val_L(e)'};
model = changeRxnBounds(model,aaList,-1000,'l');

%% Adjust substrate uptake
% carbon substrates:
% - glucose (glc)

substrateList = {'EX_glc(e)'};
model = changeRxnBounds(model,substrateList,-1000,'l');

%% Optimize model
clc
changeCobraSolver('glpk');
FBAsolution = optimizeCbModel(model)
printFluxVector(model,FBAsolution.x,1,1)

%% Add temporary exchanges
tmpList = {'DM_cmp[c]','DM_dgmp[c]','DM_gmp[c]','DM_pe-hs[c]',...
    'DM_xolest-hs[c]'};
modelTest = changeRxnBounds(modelTest,tmpList,0,'l');
tmpGrow = false(numel(tmpList),1);
for i = 1:numel(tmpList)
    modelTest = changeRxnBounds(modelTest,tmpList{i},-1000,'l');
    FBAsolution = optimizeCbModel(modelTest);
    tmpGrow(i) = abs(FBAsolution.f) > 1e-6;
%     modelTest = changeRxnBounds(modelTest,tmpList{i},0,'l');
end

%%
modelTest = changeObjective(modelTest,'biomass',1);
FBAsolution = optimizeCbModel(modelTest)
printFluxVector(modelTest,FBAsolution.x,1,1)

%% Check exchanges
excInd = findExcRxns(modelTest);
excRxns = modelTest.rxns(excInd);
excLBs = modelTest.lb(excInd);
excUBs = modelTest.ub(excInd);

%% Check biomass components
clc
[modelTest,rxnNames] = addDemandReaction(model,biomassMets);

blockedBiomass = false(numel(rxnNames),1);
for i = 33%:numel(rxnNames)
    modelTest = changeObjective(modelTest,rxnNames{i},1);
    FBAsolution = optimizeCbModel(modelTest);
    blockedBiomass(i) = abs(FBAsolution.f) < 1e-6;
end

blockedBiomassMets = biomassMets(blockedBiomass);

%% Check biomass precursors
clc
[modelTest,precursorRxn] = addDemandReaction(modelTest,'pmtcoa[c]');
modelTest = changeObjective(modelTest,precursorRxn,1);
FBAsolution = optimizeCbModel(modelTest)

%% Check biomass components with precursors added
modelTest = changeRxnBounds(modelTest,precursorRxn,-1000,'l');
modelTest = changeObjective(modelTest,'DM_R1coa-hs[c]',1);
FBAsolution = optimizeCbModel(modelTest)
modelTest = changeRxnBounds(modelTest,precursorRxn,0,'l');

%% Turn off biomass component demands
dmndList = {'DM_asp-DASH-L[c]','DM_gly[c]','DM_h[c]','DM_clpn-hs[c]'};
dmndGrow = true(numel(dmndList),1);
for i = 1:numel(dmndList)
    modelTest = changeRxnBounds(modelTest,dmndList{i},0,'u');
    FBAsolution = optimizeCbModel(modelTest);
    dmndGrow(i) = abs(FBAsolution.f) < 1e-6;
    modelTest = changeRxnBounds(modelTest,dmndList{i},1000,'u');
end
    
%% Find inactive reactions / dead ends
inactiveRxns = checkModelConsistency(modelTest,[],'glpk');
deadEndMets = detectDeadEnds_fast(modelTest);
deadEndMets = modelTest.mets(deadEndMets);
% tic;[nonConsumed,nonProduced] = findDeadEnds(modelTest,modelTest.mets);toc
tic;[nonConsumed,nonProduced] = findDeadEnds(modelTest,deadEndMets);toc
nonProdOnly = setdiff(nonProduced,nonConsumed);

%% Examine reaction participation for a metabolite
metSubSystems = {''};
for m = 1:numel(blockedBiomassMets)
    met = blockedBiomassMets{m};
    rxnInfo = modelSubset(modelTest,met);
    metSubSystems = union(metSubSystems,unique(rxnInfo(:,2)));
end
del = union(strmatch('Transport',metSubSystems),...
    strmatch('Demand',metSubSystems));
metSubSystems([1;del(:)]) = [];

%% Check all subsystems for dead ends
subSysDead = false(numel(metSubSystems),1);
for i = 1:numel(metSubSystems)
    display([num2str(i),'/',num2str(numel(metSubSystems))]);
    subSys = metSubSystems{i};
    subModel = extractSubSysModel(modelTest,subSys);
    [nonConsumed,nonProduced] = findDeadEnds(modelTest,subModel.mets);
    subSysDead(i) = numel(nonConsumed) | numel(nonProduced);
end

%% Find dead ends within a subsystem
clc
subSys = 'Sphingolipid Metabolism';
subModel = extractSubSysModel(modelTest,subSys);
[nonConsumed,nonProduced] = findDeadEnds(modelTest,subModel.mets);
deadEndMets = union(nonConsumed,nonProduced);
subSys