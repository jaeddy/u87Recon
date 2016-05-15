
% U87 growth

%% Load model
clear
clc
load HR2plusA_CbModel
model = Recon2_plusA;

%% Add missing reactions

% This transport is required for growth

rxnName = 'GM1t';
rxnFormula = 'gm1_hs[g] -> gm1_hs[c]';
model = addReaction(model,rxnName,rxnFormula);
model = changeRxnBounds(model,rxnName,0,'l');

metIdx = strmatch('gm1_hs[c]',model.mets);
model.metCharge(metIdx) = -1;
model.metNames{metIdx} = 'ganglioside GM1';
model.metFormulas{metIdx} = 'C55H95N3O30FULLRCO';
model.metChEBIID{metIdx} = '18216';
model.metKEGGID{metIdx} = '';
model.metPubChemID{metIdx} = '';
model.metHMDBID{metIdx} = '';
model.metInChIString{metIdx} = '';
model.metEHMNAbbrevs{metIdx} = '';
model.metHepatonetAbbrevs{metIdx} = '';
% 
rxnIdx = strmatch(rxnName,model.rxns);
rxnNotes{rxnIdx} = 'NOTES:cytosolic GM1 needed for growth; mechanism unknown JE';
rxnNames{rxnIdx} = 'ganglioside GM1 transport';

% This transport is required for growth

rxnName = 'Tyr_ggnt';
rxnFormula = 'Tyr_ggn[e] -> Tyr_ggn[c]';
model = addReaction(model,rxnName,rxnFormula);
model = changeRxnBounds(model,rxnName,0,'l');

% This reaction is required for growth?

% rxnName = 'ALKRS';
% grRule = '';
% subSystem = '';
% rxnFormula = 'Rtotal[x] -> alkylR1oh[x]';
% model = addReaction(model,rxnName,rxnFormula);
% model = changeRxnBounds(model,rxnName,0,'l');

%% Define & add U87 biomass equation
biomass_spec = read_biomass(model, 'U87_biomass.txt');
biomassMets = model.mets(biomass_spec ~= 0);
biomassCoeff = biomass_spec(biomass_spec ~= 0);

% dak2gpe(11), pchol_hs(24), glygn2(35), ps_hs(45), xolest_hs(48)
biomassMets([11]) = [];
biomassCoeff([11]) = [];

% Add biomass
model = addReaction(model,'biomass_U87',biomassMets,biomassCoeff,0);
model = changeObjective(model,'biomass_U87',1);
model = changeRxnBounds(model,'biomass_reaction',0,'b');

%% Add cytosolic ATP maintenance reaction
% ~~ Recon2 currently has mitochondrial ATP maintenance, so skip for now
% model = addReaction(model,'DM_atp_c_','atp[c] + h2o[c] -> adp[c] + h[c] + pi[c]');
% model = changeRxnBounds(model,'DM_atp_c_',0,'l');

%% Select model inputs
excInd = findExcRxns(model);
excRxns = model.rxns(excInd);
excLBs = model.lb(excInd);
excUBs = model.ub(excInd);

excRxns_nonFixed = excRxns(rem(excLBs,1)==0 & excLBs~=0);
excRxns_fixed = excRxns(rem(excLBs,1)~=0);
model = changeRxnBounds(model,excRxns_fixed,0,'b');

excRxns_large = excRxns(rem(excLBs,10)==0);
model = changeRxnBounds(model,excRxns_large,0,'l');

%% Clear model inputs 
model = changeRxnBounds(model,excRxns_nonFixed,0,'l');
model = changeRxnBounds(model,'EX_Tyr_ggn(e)',-1000,'l');

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

% aaList = {'EX_gln_L(e)'};

% DMEM list (restrictive)
% aaList = {'EX_gln_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
%     'EX_tyr_L(e)'};

% DMEM list (less restrictive)
% aaList = {'EX_gln_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
%     'EX_tyr_L(e)','EX_arg_L(e)','EX_thr_L(e)','EX_val_L(e)'};

% DMEM list (full)
aaList = {'EX_arg_L(e)','EX_cys_L(e)','EX_gln_L(e)','EX_gly(e)',...
    'EX_his_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
    'EX_met_L(e)','EX_phe_L(e)','EX_ser_L(e)','EX_thr_L(e)',...
    'EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)'};

% Essential AAs
% aaList = {'EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
%     'EX_met_L(e)','EX_phe_L(e)','EX_thr_L(e)','EX_trp_L(e)',...
%     'EX_val_L(e)','EX_his_L(e)'};
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
% printFluxVector(model,FBAsolution.x,1,1)

% %% Check biomass components
% clc
% [modelTest,rxnNames] = addDemandReaction(model,biomassMets);
% rxnNames{23} = 'DM_pchol_hs';
% 
% blockedBiomass = false(numel(rxnNames),1);
% for i = 1:numel(rxnNames)
%     modelTest = changeObjective(modelTest,rxnNames{i},1);
%     FBAsolution = optimizeCbModel(modelTest);
%     blockedBiomass(i) = abs(FBAsolution.f) < 1e-6;
% end
% 
% blockedBiomassMets = biomassMets(blockedBiomass);
% 
% %% Trace precursors
% clc
% met = 'Tyr_ggn[e]';
% rxns = findRxnsFromMets(model,met);
% printRxnFormula(model,rxns);
% 
% %% Check biomass precursors
% met = 'Tyr_ggn[c]';
% [modelTest,precursorRxn] = addDemandReaction(model,met);
% modelTest = changeObjective(modelTest,precursorRxn,1);
% FBAsolution = optimizeCbModel(modelTest);
% result = {'Success','Fail'};
% display([met,'? ',result{2-double(FBAsolution.f > 0)},': ',...
%     num2str(FBAsolution.f)])
% 
% %% Add false uptakes
% checkList = {'Ex_galdag_hs[c]','Ex_dak2gpc_hs[c]','Ex_dak2gpe_hs[c]'};
% model = changeRxnBounds(model,checkList,-1000,'l');
% % model = changeRxnBounds(model,checkList,1000,'u');

