

% prime Recon 2 for U87 mCADRE
clear
clc

addRxns = 1;
useU87biomass = 1;
resetExcs = 1;
checkBiomass = 0;

% Load model
% load HR2_CbModel_Dec2012
% load Recon2.v02.mat
% model = modelRecon2beta121114_fixed;
load HR2v03_CbModel_Jan2014

%% Add missing reactions
if addRxns
    
% This reaction is needed to syntheisze NAD from Trp (was in Recon 1, but
% has been removed for some reason)

% rxnName = 'NNDPR';
% rxnFormula = '2 h[c] + prpp[c] + quln[c] -> co2[c] + nicrnt[c] + ppi[c]';
% metaboliteList = {'h[c]','prpp[c]','quln[c]','co2[c]','nicrnt[c]','ppi[c]'};
% stoichCoeffList = [-2,-1,-1,1,1,1];
% revFlag = 0;
% lowerBound = 0;
% upperBound = 1000;
% objCoeff = 0;
% subSystem = 'NAD metabolism';
% grRule = '(23474.1)';
% model = addReaction(model,rxnName,metaboliteList,stoichCoeffList,revFlag,...
%     lowerBound,upperBound,objCoeff,subSystem,grRule);
% model.ExchRxnBool(end+1) = 0;
% model.EXRxnBool(end+1) = 0;
% model.DMRxnBool(end+1) = 0;
% model.SinkRxnBool(end+1) = 0;
% model.SIntRxnBool(end+1) = 1;
% 
% % This reaction is needed to synthesize mitochondrial NADP from NAD
% 
% rxnName = 'NADKm';
% rxnFormula = 'atp[m] + nad[m] -> adp[m] + h[m] + nadp[m]';
% metaboliteList = {'atp[m]','nad[m]','adp[m]','h[m]','nadp[m]'};
% stoichCoeffList = [-1,-1,1,1,1];
% revFlag = 0;
% lowerBound = 0;
% upperBound = 1000;
% objCoeff = 0;
% subSystem = 'NAD metabolism';
% grRule = '(133686.1)'; % C5orf33
% model = addReaction(model,rxnName,metaboliteList,stoichCoeffList,revFlag,...
%     lowerBound,upperBound,objCoeff,subSystem,grRule);
% rxnIdx = strmatch(rxnName,model.rxns);
% model.rxnNotes{rxnIdx} = 'NOTES:See PMID: 23212377 - human mitochondrial NAD kinase JE';
% model.ExchRxnBool(end+1) = 0;
% model.EXRxnBool(end+1) = 0;
% model.DMRxnBool(end+1) = 0;
% model.SinkRxnBool(end+1) = 0;
% model.SIntRxnBool(end+1) = 1;

% This transport (primed glycogenin) is required for growth when glycogen
% is in biomass

rxnName = 'Tyr_ggnt';
rxnFormula = 'Tyr_ggn[e] -> Tyr_ggn[c]';
model = addReaction(model,rxnName,rxnFormula);
model = changeRxnBounds(model,rxnName,0,'l');

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
rxnIdx = strmatch(rxnName,model.rxns);
rxnNotes{rxnIdx} = 'NOTES:cytosolic GM1 needed for growth; mechanism unknown JE';
rxnNames{rxnIdx} = 'ganglioside GM1 transport';
model.ExchRxnBool(end+1) = 0;
model.EXRxnBool(end+1) = 0;
model.DMRxnBool(end+1) = 0;
model.SinkRxnBool(end+1) = 0;
model.SIntRxnBool(end+1) = 1;

end

%% Define & add U87 biomass equation
if useU87biomass
    biomass_spec = read_biomass(model, 'U87_biomass.txt');
    biomassMets = model.mets(biomass_spec ~= 0);
    biomassCoeff = biomass_spec(biomass_spec ~= 0);
    
    % Remove components
    biomassCoeff(strmatch('dak2gpe_hs[c]',biomassMets)) = [];
    biomassMets(strmatch('dak2gpe_hs[c]',biomassMets)) = [];
    
    biomassCoeff(strmatch('tag_hs[c]',biomassMets)) = [];
    biomassMets(strmatch('tag_hs[c]',biomassMets)) = [];
%     
%     biomassCoeff(strmatch('xolest_hs[c]',biomassMets)) = [];
%     biomassMets(strmatch('xolest_hs[c]',biomassMets)) = [];
%     
%     biomassCoeff(strmatch('gm1_hs[c]',biomassMets)) = [];
%     biomassMets(strmatch('gm1_hs[c]',biomassMets)) = [];
% 
%     biomassCoeff(strmatch('mag_hs[c]',biomassMets)) = [];
%     biomassMets(strmatch('mag_hs[c]',biomassMets)) = [];

    % Add biomass
    model = addReaction(model,'biomass_U87',biomassMets,biomassCoeff,0);
    model = changeObjective(model,'biomass_U87',1);
    model = changeRxnBounds(model,'biomass_reaction',0,'b');
else biomassMets = model.mets(...
        model.S(:,strmatch('biomass_reaction',model.rxns))~=0);
end

%% Open all inactivated reactions
revOff = model.lb == 0 & model.ub == 0 & model.rev == 1;
irrevOff = model.lb == 0 & model.ub == 0 & model.rev == 0;

model.lb(revOff) = -1000;
model.ub(revOff | irrevOff) = 1000;

%% Add cytosolic ATP maintenance reaction
% ~~ Recon2 currently has mitochondrial ATP maintenance, so skip for now
% model = addReaction(model,'DM_atp_c_','atp[c] + h2o[c] -> adp[c] + h[c] + pi[c]');
% model = changeRxnBounds(model,'DM_atp_c_',0,'l');

%%
if resetExcs > 0
    %~~~ Select model inputs
    % Recon2.02 has a few handy boolean vectors that make this step easy,
    % without relying on the 'findExcRxns' function
    excInd = ~model.SIntRxnBool; %findExcRxns(model);
    excInd(950) = 0; % ATP maintenance is cytosolic
    excRxns = model.rxns(excInd);
    excLBs = model.lb(excInd);
    excUBs = model.ub(excInd);
    
%     excFormulas = printRxnFormula(model,excRxns,0);
%     excMets = strtrim(regexprep(excFormulas,'(->|<=>)',''));
%     excMetFormulas = model.metFormulas(directIntersect(model.mets,excMets));
%     hasX = ~cellfun('isempty',strfind(excMetFormulas,'X'));
%     excRxns_hasX = excRxns(hasX);
%     excLBs_hasX = excLBs(hasX);
    
    %~~~ Clear model inputs 
    % I'm really clearing all non-internal reactions here, including
    % exchanges, sinks, and demands.
    
    model = changeRxnBounds(model,excRxns,0,'l');
        
    % Reset bounds of exchanges for X-containing metabolites
%     excLBs_hasX = -1000*ones(size(excLBs_hasX));
%     model = changeRxnBounds(model,excRxns_hasX,excLBs_hasX,'l');
    
%     model = changeRxnBounds(model,'EX_Tyr_ggn(e)',-1000,'l');
    
    %~~~ Define media
    
    % Essential minerals:
    % - oxygen (o2)
    % x- sodium (na1)
    % x- potassium (k)
    % x- calcium (ca2)
    % x- iron (fe2)
    % x- chlorine (cl)
    % x- phosphate (pi)
    % x- iodine (iodine)
    % - zinc (HC02172)

    mineralList = {'EX_o2(e)','EX_na1(e)','EX_k(e)','EX_ca2(e)',...
        'EX_fe2(e)','EX_cl(e)','EX_pi(e)'};
    model = changeRxnBounds(model,mineralList,-10,'l');

    % Essential amino acids:
    % - L-glutamine (gln)
    % - L-isoleucine (ile)
    % - L-leucine (leu)
    % - L-lysine (lys)
    % - L-tyrosine (tyr)
    % * L-arginine (arg)
    % * L-threonine (thr)
    % * L-valine (val)

    aaList = {'EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
        'EX_met_L(e)','EX_phe_L(e)','EX_thr_L(e)','EX_trp_L(e)',...
        'EX_val_L(e)','EX_his_L(e)'}; 
%     aaList = {'EX_arg_L(e)','EX_cys_L(e)','EX_gln_L(e)','EX_gly(e)',...
%         'EX_his_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
%         'EX_met_L(e)','EX_phe_L(e)','EX_ser_L(e)','EX_thr_L(e)',...
%         'EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)'};
    model = changeRxnBounds(model,aaList,-10,'l');

    % Essential fatty acids:
    % - linoleic acid (lnlc)
    % - alpha-linolenic acid (lnlnca)
    
    faList = {'EX_lnlc(e)','EX_lnlnca(e)'};
    model = changeRxnBounds(model,faList,-10,'l');
    
    % Essential vitamins:
    % - retinol / vitamin A (retinol)
    % - pantohenate / vitamin B5 (pnto_R)
    % - riboflavin / vitamin B2 (ribflv)
    % - choline / vitamin Bp (chol)
    % - thiamin / vitamin B1 (thm)
    % - niacin / vitamin B3 (NA)
    % - pyridoxine / vitamin B6 (pydxn)
    % - biotin / vitamin B7 (btn)
    % - folate / vitamin B9 (fol)
    % - adenosylcobalamin / vitamin B12 (adocbl) ** NO EXCHANGE **
    % - ascorbate / vitamin C (ascb_L)
    % - ergocalciferol / vitamin D (vitd2)
    % - alpha-tocopherol / vitamin E (avite1)
    % - beta-tocopherol / vitamin E (bvite)
    % - phylloquinone / vitamin K (phyQ)

    vitList = {'EX_pnto_R(e)','EX_ribflv(e)',...
        'EX_retinol(e)','EX_chol(e)','EX_thm(e)','EX_pydxn(e)',...
        'EX_btn(e)','EX_fol(e)','EX_adocbl(e)','EX_ascb_L(e)',...
        'EX_vitd2(e)','EX_avite1(e)','EX_bvite(e)','EX_phyQ(e)'};
    model = changeRxnBounds(model,vitList,-10,'l');
    
    %~~~ Adjust substrate uptake
    % carbon substrates:
    % - glucose (glc)

    substrateList = {'EX_glc(e)'};
%     substrateList = [substrateList,...
%         {'EX_tag_hs(e)'}];
    model = changeRxnBounds(model,substrateList,-10,'l');    
end

%% Optimize model

changeCobraSolver('glpk');
% FBAsolution = optimizeCbModel(model,'','one');
FBAsolution = optimizeCbModel(model);
printFluxVector(model,FBAsolution.x,1,1,-1,[],[],0)

%% Check biomass components
% checkBiomass = 1;
if checkBiomass
    clc
    [modelTest,rxnNames] = addDemandReaction(model,biomassMets);
    exMet = 'Rtotal2coa[c]';
    modelTest = addExchangeRxn(modelTest,{exMet},{['EX_',exMet]},-1000,1000);
%     rxnNames{strmatch('DM_pchol_hs[c]',rxnNames)} = 'DM_pchol_hs';

    if ~useU87biomass
        rxnNames{strmatch('DM_dgtp[n]',rxnNames)} = 'DM_dgtp_n_';
        rxnNames{strmatch('DM_dctp[n]',rxnNames)} = 'DM_dctp_n_';
        rxnNames{strmatch('DM_datp[n]',rxnNames)} = 'DM_datp_n_';
        rxnNames{strmatch('DM_dttp[n]',rxnNames)} = 'DM_dttp_n_';
    end

    blockedBiomass = false(numel(rxnNames),1);
    for i = 1:numel(rxnNames)
        modelTest = changeObjective(modelTest,rxnNames{i},1);
        FBAsolution = optimizeCbModel(modelTest);
        blockedBiomass(i) = abs(FBAsolution.f) < 1e-6;
    end

    blockedBiomassMets = biomassMets(blockedBiomass);
end




