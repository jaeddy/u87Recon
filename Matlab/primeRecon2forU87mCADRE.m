

% prime Recon 2 for U87 mCADRE
clear
clc

addRxns = 1;
useU87biomass = 1;
resetExcs = 3;
checkBiomass = 0;

% Load model
% load HR2_CbModel_Dec2012
% load Recon2.v02.mat
% model = modelRecon2beta121114;
load HR2v03_CbModel_Jan2014

%% Add missing reactions
if addRxns
    
    % This reaction is needed to syntheisze NAD from Trp (was in Recon 1, but
    % has been removed for some reason)

%     rxnName = 'NNDPR';
%     % rxnFormula = '2 h[c] + prpp[c] + quln[c] -> co2[c] + nicrnt[c] + ppi[c]';
%     metaboliteList = {'h[c]','prpp[c]','quln[c]','co2[c]','nicrnt[c]','ppi[c]'};
%     stoichCoeffList = [-2,-1,-1,1,1,1];
%     revFlag = 0;
%     lowerBound = 0;
%     upperBound = 1000;
%     objCoeff = 0;
%     subSystem = 'NAD metabolism';
%     grRule = '(23474.1)';
%     model = addReaction(model,rxnName,metaboliteList,stoichCoeffList,revFlag,...
%         lowerBound,upperBound,objCoeff,subSystem,grRule);
% 
%     % This reaction is needed to synthesize mitochondrial NADP from NAD
% 
%     rxnName = 'NADKm';
%     % rxnFormula = 'atp[m] + nad[m] -> adp[m] + h[m] + nadp[m]';
%     metaboliteList = {'atp[m]','nad[m]','adp[m]','h[m]','nadp[m]'};
%     stoichCoeffList = [-1,-1,1,1,1];
%     revFlag = 0;
%     lowerBound = 0;
%     upperBound = 1000;
%     objCoeff = 0;
%     subSystem = 'NAD metabolism';
%     grRule = '(133686.1)'; % C5orf33
%     model = addReaction(model,rxnName,metaboliteList,stoichCoeffList,revFlag,...
%         lowerBound,upperBound,objCoeff,subSystem,grRule);
%     rxnIdx = strmatch(rxnName,model.rxns);
%     model.rxnNotes{rxnIdx} = 'NOTES:See PMID: 23212377 - human mitochondrial NAD kinase JE';

    % This transport (primed glycogenin) is required for growth when glycogen
    % is in biomass

    rxnName = 'Tyr_ggnt';
    rxnFormula = 'Tyr_ggn[e] -> Tyr_ggn[c]';
    model = addReaction(model,rxnName,rxnFormula);
    model = changeRxnBounds(model,rxnName,0,'l');

    % This transport is required for growth of U87

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

end

%% Define & add U87 biomass equation
if useU87biomass
    biomass_spec = read_biomass(model, 'U87_biomass.txt');
    biomassMets = model.mets(biomass_spec ~= 0);
    biomassCoeff = biomass_spec(biomass_spec ~= 0);
    
    % Remove components
    biomassCoeff(strmatch('dak2gpe_hs[c]',biomassMets)) = [];
    biomassMets(strmatch('dak2gpe_hs[c]',biomassMets)) = [];
    
%     biomassCoeff(strmatch('tag_hs[c]',biomassMets)) = [];
%     biomassMets(strmatch('tag_hs[c]',biomassMets)) = [];
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

%% Specify fatty acid composition
faFlux = 1000;

model = changeRxnBounds(model,'ARTFR11',faFlux,'u'); % pmt (C16:0)
model = changeRxnBounds(model,'ARTFR12',faFlux,'u'); % hd (C16:1)
model = changeRxnBounds(model,'ARTFR13',faFlux,'u'); % td (C14:0)
model = changeRxnBounds(model,'ARTFR202',faFlux,'u'); % lnlnca (C18:3 n-3)
model = changeRxnBounds(model,'ARTFR203',faFlux,'u'); % lnlncg (C18:3 n-6)
model = changeRxnBounds(model,'ARTFR204',faFlux,'u'); % strdnc (C18:4 n-3)
model = changeRxnBounds(model,'ARTFR205',faFlux,'u'); % dlnlcg (C20:3 n-6)
model = changeRxnBounds(model,'ARTFR206',faFlux,'u'); % arachd (C20:4)
model = changeRxnBounds(model,'ARTFR207',faFlux,'u'); % arach (C20:0)
model = changeRxnBounds(model,'ARTFR208',faFlux,'u'); % tmndnc (C20:5 n-3)
model = changeRxnBounds(model,'ARTFR209',faFlux,'u'); % adrn (C22:4)
model = changeRxnBounds(model,'ARTFR210',faFlux,'u'); % lnlc (C18:2 n-6)
model = changeRxnBounds(model,'ARTFR211',faFlux,'u'); % clpnd (C22:5, n-3)
model = changeRxnBounds(model,'ARTFR212',faFlux,'u'); % dcsptnl (C22:5 n-6)
model = changeRxnBounds(model,'ARTFR213',faFlux,'u'); % c226 (C22:6)
model = changeRxnBounds(model,'ARTFR31',faFlux,'u'); % st (C18:0)
model = changeRxnBounds(model,'ARTFR32',faFlux,'u'); % ode (C18:1 cis-9)
model = changeRxnBounds(model,'ARTFR33',faFlux,'u'); % vacc (C18:1 trans-11)
model = changeRxnBounds(model,'ARTFR34',faFlux,'u'); % lneld (C18:2)
model = changeRxnBounds(model,'ARTFR41',faFlux,'u'); % hd (C16:1)
model = changeRxnBounds(model,'ARTFR42',faFlux,'u'); % ode (C18:1 cis-9)
model = changeRxnBounds(model,'ARTFR43',faFlux,'u'); % vacc (C19:1 trans-11)
model = changeRxnBounds(model,'ARTFR44',faFlux,'u'); % lneld (C18:2)
model = changeRxnBounds(model,'ARTFR45',faFlux,'u'); % nrvnc (C24:1)
model = changeRxnBounds(model,'ARTFR46',faFlux,'u'); % od2 (C18:1 trans-9)
model = changeRxnBounds(model,'ARTFR51',faFlux,'u'); % lgnc (C24:0)
model = changeRxnBounds(model,'ARTFR52',faFlux,'u'); % hexc (C26:0)
model = changeRxnBounds(model,'ARTFR53',faFlux,'u'); % eicostet (C20:4 n-3)
model = changeRxnBounds(model,'ARTFR54',faFlux,'u'); % tetpent6 (C24:5 n-6)
model = changeRxnBounds(model,'ARTFR55',faFlux,'u'); % tetpent3 (C24:5 n-3)
model = changeRxnBounds(model,'ARTFR56',faFlux,'u'); % tettet6 (C24:4 n-6)
model = changeRxnBounds(model,'ARTFR57',faFlux,'u'); % tethex3 (C24:6 n-3)
model = changeRxnBounds(model,'ARTFR61',faFlux,'u'); % hdd2 (C16:1 trans-2)

%% Enable R group synthesis
rFlux = 1000;

model = changeRxnBounds(model,'RTOT_2',rFlux,'u'); 
model = changeRxnBounds(model,'RTOT_3',rFlux,'u'); 
model = changeRxnBounds(model,'RTOT1',rFlux,'u'); % C16 groups
model = changeRxnBounds(model,'RTOT2',rFlux,'u'); % essental + derivatives
model = changeRxnBounds(model,'RTOT3',rFlux,'u'); % non-essential
model = changeRxnBounds(model,'RTOT4',rFlux,'u'); % monounsaturated
model = changeRxnBounds(model,'RTOT5',rFlux,'u'); % other (C26, C24)
model = changeRxnBounds(model,'RTOT6',rFlux,'u'); % oxidase products

%% Add cytosolic ATP maintenance reaction
% ~~ Recon2 currently has mitochondrial ATP maintenance, so skip for now
% model = addReaction(model,'DM_atp_c_','atp[c] + h2o[c] -> adp[c] + h[c] + pi[c]');
% model = changeRxnBounds(model,'DM_atp_c_',0,'l');

if resetExcs > 0
    %% Select model inputs
    excInd = findExcRxns(model);
    excRxns = model.rxns(excInd);
    excLBs = model.lb(excInd);
    excUBs = model.ub(excInd);

    excRxns_small = excRxns(rem(excLBs,10)~=0 & excLBs~=0);
    excRxns_large = excRxns(rem(excLBs,10)==0 & excLBs~=0);
    
    excFormulas = printRxnFormula(model,excRxns,0);
    excMets = strtrim(regexprep(excFormulas,'(->|<=>)',''));
    excMetFormulas = model.metFormulas(directIntersect(model.mets,excMets));
    hasX = ~cellfun('isempty',strfind(excMetFormulas,'X'));
    excRxns_hasX = excRxns(hasX);
    excLBs_hasX = excLBs(hasX);
    
    %% Clear model inputs 
    switch resetExcs
        case 1
            model = changeRxnBounds(model,excRxns_large,0,'l');
        case 2
            model = changeRxnBounds(model,excRxns_small,0,'l');
        case 3
            model = changeRxnBounds(model,excRxns,0,'l');
    end
        
    % Reset bounds of exchanges for X-containing metabolites
    excLBs_hasX = -1000*ones(size(excLBs_hasX));
    model = changeRxnBounds(model,excRxns_hasX,excLBs_hasX,'l');
%     model = changeRxnBounds(model,'EX_Tyr_ggn(e)',-1000,'l');
    
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
%     aaList = {'EX_gln_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
%         'EX_tyr_L(e)'};

    % DMEM list (less restrictive)
%     aaList = {'EX_gln_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
%         'EX_tyr_L(e)','EX_arg_L(e)','EX_thr_L(e)','EX_val_L(e)'};

    % DMEM list (full)
%     aaList = {'EX_arg_L(e)','EX_cys_L(e)','EX_gln_L(e)','EX_gly(e)',...
%         'EX_his_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
%         'EX_met_L(e)','EX_phe_L(e)','EX_ser_L(e)','EX_thr_L(e)',...
%         'EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)'};

    % Essential AAs
    aaList = {'EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
        'EX_met_L(e)','EX_phe_L(e)','EX_thr_L(e)','EX_trp_L(e)',...
        'EX_val_L(e)','EX_his_L(e)'};
    model = changeRxnBounds(model,aaList,-1000,'l');

    
    %% Add fatty acids
    % fatty acids:
    % - linoleic acid (lnlc)
    % - alpha-linolenic acid (lnlnca)
    
    % Essential FAs
    faList = {'EX_lnlc(e)','EX_lnlnca(e)'};
    model = changeRxnBounds(model,faList,-1000,'l');
    
    %% Add vitamins
    % DMEM vitamins:
    % - pantohenate / vitamin B5 (pnto_R)
    % - riboflavin (ribflv)
    
    % Essential vitamins
    vitList = {'EX_pnto_R(e)','EX_ribflv(e)'};
    model = changeRxnBounds(model,vitList,-1000,'l');
    
    %% Add other components
    % B27 supplement:
    % - reduced glutathione (gthrd)
    % - corticosterone (crtstrn)
    % - D-galactose (gal)
    % - ethanolamine (etha)
    % - L-carnitine (crn)
    % - progesterone (prgstrn)
    % - putrescine (ptrc)
    % - triiodothyronine (triodthy)
    
%     otherList = {'EX_gthrd(e)'};
%     model = changeRxnBounds(model,otherList,-1000,'l');
    
    %% Adjust substrate uptake
    % carbon substrates:
    % - glucose (glc)

    substrateList = {'EX_glc(e)'};
%     substrateList = [substrateList,...
%         {'EX_tag_hs(e)'}];
    model = changeRxnBounds(model,substrateList,-10,'l');    
end

%% Optimize model
% clc
changeCobraSolver('glpk');
FBAsolution = optimizeCbModel(model);
printFluxVector(model,FBAsolution.x,1,1,-1,[],[],0)

%% Check biomass components
if checkBiomass
    clc
    [modelTest,rxnNames] = addDemandReaction(model,biomassMets);
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

%% Trace precursors

% broken = {'pail_hs[c]';...
%     'pchol_hs[c]';...
%     'pe_hs[c]';...
%     'pglyc_hs[c]';...
%     'clpn_hs[c]';...
%     'ps_hs[c]';...
%     'sphmyln_hs[c]';};
% clc
% 
% % met = broken{7};
% % met = 'R2coa_hs[c]';
% met = 'nadp[m]';
% rxns = findRxnsFromMets(model,met);
% printRxnFormula(model,rxns);
% 
% % Check biomass precursors
% % met = 'pail3p_hs[c]';
% [modelTest,precursorRxn] = addDemandReaction(model,met);
% % [modelTest] = addDemandReaction(modelTest,'fad[m]');
% % [modelTest] = addDemandReaction(modelTest,'nadp[m]');
% 
% % modelTest = changeRxnBounds(modelTest,'DM_pnto_R',-1000,'l');
% % modelTest = changeRxnBounds(modelTest,'EX_ribflv(e)',-1000,'l');
% 
% % (add false exchanges)
% % exMet = 'R2coa_hs[c]';
% % modelTest = addExchangeRxn(modelTest,{exMet},{['EX_',exMet]},-1000,1000);
% % exMet = 'nadph[m]';
% % modelTest = addExchangeRxn(modelTest,{exMet},{['EX_',exMet]},-1000,1000);
% % exMet = 'fadh2[m]';
% % modelTest = addExchangeRxn(modelTest,{exMet},{['EX_',exMet]},-1000,1000);
% % exMet = 'lnlncacoa[c]';
% % modelTest = addExchangeRxn(modelTest,{exMet},{['EX_',exMet]},-1000,1000);
% 
% 
% modelTest = changeObjective(modelTest,precursorRxn,1);
% % modelTest = changeObjective(modelTest,'EX_pnto_R(e)',1);
% FBAsolution = optimizeCbModel(modelTest);
% result = {'Success','Fail'};
% display([met,'? ',result{2-double(FBAsolution.f > 1e-6)},': ',...
%     num2str(FBAsolution.f)])


