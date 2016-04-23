function model = formulateDMEM(model,DMEMopt)
resetExcs = 3;

%%
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

switch DMEMopt
    case 1
    % DMEM list (restrictive)
    aaList = {'EX_gln_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
        'EX_tyr_L(e)'};
    case 2
    % DMEM list (less restrictive)
    aaList = {'EX_gln_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
        'EX_tyr_L(e)','EX_arg_L(e)','EX_thr_L(e)','EX_val_L(e)'};
    case 3
    % DMEM list (full)
    aaList = {'EX_arg_L(e)','EX_cys_L(e)','EX_gln_L(e)','EX_gly(e)',...
        'EX_his_L(e)','EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
        'EX_met_L(e)','EX_phe_L(e)','EX_ser_L(e)','EX_thr_L(e)',...
        'EX_trp_L(e)','EX_tyr_L(e)','EX_val_L(e)'};
end
%     % Essential AAs
%     aaList = {'EX_ile_L(e)','EX_leu_L(e)','EX_lys_L(e)',...
%         'EX_met_L(e)','EX_phe_L(e)','EX_thr_L(e)','EX_trp_L(e)',...
%         'EX_val_L(e)','EX_his_L(e)'};
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
    model = changeRxnBounds(model,substrateList,-1000,'l');    
end
