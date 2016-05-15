function model = convertSBMLToCobra_Recon2(modelSBML,defaultBound,compSymbolList,compNameList)
%convertSBMLToCobra Convert SBML format model (created using SBML Toolbox)
%to Cobra format
%
% model = convertSBMLToCobra(modelSBML,defaultBound)
%
%INPUTS
% modelSBML         SBML model structure
%
%OPTIONAL INPUTS
% defaultBound      Maximum bound for model (Default = 1000)
% compSymbolList    List of compartment symbols
% compNameList      List of compartment names corresponding to compSymbolList
%
%OUTPUT
% model             COBRA model structure
% Markus Herrgard 1/25/08
%
% Ines Thiele 01/27/2010 - I added new field to be read-in from SBML file
% if provided in file (e.g., references, comments, metabolite IDs, etc.)
%
% Richard Que 02/08/10 - Properly format reaction and metabolite fields
%                        from SBML.
%

if (nargin < 2)
    defaultBound = 1000;
end

if nargin < 3
    compSymbolList = {};
    compNameList = {};
end

%%
nMetsTmp = length(modelSBML.species);
nRxns = length(modelSBML.reaction);

%% Construct initial metabolite list
formulaCount = 0;
speciesList = {};
chargeList = [];
metFormulas = {};
haveFormulasFlag = false;
tmpSpecies = [];
metInChi = {};
metNeutralFormulas = {};
metOrigins = {};
ehmnAbbrevs = {};
hepatonetAbbrevs = {};

for i = 1:nMetsTmp
    % Ignore exchange species
    if (isempty(regexp(modelSBML.species(i).id,'_b$','ONCE')));
        tmpSpecies = [ tmpSpecies  modelSBML.species(i)];
        speciesList{end+1} = modelSBML.species(i).id;
        if isfield(modelSBML.species(i),'charge')
            tmpCharge = double(modelSBML.species(i).charge);
        else
            tmpCharge = 0;
        end
        notesField = modelSBML.species(i).notes;
        % Get formula (and other stuff) if in notes field
        if (~isempty(notesField))
            [tmp,tmp,tmp,tmp,formula,tmp,tmp,tmp,tmp,charge,inchi,nFormula,origin,ehmnAbbrev,hepatoAbbrev] = parseSBMLNotesField(notesField);
            if numel(charge)
                tmpCharge = charge;
            end
            metFormulas{end+1} = formula;
            metInChi{end+1} = inchi;
            metNeutralFormulas{end+1} = nFormula;
            metOrigins{end+1} = origin;
            ehmnAbbrevs{end+1} = ehmnAbbrev;
            hepatonetAbbrevs{end+1} = hepatoAbbrev;
            formulaCount = formulaCount + 1;
            haveFormulasFlag = true;
        else
            metFormulas{end+1} = '';
            metInChi{end+1} = '';
            metNeutralFormulas{end+1} = '';
            metOrigins{end+1} = '';
            ehmnAbbrevs{end+1} = '';
            hepatonetAbbrevs{end+1} = '';
        end
        chargeList= [chargeList tmpCharge];
    end
end

nMets = length(speciesList);

%% Construct stoichiometric matrix and reaction list
S = sparse(nMets,nRxns);
rev = zeros(nRxns,1);
lb = zeros(nRxns,1);
ub = zeros(nRxns,1);
c = zeros(nRxns,1);
rxns = cell(nRxns,1);
rules = cell(nRxns,1);
genes = cell(nRxns,1);
allGenes = {};
h = waitbar(0,'Reading SBML file ...');
hasNotesField = false;
for i = 1:nRxns
    if mod(i,10) == 0
        waitbar(i/nRxns,h);
    end
    % Read the gpra from the notes field
    notesField = modelSBML.reaction(i).notes;
    if (~isempty(notesField))
        [geneList,rule,subSystem,grRule,formula,confidenceScore,citation,comment,ecNumber,tmp,tmp,tmp,origin] = parseSBMLNotesField(notesField);
        subSystems{i} = subSystem;
        genes{i} = geneList;
        allGenes = [allGenes geneList];
        rules{i} = rule;
        grRules{i} = grRule;
        hasNotesField = true;
        confidenceScores{i}= confidenceScore;
        citations{i} = citation;
        comments{i} = comment;
        ecNumbers{i} = ecNumber;
        rxnOrigins{i} = origin;
    end
    rev(i) = modelSBML.reaction(i).reversible;
    rxnNameTmp = regexprep(modelSBML.reaction(i).name,'^R_','');
    rxnNames{i} = regexprep(rxnNameTmp,'_+',' ');
    rxnsTmp = regexprep(modelSBML.reaction(i).id,'^R_','');
    rxns{i} = cleanUpFormatting(rxnsTmp);
    % Construct S-matrix
    reactantStruct = modelSBML.reaction(i).reactant;
    for j = 1:length(reactantStruct)
        speciesID = find(strcmp(reactantStruct(j).species,speciesList));
        if (~isempty(speciesID))
            stoichCoeff = reactantStruct(j).stoichiometry;
            S(speciesID,i) = -stoichCoeff;
        end
    end
    productStruct = modelSBML.reaction(i).product;
    for j = 1:length(productStruct)
        speciesID = find(strcmp(productStruct(j).species,speciesList));
        if (~isempty(speciesID))
            stoichCoeff = productStruct(j).stoichiometry;
            S(speciesID,i) = stoichCoeff;
        end
    end
    if isfield(modelSBML.reaction(i).kineticLaw,'parameter')
        parameters = modelSBML.reaction(i).kineticLaw.parameter;
    else
        parameters =[];
    end
    if (~isempty(parameters))
        for j = 1:length(parameters)
            paramStruct = parameters(j);
            switch paramStruct.id
                case 'LOWER_BOUND'
                    lb(i) = paramStruct.value;
                    if (lb(i) < -defaultBound)
                        lb(i) = -defaultBound;
                    end
                case 'UPPER_BOUND'
                    ub(i) = paramStruct.value;
                    if (ub(i) > defaultBound)
                        ub(i) = defaultBound;
                    end
                case 'OBJECTIVE_COEFFICIENT'
                    c(i) = paramStruct.value;
            end
        end
    else
        ub(i) = defaultBound;
        if (rev(i) == 1)
            lb(i) = -defaultBound;
        else
            lb(i) = 0;
        end
    end
end
%close the waitbar if this is matlab
if ( regexp( version, 'R20') )
    close(h);
end
allGenes = unique(allGenes);

%% Construct gene to rxn mapping
if (hasNotesField)
    
    rxnGeneMat = sparse(nRxns,length(allGenes));
    h = waitbar(0,'Constructing GPR mapping ...');
    for i = 1:nRxns
        if mod(i,10) == 0
            waitbar(i/nRxns,h);
        end
        if iscell(genes{i})
            [tmp,geneInd] = ismember(genes{i},allGenes);
        else
            [tmp,geneInd] = ismember(num2cell(genes{i}),allGenes);
        end
        
        rxnGeneMat(i,geneInd) = 1;
        for j = 1:length(geneInd)
            %rules{i} = strrep(rules{i},['x(' num2str(j) ')'],['x(' num2str(geneInd(j)) ')']);
            rules{i} = strrep(rules{i},['x(' num2str(j) ')'],['x(' num2str(geneInd(j)) '_TMP_)']);
        end
        rules{i} = strrep(rules{i},'_TMP_','');
    end
    %close the waitbar if this is matlab
    if ( regexp( version, 'R20') )
        close(h);
    end
    
end

%% Construct metabolite list
mets = cell(nMets,1);
compartmentList = cell(length(modelSBML.compartment),1);
if isempty(compSymbolList), useCompList = true; else useCompList = false; end
for i=1:length(modelSBML.compartment)
    compartmentList{i} = modelSBML.compartment(i).id;
end

h = waitbar(0,'Constructing metabolite lists ...');
hasAnnotationField = 0;
for i = 1:nMets
    if mod(i,10) == 0
        waitbar(i/nMets,h);
    end
    % Parse metabolite id's
    % Get rid of the M_ in the beginning of metabolite id's
    metID = regexprep(speciesList{i},'^M_','');
    metID = regexprep(metID,'^_','');
    % Find compartment id
    tmpCell = {};
    if useCompList
        for j=1:length(compartmentList)
            tmpCell = regexp(metID,['_(' compartmentList{j} ')$'],'tokens');
            if ~isempty(tmpCell), break; end
        end
        if isempty(tmpCell), useCompList = false; end
    elseif ~isempty(compSymbolList)
        for j = 1: length(compSymbolList)
            tmpCell = regexp(metID,['_(' compSymbolList{j} ')$'],'tokens');
            if ~isempty(tmpCell), break; end
        end
    end
    if isempty(tmpCell), tmpCell = regexp(metID,'_(.)$','tokens'); end
    if ~isempty(tmpCell)
        compID = tmpCell{1};
        metTmp = [regexprep(metID,['_' compID{1} '$'],'') '[' compID{1} ']'];
%         metTmp = regexprep(metTmp,'_+','-');
%         metTmp = regexprep(metTmp,'-DASH-','-');%%%
    else
        metTmp = metID;
    end
    %Clean up met ID
    mets{i} = cleanUpFormatting(metTmp);
    % Parse metabolite names
    % Clean up some of the weird stuff in the sbml files
    metNamesTmp = regexprep(tmpSpecies(i).name,'^M_','');
    metNamesTmp = cleanUpFormatting(metNamesTmp);
    metNamesTmp = regexprep(metNamesTmp,'^_','');
%     metNamesTmp = strrep(metNamesTmp,'_','-');
    metNamesTmp = regexprep(metNamesTmp,'-+','-');
    metNamesTmp = regexprep(metNamesTmp,'-$','');
    metNamesAlt{i} = metNamesTmp;
    % Separate formulas from names
    %[tmp,tmp,tmp,tmp,tokens] = regexp(metNamesTmp,'(.*)-((([A(Ag)(As)C(Ca)(Cd)(Cl)(Co)(Cu)F(Fe)H(Hg)IKLM(Mg)(Mn)N(Na)(Ni)OPRS(Se)UWXY(Zn)]?)(\d*)))*$');
    if (~haveFormulasFlag)
        [tmp,tmp,tmp,tmp,tokens] = regexp(metNamesTmp,'(.*)_((((A|Ag|As|C|Ca|Cd|Cl|Co|Cu|F|Fe|H|Hg|I|K|L|M|Mg|Mn|Mo|N|Na|Ni|O|P|R|S|Se|U|W|X|Y|Zn)?)(\d*)))*$');
        if (isempty(tokens))
            if length(metFormulas)<i||(metFormulas{i}=='')
                metFormulas{i} = '';
            end
            metNames{i} = metNamesTmp;
        else
            formulaCount = formulaCount + 1;
            metFormulas{i} = tokens{1}{2};
            metNames{i} = tokens{1}{1};
        end
    else
        metNames{i} = metNamesTmp;
    end
    if isfield(modelSBML.species(i),'annotation')
        hasAnnotationField = 1;
        [metCHEBI,metKEGG,metPubChem,metInChI,metHMDB] = parseSBMLAnnotationField(modelSBML.species(i).annotation);
        metCHEBIID{i} = metCHEBI;
        metKEGGID{i} = metKEGG;
        metPubChemID{i} = metPubChem;
        metInChIString{i} = metInChI;
        metHMDBID{i} = metHMDB;
    end
end
if ( regexp( version, 'R20') )
    close(h);
end

%% Collect everything into a structure
model.rxns = rxns;
model.mets = mets;
model.S = S;
model.rev = rev;
model.lb = lb;
model.ub = ub;
model.c = c;
model.metCharge = transpose(chargeList);
if (hasNotesField)
    model.rules = rules;
    model.genes = columnVector(allGenes);
    model.rxnGeneMat = rxnGeneMat;
    model.grRules = columnVector(grRules);
    model.subSystems = columnVector(subSystems);
    model.rxnOrigins = columnVector(rxnOrigins);
    model.confidenceScores = columnVector(confidenceScores);
    model.rxnReferences = columnVector(citations);
    model.rxnECNumbers = columnVector(ecNumbers);
    model.rxnNotes = columnVector(comments);
end
model.rxnNames = columnVector(rxnNames);

model.metNames = columnVector(metNames);
model.metFormulas = columnVector(metFormulas);
model.metNeutralFormulas = columnVector(metNeutralFormulas);
    
if (hasAnnotationField)
    model.metOrigins = columnVector(metOrigins);
    model.metHMDBID = columnVector(metHMDBID);
    model.metChEBIID = columnVector(metCHEBIID);
    model.metKEGGID = columnVector(metKEGGID);
    model.metPubChemID = columnVector(metPubChemID);
    model.metInChIString = columnVector(metInChIString);
    model.metEHMNAbbrevs = columnVector(ehmnAbbrevs);
    model.metHepatonetAbrrevs = columnVector(hepatonetAbbrevs);
end

%%
function [genes,rule,subSystem,grRule,formula,confidenceScore,citation,comment,ecNumber,charge,inchi,nFormula,origin,ehmnAbbrev,hepatoAbbrev] = parseSBMLNotesField(notesField)
%parseSBMLNotesField Parse the notes field of an SBML file to extract
%gene-rxn associations
%
% [genes,rule] = parseSBMLNotesField(notesField)
%
% Markus Herrgard 8/7/06
% Ines Thiele 1/27/10 Added new fields
% Handle different notes fields
if isempty(regexp(notesField,'html:p', 'once'))
    tag = 'p';
else
    tag = 'html:p';
end

subSystem = '';
grRule = '';
genes = {};
rule = '';
formula = '';
confidenceScore = '';
citation = '';
ecNumber = '';
comment = '';
charge = [];
Comment = 0;
inchi = '';
nFormula = '';
origin = '';
ehmnAbbrev = '';
hepatoAbbrev = '';

[tmp,fieldList] = regexp(notesField,['<' tag '>.*?</' tag '>'],'tokens','match');

for i = 1:length(fieldList)
    fieldTmp = regexp(fieldList{i},['<' tag '>(.*)</' tag '>'],'tokens');
    fieldStr = fieldTmp{1}{1};
    if (regexp(fieldStr,'GENE_ASSOCIATION'))
        gprStr = regexprep(strrep(fieldStr,'GENE_ASSOCIATION:',''),'^(\s)+','');
        grRule = gprStr;
        [genes,rule] = parseBoolean(gprStr);
    elseif (regexp(fieldStr,'GENE ASSOCIATION'))
        gprStr = regexprep(strrep(fieldStr,'GENE ASSOCIATION:',''),'^(\s)+','');
        grRule = gprStr;
        [genes,rule] = parseBoolean(gprStr);
    elseif (regexp(fieldStr,'SUBSYSTEM'))
        subSystem = regexprep(strrep(fieldStr,'SUBSYSTEM:',''),'^(\s)+','');
        subSystem = strrep(subSystem,'S_','');
        subSystem = regexprep(subSystem,'_+',' ');
        if (isempty(subSystem))
            subSystem = 'Exchange';
        end
    elseif (regexp(fieldStr,'EC Number'))
        ecNumber = regexprep(strrep(fieldStr,'EC Number:',''),'^(\s)+','');
    elseif (regexp(fieldStr,'NEUTRAL_FORMULA'))
        nFormula = regexprep(strrep(fieldStr,'NEUTRAL_FORMULA:',''),'^(\s)+','');
    elseif (regexp(fieldStr,'FORMULA'))
        formula = regexprep(strrep(fieldStr,'FORMULA:',''),'^(\s)+','');
    elseif (regexp(fieldStr,'CHARGE'))
        charge = str2num(regexprep(strrep(fieldStr,'CHARGE:',''),'^(\s)+',''));
    elseif (regexp(fieldStr,'AUTHORS'))
        if isempty(citation)
            citation = strcat(regexprep(strrep(fieldStr,'AUTHORS:',''),'^(\s)+',''));
        else
            citation = strcat(citation,';',regexprep(strrep(fieldStr,'AUTHORS:',''),'^(\s)+',''));
        end
    elseif (regexp(fieldStr,'INCHI'))
        inchi = regexprep(strrep(fieldStr,'INCHI:',''),'^(\s)+','');
    elseif (regexp(fieldStr,'ORIGIN'))
        origin = regexprep(strrep(fieldStr,'ORIGIN:',''),'^(\s)+','');
    elseif (regexp(fieldStr,'EHMN_ABBREVIATION'))
        ehmnAbbrev = regexprep(strrep(fieldStr,'EHMN_ABBREVIATION:',''),'^(\s)+','');
    elseif (regexp(fieldStr,'HEPATONET_1.0_ABBREVIATION'))
        hepatoAbbrev = regexprep(strrep(fieldStr,'HEPATONET_1.0_ABBREVIATION:',''),'^(\s)+','');
    elseif Comment == 1 && isempty(regexp(fieldStr,'genes:', 'once'))
        Comment = 0;
        comment = fieldStr;
    elseif (regexp(fieldStr,'Confidence'))
        confidenceScore = regexprep(strrep(fieldStr,'Confidence Level:',''),'^(\s)+','');
        Comment = 1;
    
    end
end
%%
function [metCHEBI,metKEGG,metPubChem,metInChI,metHMDB] = parseSBMLAnnotationField(annotationField)
%parseSBMLAnnotationField Parse the annotation field of an SBML file to extract
%metabolite information associations
%
% [genes,rule] = parseSBMLAnnotationField(annotationField)
%
% Ines Thiele 1/27/10 Added new fields
% Handle different notes fields


metCHEBI = '';
metKEGG = '';
metPubChem = '';
metInChI = '';
metHMDB = '';

[tmp,fieldList] = regexp(annotationField,'<rdf:li rdf:resource="urn:miriam:(\w+).*?"/>','tokens','match');

%fieldTmp = regexp(fieldList{i},['<' tag '>(.*)</' tag '>'],'tokens');
for i = 1:length(fieldList)
    fieldTmp = regexp(fieldList{i},['<rdf:li rdf:resource="urn:miriam:(.*)"/>'],'tokens');
    fieldStr = fieldTmp{1}{1};
    if (regexp(fieldStr,'obo.chebi'))
        metCHEBI = strrep(fieldStr,'obo.chebi:CHEBI%','');
    elseif (regexp(fieldStr,'kegg.compound'))
        metKEGG = strrep(fieldStr,'kegg.compound:','');
    elseif (regexp(fieldStr,'pubchem.substance'))
        metPubChem = strrep(fieldStr,'pubchem.substance:','');
    elseif (regexp(fieldStr,'hmdb'))
        metHMDB = strrep(fieldStr,'hmdb:','');
    end
end
% get InChI string
fieldTmp = regexp(annotationField,'<in:inchi xmlns:in="http://biomodels.net/inchi" metaid="(.*?)">(.*?)</in:inchi>','tokens');
if ~isempty(fieldTmp)
fieldStr = fieldTmp{1}{2};
if (regexp(fieldStr,'InChI'))
    metInChI = strrep(fieldStr,'InChI=','');
end
end

%% Cleanup Formatting
function str = cleanUpFormatting(str)
str = strrep(str,'-DASH-','-');
str = strrep(str,'_DASH_','-');
str = strrep(str,'_FSLASH_','/');
str = strrep(str,'_BSLASH_','\');
str = strrep(str,'_LPAREN_','(');
str = strrep(str,'_LSQBKT_','[');
str = strrep(str,'_RSQBKT_',']');
str = strrep(str,'_RPAREN_',')');
str = strrep(str,'_COMMA_',',');
str = strrep(str,'_PERIOD_','.');
str = strrep(str,'_APOS_','''');
str = regexprep(str,'_e_$','(e)');
str = regexprep(str,'_e$','(e)');
str = strrep(str,'&amp;','&');
str = strrep(str,'&lt;','<');
str = strrep(str,'&gt;','>');
str = strrep(str,'&quot;','"');