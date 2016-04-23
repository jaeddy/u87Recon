
% update metabolite annotation in Recon2.03

model = modelRecon2beta121114_fixed_updatedID;

modelMets = model.mets;
modelMets = regexprep(modelMets,'\[\w\]','');

annotationTable = textArrayRead('1317714181111037_add2 2.tsv',9,'\t');

[hasAnnotation,annotationIdx] = ismember(modelMets,annotationTable(:,1));

annotationFields = {'metNames';...
    'metInchiKey';...
    'metCHEBIID';...
    'metHMDB';...
    'metKeggID';...
    'metPubChemID';...
    'metLMSD';};

%% Reconcile with annotation updates from Haraldsdottir publication
for i = 1:numel(annotationFields)
    if isfield(model,annotationFields{i})
        modelField = getfield(model,annotationFields{i});
        annotationField = annotationTable(annotationIdx,i+1);
        
        mismatches = ~strcmp(modelField,annotationField);
        [i, sum(mismatches)]
        modelField(mismatches) = annotationField(mismatches);
    else
        modelField = annotationTable(annotationIdx,i+1);
    end
    model = setfield(model,annotationFields{i},modelField);
end

%% Fill in missing annotation for multi-compartment metabolites
[metTable,discrepancyMat,problemMets] = metAnalyzer(model);

updateFields = find(sum(abs(discrepancyMat))>0);
for i = 1:numel(updateFields)
    modelField = metTable(2:end,updateFields(i));
    model = setfield(model,metTable{1,updateFields(i)},modelField);
end
