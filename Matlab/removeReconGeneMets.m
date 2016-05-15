
% remove gene species from mets

lastMet = 5056;

model2 = removeMetabolites(model,model.mets(lastMet+1:end));

model2.metCharge(lastMet+1:end) = [];
model2.metCHEBIID(lastMet+1:end) = [];
model2.metKeggID(lastMet+1:end) = [];
model2.metPubChemID(lastMet+1:end) = [];
model2.metHMDB(lastMet+1:end) = [];
model2.metInchiString(lastMet+1:end) = [];
model2.metEHMNID(lastMet+1:end) = [];
model2.metHepatoNetID(lastMet+1:end) = [];

% model = model2;
% save HR2_CbModel_Dec2012 model