
%% Read input files for core gene lists
C_H_file = 'U-87_MG_genes.txt';

fid = fopen(C_H_file);
C_H_genes = textscan(fid,'%s','delimiter','\n'); C_H_genes = C_H_genes{:};
fclose(fid);

%% GEO evidence

load U87_GEO_evidence

% Convert to ubiquity scores
U = sum(X_G >= 2,2)/size(X_G,2);

%% RNA-seq evidence
clc
load U87_fpkm
load symbols2id

X_mg = fpkm(:,1:6);
% X_mg = double(fpkm_mg >= 1);
U_mg = sum(X_mg >= 1,2)/size(X_mg,2);

% Replace gene symbols with Entrez gene IDs
[symbolsInt,mappingIdx,genesIdx] = intersect(genesymbol2geneid(:,3),genes);
geneIDs = repmat({''},size(genes));
geneIDs(genesIdx) = genesymbol2geneid(mappingIdx,2);

% Sort ubiquity scores from RNA-seq corresponding to those from GEO
[idsInt,mgIdx,geoIdx] = intersect(geneIDs,G);
U_mg_matched = zeros(size(U));

% For any genes expressed in the majority of RNA-seq samples, replace GEO
% ubiquity scores with any higher RNA-seq ubiquity scores
U_mg_matched(geoIdx) = U_mg(mgIdx);

majority_idx = U_mg_matched > 0.5;
greater_idx = U_mg_matched > U;
sum(majority_idx & greater_idx)
replace_idx = majority_idx & greater_idx;

U(replace_idx) = U_mg_matched(replace_idx);

%% Save input evidence
load Recon201_U87

save('U87mCADREInputs','C_H_genes','U','G','model')
