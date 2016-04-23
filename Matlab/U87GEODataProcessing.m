
load all_datasets
load Affy_hg19_gene_mapping

% Note: the following code (commented out) was used to preprocess the
% microarray data

%[GSE22385_expression] = preprocessConsensusPlatform({'GSE22385'},{'HG-U133A_2'});
%save('all_datasets','GSE22385_expression')

%[GSE1923_expression] = preprocessConsensusPlatform({'GSE1923'},{'HG_U95Av2'});
%save('all_datasets','GSE1923_expression','-append')

%[GSE4412_A_expression] = preprocessConsensusPlatform({'GSE4412(A)'},{'HG-U133A'});
%save('all_datasets','GSE4412_A_expression','-append')

%[GSE4412_B_expression] = preprocessConsensusPlatform({'GSE4412(B)'},{'HG-U133B'});
%save('all_datasets','GSE4412_B_expression','-append')

%[GSE8692_expression] = preprocessConsensusPlatform({'GSE8692'},{'HG-U133A'});
%save('all_datasets','GSE8692_expression','-append')

% [GSE4290_expression] = preprocessConsensusPlatform({'GSE4290'},{'HG-U133_Plus_2'});
% save('all_datasets','GSE4290_expression','-append')

%%
load probes2id
[probesInt,mappingIdx,probesIdx] = intersect(updated_mapping(:,1),probes);
geneIDs = strtrim(cellstr(num2str(probe2geneid)));
updated_mapping(:,3) = repmat({''},size(updated_mapping(:,1)));
updated_mapping(mappingIdx,3) = geneIDs(probesIdx);

%%
% GSE22385
d = 1;

calls = double(GSE22385_expression.Calls);
probes = rownames(GSE22385_expression.Data);
genes = repmat({''},size(probes));

[AffyDataInt,AffyIdx,DataIdx] = intersect(updated_mapping(:,1),probes);
genes(DataIdx) = updated_mapping(AffyIdx,3);

dataStruct(d).genes = genes;
dataStruct(d).calls = calls;

% GSE1923
d = 2;

calls = double(GSE1923_expression.Calls);
probes = rownames(GSE1923_expression.Data);
genes = repmat({''},size(probes));

[AffyDataInt,AffyIdx,DataIdx] = intersect(updated_mapping(:,1),probes);
genes(DataIdx) = updated_mapping(AffyIdx,3);

dataStruct(d).genes = genes;
dataStruct(d).calls = calls;

% GSE4412_A
d = 3;

calls = double(GSE4412_A_expression.Calls);
probes = rownames(GSE4412_A_expression.Data);
genes = repmat({''},size(probes));

[AffyDataInt,AffyIdx,DataIdx] = intersect(updated_mapping(:,1),probes);
genes(DataIdx) = updated_mapping(AffyIdx,3);

dataStruct(d).genes = genes;
dataStruct(d).calls = calls;

% GSE4412_B
d = 4;

calls = double(GSE4412_B_expression.Calls);
probes = rownames(GSE4412_B_expression.Data);
genes = repmat({''},size(probes));

[AffyDataInt,AffyIdx,DataIdx] = intersect(updated_mapping(:,1),probes);
genes(DataIdx) = updated_mapping(AffyIdx,3);

dataStruct(d).genes = genes;
dataStruct(d).calls = calls;

% GSE8692
d = 5;

calls = double(GSE8692_expression.Calls);
probes = rownames(GSE8692_expression.Data);
genes = repmat({''},size(probes));

[AffyDataInt,AffyIdx,DataIdx] = intersect(updated_mapping(:,1),probes);
genes(DataIdx) = updated_mapping(AffyIdx,3);

dataStruct(d).genes = genes;
dataStruct(d).calls = calls;

% GSE4290
d = 6;

calls = double(GSE4290_expression.Calls);
probes = rownames(GSE4290_expression.Data);
genes = repmat({''},size(probes));

[AffyDataInt,AffyIdx,DataIdx] = intersect(updated_mapping(:,1),probes);
genes(DataIdx) = updated_mapping(AffyIdx,3);

dataStruct(d).genes = genes;
dataStruct(d).calls = calls;

%% Combine data

X_G = [];
G = [];

for i = 1:numel(dataStruct)
    G_i = dataStruct(i).genes;
    keepGenes = cellfun('isempty',strfind(G_i,'NaN')) & ~cellfun('isempty',G_i);
    X_i = dataStruct(i).calls(keepGenes,:); G_i = G_i(keepGenes);
    
    G_i_u = unique(G_i);
    X_i_u = zeros(numel(G_i_u),size(X_i,2));
    for j = 1:numel(G_i_u)
        G_j_idx = strmatch(G_i_u{j},G_i,'exact');
        X_j = X_i(G_j_idx,:);
        X_i_u(j,:) = max(X_j,[],1);
    end
    
    G_plus_i = union(G,G_i_u);
    X_G_plus_i = zeros(numel(G_plus_i),size(X_G,2)+size(X_i_u,2));
    
    [Gint,G_plus_i_idx,G_idx] = intersect(G_plus_i,G);
    X_G_plus_i(G_plus_i_idx,1:size(X_G,2)) = X_G(G_idx,:);
    
    [Gint,G_plus_i_idx,i_idx] = intersect(G_plus_i,G_i_u);
    X_G_plus_i(G_plus_i_idx,size(X_G,2)+(1:size(X_i_u,2))) = X_i_u(i_idx,:);
    
    G = G_plus_i;
    X_G = X_G_plus_i;
end

%% Convert to ubiquity scores
U = sum(X_G >= 2,2)/size(X_G,2);




