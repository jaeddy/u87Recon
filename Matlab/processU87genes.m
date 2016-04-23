U87genes = findGenesFromRxns(PM,PM.rxns);
numgenes = sum(cellfun('length',U87genes));
geneList = repmat({''},numgenes,1);

count = 0;
for i = 1:numel(U87genes)
    genes_i = U87genes{i};
    geneList(count+1:count+length(genes_i)) = genes_i;
    count = count+length(genes_i);
end

U87genes_unique = unique(geneList);

%%
U87genes_unique = regexprep(U87genes_unique,'\.[0-9]*','');

textArrayWrite('U87genes.txt',U87genes_unique,'\t');

U87symbols = textArrayRead('U87symbols.txt',1,'\t');
U87symbols_unique = unique(U87symbols);

[geneRxnStruct,geneRxns] = findRxnsFromGenes(PM,U87genes_unique,1,1);

U87info = textArrayRead('U87info.txt',15,'\t');
colHeaders = {'tax_id','GeneID','Symbol','LocusTag','Synonyms','dbXrefs','chromosome','map location','description','type of gene','Symbol from nomenclature authority','Full name from nomenclature authority','Nomenclature status','Other designations','Modification date'};

[U87symbols_unique,u_idx] = unique(U87symbols);
U87info_unique = U87info(u_idx,:);
U87info_unique(:,1) = [];
U87info_unique(sum(cellfun('isempty',U87info_unique),2)==14,:) = [];

U87genes_hpa = textArrayRead('U-87_MG_genes.txt',1,'\t');
U87genes_hpa(~cellfun('isempty',regexp(lower(U87genes_hpa),'[a-z]')));

U87genes_hpa_strong = U87genes_hpa(1:1961);
U87genes_hpa_moderate = U87genes_hpa(1962:end);
U87genes_hpa_strong(~cellfun('isempty',regexp(lower(U87genes_hpa_strong),'[a-z]'))) = [];
U87genes_hpa_moderate(~cellfun('isempty',regexp(lower(U87genes_hpa_moderate),'[a-z]'))) = [];

U87info_unique = [repmat({''},1021,1),U87info_unique];
is_strong = ismember(U87info_unique(:,2),U87genes_hpa_strong);
is_moderate = ismember(U87info_unique(:,2),U87genes_hpa_moderate);
U87info_unique(is_moderate,1) = {'*'};
U87info_unique(is_strong,1) = {'**'};