
clear 
clc

% clear files
dirInfo = dir;
[files{1:length(dirInfo)}] = deal(dirInfo.name); files = files(:);
dataFiles = files(~cellfun('isempty',regexp(files,'[0-9]\.txt','match')));
mapFiles = files(~cellfun('isempty',regexp(files,'map.txt','match')));

sampleSize = [...
    3;... % GSE1923
    2;... % GSE21217
    3;... % GSE22385
    3;... % GSE32100
    3;... % GSE35208
    1;... % GSE45161
%     2;... % GSE4536 % need to remove this dataset (xenograft)
    1];   % GSE9171

for d = 1:numel(dataFiles)
% First, process data from individual data sets

% Import GEO data from R-generated txt files
numSamples = str2double(cell2mat(...
    regexp(dataFiles{d},'(?<=_)[0-9]','match')));
dataset_PA = textArrayRead(dataFiles{d},numSamples+1,'\t');
dataset_PA(1,:) = [dataset_PA(1,end),dataset_PA(1,1:end-1)];
dataset_PA = regexprep(dataset_PA,'"','');

dataset_map = textArrayRead(mapFiles{d},2,'\t');
dataset_map = regexprep(dataset_map,'"','');

% Map probe IDs to gene IDs
[mapProbes,firstIdx] = unique(dataset_map(:,1),'first');
[~,lastIdx] = unique(dataset_map(:,1),'last');
multiples = mapProbes(lastIdx-firstIdx > 0);

% Remove probes mapping to multiple genes
dataset_map_filtered = ...
    dataset_map(~ismember(dataset_map(:,1),multiples),:);

[~,gseIdx,mapIdx] = intersect(dataset_PA,dataset_map_filtered);
dataset_PA_filtered = [dataset_map_filtered([1;mapIdx],2),...
    dataset_PA([1;gseIdx],2:end)];

dataset_PA_filtered(2:end,:) = sortrows(dataset_PA_filtered(2:end,:),1);

% Compress presence/absence calls for genes with multiple probes
[geneIDs,firstIdx] = unique(dataset_PA_filtered(2:end,1),'first');
[~,lastIdx] = unique(dataset_PA_filtered(2:end,1),'last');

% For each sample, if gene is present according to at least one probe,
% designate the gene as present
calls = repmat({'A'},numel(geneIDs),size(dataset_PA_filtered,2)-1);
for i = 1:numel(geneIDs)
    if lastIdx(i)-firstIdx(i) > 0
        range_i = (firstIdx(i):lastIdx(i))+1;
        P_i = ~cellfun('isempty',...
            regexp(dataset_PA_filtered(range_i,2:end),'P'));
      calls(i,sum(P_i)>0) = {'P'};
    else
        calls(i,:) = dataset_PA_filtered(firstIdx(i)+1,2:end);
    end
end

dataset_PA_d = [dataset_PA_filtered(1,:);[geneIDs,calls]];

% Second, combine all datasets based on common gene IDs
if d == 1
    combined_PA = dataset_PA_d;
else
    % Find union list of gene IDs
    G = union(dataset_PA_d(2:end,1),combined_PA(2:end,1));
    numGenes = numel(G);
    sizeDataset = size(dataset_PA_d(:,2:end),2);
    sizeCombined = size(combined_PA(:,2:end),2);
    
    % Initialize updated matrix
    updated_PA = repmat({''},numGenes,sizeDataset+sizeCombined);
    
    % Add original data from combined_PA to updated matrix
    [~,updatedIdx,combinedIdx] = intersect(G,combined_PA(:,1));
    updated_PA(updatedIdx,1:sizeCombined) = ...
        combined_PA(combinedIdx,2:end);
    
    % Add new dataset to updated matrix
    [~,updatedIdx,datasetIdx] = intersect(G,dataset_PA_d(:,1));
    updated_PA(updatedIdx,sizeCombined+1:sizeCombined+sizeDataset) = ...
        dataset_PA_d(datasetIdx,2:end);
    
    % Re-add row names and column headers
    combined_PA = [[combined_PA(1,:),dataset_PA_d(1,2:end)];...
        [G,updated_PA]];
end

end

%% Convert presence/absence calls to ubiquity scores

X_G = zeros(size(combined_PA(2:end,2:end)));
X_G(~cellfun('isempty',strfind(combined_PA(2:end,2:end),'P'))) = 1;
X_G(cellfun('isempty',combined_PA(2:end,2:end))) = 0.5;

U = sum(X_G,2)./size(X_G,2);


