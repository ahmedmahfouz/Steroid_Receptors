%%% 15 Feb. 2013
%%% calculate the correlation between two genes within a
%%% specific structure (defined by the ABA acronym)

% define the files directory
filesDirectory = 'files/';
matObj = matfile([filesDirectory 'fullExpressionMatrix.mat']);

% load gene names
load([filesDirectory 'allGenes.mat']);
% load experiments' numbers
load([filesDirectory 'allExpNumbers.mat']);
% load voxels annotations
load([filesDirectory 'voxelAnnotation.mat']);

geneA = 'Tph2';
geneB = 'Maob';
strOfInterest = 'MB';

geneA_index = find(strcmp(allGenes, geneA) == 1);
geneA_experimentNo = allExpNumbers(geneA_index);
geneB_index = find(ismember(allGenes, geneB) == 1);
geneB_experimentNo = allExpNumbers(geneB_index);
indx = [geneA_index geneB_index];
for i = 1 : length(indx)
    if i == 1
        miniExpMat = matObj.fullExpressionMatrix(:,indx(i));
    else
        miniExpMat = [miniExpMat matObj.fullExpressionMatrix(:,indx(i))];
    end
end

inclVoxAnnotations = indicateSubStr(strOfInterest, 'files/');
inclVoxelsInd = find(ismember(voxelAnnotation, inclVoxAnnotations) == 1);

expMat_logTransformed = miniExpMat;

dataMat = expMat_logTransformed(inclVoxelsInd, :);
c = corr(dataMat, 'type', 'Spearman');




