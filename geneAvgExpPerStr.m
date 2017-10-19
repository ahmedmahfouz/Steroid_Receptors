%%% 25 Sep. 2013
%%% For a given gene, generate an excel sheet with the average expression
%%% across all structures

function geneAvgExpPerStr(geneOfInterest, strOfInterest)

% define the files and results directory
filesDirectory = 'C:/Users/amahfouz/Documents/MATLAB/Data/NuclearReceptors/files/';
resultsDir = 'C:/Users/amahfouz/Documents/MATLAB/Results/NuclearReceptors/';
% load the average expression matrix
load([resultsDir 'strAvgExp/allGenes/strAvgExp.mat']);
% load the structure size matrix
load([resultsDir 'strAvgExp/allGenes/strSize.mat']);
% retrieve structures' names and ids
[num txt] = xlsread([filesDirectory 'ARAontology.xls']);
strNames = txt(2:end,5);
strIds = num(:,3);
clear num; clear txt;
% % load voxels annotations
load([filesDirectory 'voxelAnnotation.mat']);
% load the list of all genes and find the gene indicies
load([filesDirectory 'allGenes.mat']);
load([filesDirectory 'allExpNumbers.mat']);
load([filesDirectory 'allExpPlanes.mat']);
gene_index = find(strcmpi(allGenes, geneOfInterest) == 1);
gene_experimentNos = allExpNumbers(gene_index);
gene_experimentPlanes = allExpPlanes(gene_index);
% get the ids of the substructures of the structure of interest
subStr = indicateSubStr(strOfInterest, filesDirectory);
strIndecies = find(ismember(strIds, subStr) == 1);
% extract average expression data for the gene experiments
geneAvgExpMat = strAvgExp(strIndecies,gene_index);
% prepare col headers
for i = 1 : length(gene_index)
    rowNames{i} = [allGenes{gene_index(i)} ' ' gene_experimentNos{i} ' ' gene_experimentPlanes{i}];
end
% define the output file
outFolder = 'C:/Users/amahfouz/Documents/MATLAB/Results/NuclearReceptors/strAvgExp/';
outFile = [outFolder 'subAvgExp_' geneOfInterest '_' strOfInterest '.xlsx'];
% write data to the output file
xlswrite(outFile, {'name', 'ID', 'size(voxels)'}, 1, 'A1');
xlswrite(outFile, rowNames, 1, 'D1');
xlswrite(outFile, strNames(strIndecies), 1, 'A2');
xlswrite(outFile, strIds(strIndecies), 1, 'B2');
xlswrite(outFile, strSize(strIndecies)', 1, 'C2');
xlswrite(outFile, geneAvgExpMat, 1, 'D2');



