%%% 6 Sep. 2013
%%% For a given gene and (large) structure, return the average expression
%%% level within all substructures.

function avgExpSubStr(geneOfInterest, strOfInterest)

% % filesDirectory = 'C:/Users/amahfouz/Documents/MATLAB/Data/NuclearReceptors/files/';
% filesDirectory = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MatlabFiles/NuclearReceptors/files/';
% 
% % load the expression matrix
% load([filesDirectory 'fullExpressionMatrix.mat']);
% % load gene names
% load([filesDirectory 'allGenes.mat']);
% % load experiments' numbers
% load([filesDirectory 'allExpNumbers.mat']);
% % load experiments' plane
% load([filesDirectory 'allExpPlanes.mat']);
% % load voxels annotations
% load([filesDirectory 'voxelAnnotation.mat']);
% 
% % get the list of voxels covering the structure of interest
% inclVoxAnnotations = indicateSubStr(strOfInterest, filesDirectory);
% % get structures' names
% [num txt] = xlsread([filesDirectory 'ARAontology.xls']);
% strNames = txt(2:end,5);
% strIDS = num(:,3);
% strIndecies = find(ismember(strIDS, inclVoxAnnotations) == 1);
% currStrIDS = strIDS(strIndecies);
% currStrName = strNames(strIndecies);
% 
% gene_index = find(strcmpi(allGenes, geneOfInterest) == 1);
% gene_experimentNos = allExpNumbers(gene_index);
% gene_experimentPlanes = allExpPlanes(gene_index);
% 
% if isempty(gene_index)
%     display('seed gene not found');
% else
%     % for each of the substructures, calculate the average expression and
%     % retain the number of voxels (size)
%     for subS = 1 : length(inclVoxAnnotations)
%         tempVoxAnnotations = indicateSubStr(currStrName{subS}, filesDirectory);
%         subStrVoxels = find(ismember(voxelAnnotation, tempVoxAnnotations) == 1);
%         subStrExpMat = fullExpressionMatrix(subStrVoxels, gene_index);
%         subStrAvgExp(subS,:) = mean(subStrExpMat);
%         subStrSize(subS) = length(subStrVoxels);
%     end
% end
% 
% resultsDir = ['/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MatlabFiles/NuclearReceptors/results/strAvgExp/'];
% save([resultsDir 'subStrSize_' geneOfInterest '_' strOfInterest '.mat'], 'subStrSize');
% save([resultsDir 'subStrAvgExp_' geneOfInterest '_' strOfInterest '.mat'], 'subStrAvgExp');
% save([resultsDir 'currStrName_' geneOfInterest '_' strOfInterest '.mat'], 'currStrName');
% save([resultsDir 'currStrIDS_' geneOfInterest '_' strOfInterest '.mat'], 'currStrIDS');
% save([resultsDir 'gene_experimentNos_' geneOfInterest '_' strOfInterest '.mat'], 'gene_experimentNos');


% write the results out in an excel sheet

resultsDir = ['C:/Users/amahfouz/Documents/MATLAB/Results/NuclearReceptors/strAvgExp/'];
load([resultsDir 'subStrSize_' geneOfInterest '_' strOfInterest '.mat']);
load([resultsDir 'subStrAvgExp_' geneOfInterest '_' strOfInterest '.mat']);
load([resultsDir 'currStrName_' geneOfInterest '_' strOfInterest '.mat']);
load([resultsDir 'currStrIDS_' geneOfInterest '_' strOfInterest '.mat']);
load([resultsDir 'gene_experimentNos_' geneOfInterest '_' strOfInterest '.mat']);

outFolder = 'C:/Users/amahfouz/Documents/MATLAB/Results/NuclearReceptors/strAvgExp/';
outFile = [outFolder 'subStrAvgExp_' geneOfInterest '_' strOfInterest];

xlswrite(outFile, {'name', 'ID', 'size(voxels)'}, 1, 'A1');
xlswrite(outFile, gene_experimentNos, 1, 'D1');
xlswrite(outFile, currStrName, 1, 'A2');
xlswrite(outFile, currStrIDS, 1, 'B2');
xlswrite(outFile, subStrSize', 1, 'C2');
xlswrite(outFile, subStrAvgExp, 1, 'D2');



