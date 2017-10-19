%%% 13 Feb. 2013
%%% Arrange all gene expressions in one big matrix

clear all;

filesDirectory = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MatlabFiles/NuclearReceptors/files/';
dataFolder = [filesDirectory 'expressionData/'];

% get annotation information of voxels using the reference atlas at
% resolution = 200µm 
annoationGridFile = [filesDirectory 'gridAnnotation/gridAnnotation.raw'];
SIZE = [67 41 58]; % 200 micron volume size
fid = fopen(annoationGridFile, 'r', 'l');
referenceVol = fread(fid, prod(SIZE), 'uint16');
fclose(fid);
referenceVol = reshape(referenceVol, SIZE);
% extract indices of labeled (i.e. brain) voxels
brainInd_CoronalOnly = find(referenceVol ~= 0);
voxelAnnotation_CoronalOnly = referenceVol(brainInd_CoronalOnly);
save([filesDirectory 'voxelAnnotation_CoronalOnly.mat'], 'voxelAnnotation_CoronalOnly');
save([filesDirectory 'brainInd_CoronalOnly.mat'], 'brainInd_CoronalOnly');
clear referenceVol;

% first load the coronal genes
coronalExperimentsFolder = [dataFolder 'coronal/'];
genesFolders = dir(coronalExperimentsFolder);
for i = 1 : length(genesFolders)-2
    geneExperimentName = genesFolders(i+2).name;
    fileName = [coronalExperimentsFolder geneExperimentName '/energy.raw'];
    fid = fopen(fileName, 'r', 'l');
    geneExperimentData = fread(fid, prod(SIZE), 'float');
    fclose(fid);
    geneExperimentData = reshape(geneExperimentData, SIZE);
    brainOnlyData = geneExperimentData(brainInd_CoronalOnly);
    clear geneExperimentData;
    expressionMatrix_CoronalOnly(:,i) = brainOnlyData';
    clear brainOnlyData;
    listOfGenes{i} = geneExperimentName(1 : strfind(geneExperimentName, '_')-1);
    listOfExperimentNumbers{i} = geneExperimentName(strfind(geneExperimentName, '_')+1 : end);
    clear geneExperimentName; clear fileName; clear geneExperimentName;
end
clear genesFolders;
save([filesDirectory 'expressionMatrix_CoronalOnly.mat'], 'expressionMatrix_CoronalOnly', '-v7.3');

