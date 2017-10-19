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
% extract indices of labeled (i.e. brain) voxels from in the left
% hemisphere from an ARA mask (sent by Lydia) - that mask also filters
% voxels that have more than 20% missing data
maskVol = metaImageRead([filesDirectory 'gridAnnotation/image_mask.mhd']);
maskVol = permute(maskVol, [2 1 3]);
maskVol(maskVol > 0) = 1;
referenceVol = referenceVol .* double(maskVol);
brainInd = find(referenceVol ~= 0);
voxelAnnotation = referenceVol(brainInd);
save([filesDirectory 'voxelAnnotation.mat'], 'voxelAnnotation');
save([filesDirectory 'brainInd.mat'], 'brainInd');
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
    brainOnlyData = geneExperimentData(brainInd);
    clear geneExperimentData;
    expressionMatrix(:,i) = brainOnlyData';
    clear brainOnlyData;
    listOfGenes{i} = geneExperimentName(1 : strfind(geneExperimentName, '_')-1);
    listOfExperimentNumbers{i} = geneExperimentName(strfind(geneExperimentName, '_')+1 : end);
    clear geneExperimentName; clear fileName; clear geneExperimentName;
end
clear genesFolders;
save([filesDirectory 'expressionMatrix.mat'], 'expressionMatrix', '-v7.3');
save([filesDirectory 'listOfGenes.mat'], 'listOfGenes');
save([filesDirectory 'listOfExperimentNumbers.mat'], 'listOfExperimentNumbers');

% second load the sagittal genes
coronalExperimentsFolder = [dataFolder 'sagittal/'];
genesFolders = dir(coronalExperimentsFolder);
for i = 1 : length(genesFolders)-2
    geneExperimentName = genesFolders(i+2).name;
    fileName = [coronalExperimentsFolder geneExperimentName '/energy.raw'];
    fid = fopen(fileName, 'r', 'l');
    geneExperimentData = fread(fid, prod(SIZE), 'float');
    fclose(fid);
    geneExperimentData = reshape(geneExperimentData, SIZE);
    brainOnlyData = geneExperimentData(brainInd);
    clear geneExperimentData;
    expressionMatrix(:,i) = brainOnlyData';
    clear brainOnlyData;
    listOfGenes{i} = geneExperimentName(1 : strfind(geneExperimentName, '_')-1);
    listOfExperimentNumbers{i} = geneExperimentName(strfind(geneExperimentName, '_')+1 : end);
    clear geneExperimentName; clear fileName; clear geneExperimentName;
end
clear genesFolders;
save([filesDirectory 'expressionMatrix_sagittal.mat'], 'expressionMatrix', '-v7.3');
save([filesDirectory 'listOfGenes_sagittal.mat'], 'listOfGenes');
save([filesDirectory 'listOfExperimentNumbers_sagittal.mat'], 'listOfExperimentNumbers');


