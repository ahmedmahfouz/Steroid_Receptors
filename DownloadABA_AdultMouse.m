%%% 12 Feb. 2013
%%% Download the ALLEN Brain Atlas :: Adult Mouse 

%% Retrieve experiments' numbers
% expMetaDataFile = 'C:\Users\amahfouz\Documents\MATLAB\Data\NuclearReceptors\files\ABA experiments IDs.csv';
expMetaDataFile = 'E:\Ahmed\HP\work\Data\NuclearReceptors\files\ABA experiments IDs.csv';
[num txt] = xlsread(expMetaDataFile);
sectionOrientation = txt(2:end,1);
genesList = txt(2:end,2);
sectionDataSet = num(:,2);
clear txt; clear num;

%% These are hard-coded paths to URLs for downloading expression volumes
API_SERVER = 'http://api.brain-map.org/';
GRID_FMT = [API_SERVER 'grid_data/download/'];

%% Download ALL coronal grid files (expression energy) - there are 4347
%%% coronal experiments
outDir = ['expressionData/coronal/'];
for expLoop = 1 : 4347
    experimentNumber = sectionDataSet(expLoop);
    fullURL = [GRID_FMT num2str(experimentNumber)];
    geneName = genesList{expLoop};
    %%% some genes have an '*' at the end -> remove
    if ~isempty(strfind(geneName, '*'));
        geneName = geneName(1:end-1);
    end
    expNumber = sectionDataSet(expLoop);
    geneFolderName = [outDir geneName '_' num2str(expNumber)];
    unzip(fullURL, geneFolderName);
    clear experimentNumber; clear geneName; clear expNumber; 
    clear geneFolderName;
end

%% Download ALL sagittal grid files (expression energy) - starting @ 4348
outDir = ['expressionData/sagittal/'];
for expLoop = 4348 : length(sectionDataSet)
    experimentNumber = sectionDataSet(expLoop);
    fullURL = [GRID_FMT num2str(experimentNumber)];
    geneName = genesList{expLoop};
    if ~isempty(strfind(geneName, '*'));
        geneName = geneName(1:end-1);
    end
    expNumber = sectionDataSet(expLoop);
    geneFolderName = [outDir geneName '_' num2str(expNumber)];
    try 
        unzip(fullURL, geneFolderName);
    catch err
    end
    clear experimentNumber; clear geneName; clear expNumber; 
    clear geneFolderName;
end

