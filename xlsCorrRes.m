%%% 19 Feb. 2013
%%% save topCorrgenes results to an excel sheet
%%% also save a volume showing the structure of interest
%%% 23 April 2013: 'TYPE' is an option to run the correlation based on the
%%% coronal experiments only (C), sagittal only (S), or on both (MIXED)


function xlsCorrRes(geneOfInterest, strOfInterest, TYPE)

% point to the files directory
filesFolder = 'E:/Ahmed/HP/work/Data/NuclearReceptors/files/';

% load experiments' plane
load([filesFolder 'allExpPlanes.mat']);
% load gene names
load([filesFolder 'allGenes.mat']);
% load experiments' numbers
load([filesFolder 'allExpNumbers.mat']);
% load the average expression, localization and fitting score and keep only
% the structure of interest info
resultsDir = 'E:\Ahmed\HP\work\Results\NuclearReceptors\';
load([resultsDir 'fitScore.mat']);
load([resultsDir 'locScore.mat']);
load([resultsDir 'strAvgExp/allGenes/strAvgExp.mat']);
[num txt] = xlsread([filesFolder 'ARAontology.xls']);
strNames = txt(2:end,5); clear txt; clear num;
strIndex = find(strcmpi(strNames, strOfInterest) == 1);

% define TYPE-specific variables
if strcmp(TYPE, 'C')
    resFolder = ['E:\Ahmed\HP\work\Results\NuclearReceptors\geneCorr/CoronalOnly/' ...
        geneOfInterest '/' geneOfInterest '_' strOfInterest '_CoronalOnly/'];
    if ~exist(resFolder, 'dir')
        display([resFolder '...does not exist']);
        return;
    else
        expDir = ['E:\Ahmed\HP\work\Results\NuclearReceptors\xlsFiles/CoronalOnly/'];
        if ~exist(expDir, 'dir')
            mkdir(expDir);
        end
        outDir = [expDir geneOfInterest '/'];
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end
        oF = [outDir geneOfInterest '_' strOfInterest '_CoronalOnly.xls'];
        allGenes(4346:end) = [];
        allExpNumbers(4346:end) = [];
        
        currStrAvgExp = strAvgExp(strIndex,1:4345);
        currStrAvgExp_BRAIN = strAvgExp(2,1:4345);
        currFitScore = fitScore(1:4345,strIndex);
        currLocScore = locScore(1:4345,strIndex);
    end
    
    
elseif strcmp(TYPE, 'S')
    resFolder = ['E:\Ahmed\HP\work\Results\NuclearReceptors/geneCorr/SagittalOnly/' ...
        geneOfInterest '/' geneOfInterest '_' strOfInterest '_SagittalOnly/'];
    if ~exist(resFolder, 'dir')
        display([resFolder '...does not exist']);
        return;
    else
        expDir = ['E:\Ahmed\HP\work\Results\NuclearReceptors\xlsFiles/SagittalOnly/'];
        if ~exist(expDir, 'dir')
            mkdir(expDir);
        end
        outDir = [expDir geneOfInterest '/'];
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end
        oF = [outDir geneOfInterest '_' strOfInterest '_SagittalOnly.xls'];
        allGenes(1:4345) = [];
        allExpNumbers(1:4345) = [];
        currStrAvgExp = strAvgExp(strIndex,4346:end);
        currStrAvgExp_BRAIN = strAvgExp(2,4346:end);
        currFitScore = fitScore(4346:end,strIndex);
        currLocScore = locScore(4346:end,strIndex);
    end
    
elseif strcmp(TYPE, 'All')
    resFolder = ['E:\Ahmed\HP\work\Results\NuclearReceptors\geneCorr/All/' ...
        geneOfInterest '/' geneOfInterest '_' strOfInterest '/'];
    if ~exist(resFolder, 'dir')
        display([resFolder '...does not exist']);
        return;
    else
        expDir = ['E:\Ahmed\HP\work\Results\NuclearReceptors\xlsFiles/All/'];
        if ~exist(expDir, 'dir')
            mkdir(expDir);
        end
        outDir = [expDir geneOfInterest '/'];
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end
        oF = [outDir geneOfInterest '_' strOfInterest '.xls'];
        currStrAvgExp = strAvgExp(strIndex,:);
        currStrAvgExp_BRAIN = strAvgExp(2,:);
        currFitScore = fitScore(:,strIndex);
        currLocScore = locScore(:,strIndex);
    end
    
else
    display('results folder not found');
    return;
end
    
% load a saved correlation vector
load([resFolder geneOfInterest '_' strOfInterest '.mat']);
% load a saved correlation p-value
load([resFolder geneOfInterest '_' strOfInterest '_PVAL.mat']);
% load the number of voxels used in each computation
load([resFolder strOfInterest '_nmVoxels.mat']);
% load the gene experiment numbers
load([resFolder geneOfInterest '_expNumbers.mat']);
load([resFolder geneOfInterest '_expPlane.mat']);

nanVals = find(isnan(c(1,:)));
if length(gene_experimentNos) > 1
    for k = 2 : length(gene_experimentNos)
        tempNanVals = find(isnan(c(k,:)));
        nanVals = [nanVals tempNanVals];
    end
end
expNAN = find(isnan(currStrAvgExp));
nanVals = [nanVals expNAN];
nanVals = unique(nanVals);
c(:,nanVals) = [];
pval(:,nanVals) = [];
nmVoxels(:,nanVals) = [];
allGenes(nanVals) = [];
allExpNumbers(nanVals) = [];
allExpPlanes(nanVals) = [];
currStrAvgExp(nanVals) = [];
currStrAvgExp_BRAIN(nanVals) = [];
currFitScore(nanVals) = [];
currLocScore(nanVals) = [];

for i = 1 : length(gene_experimentNos)
    
    % remove all useless values
    vToKeep = find(nmVoxels(i,:) >= 5);
    if ~isempty(vToKeep)
        cNEW = c(i,vToKeep);
        pvalNEW = pval(i,vToKeep);
        nmVoxelsNEW = nmVoxels(i,vToKeep);
        allGenesNEW = allGenes(vToKeep);
        allExpNumbersNEW = allExpNumbers(vToKeep);
        allExpPlanesNEW = allExpPlanes(vToKeep);
        currStrAvgExpNEW = currStrAvgExp(vToKeep);
        currcStrAvgExp_BRAIN_NEW = currStrAvgExp_BRAIN(vToKeep);
        currFitScoreNEW = currFitScore(vToKeep);
        currLocScoreNEW = currLocScore(vToKeep);

        % sort the correlation values from largest to smallest
        [sortedCorrVals sortingInd] = sort(cNEW, 'descend');
        % calculate expression rank and rank product
        [expRanks, ~] = tiedrank((currStrAvgExpNEW(sortingInd)') ./ currcStrAvgExp_BRAIN_NEW(sortingInd)');
        relExpRank = -log10((length(allGenesNEW)-expRanks+1)/length(expRanks));
        relcorrRank = -log10([1:length(allGenesNEW)]'/length(allGenesNEW));
        rankProd = relExpRank .* relcorrRank;

        % write all the headers in the excel file
        xlswrite(oF, {['seed gene: ' geneOfInterest '_' gene_experimentNos{i} '_' gene_experimentPlanes{i}], ['structure oif interest: ' strOfInterest]}, i, 'A1');   
        xlswrite(oF, {'rank', 'gene symbol', 'experiment no.', 'plane', ...
            'corr.', 'p-value', 'avg. exp.', 'norm. avg. exp.', 'loc. score' ...
            'number of voxels used in the computation', 'exp. rank', ...
            'rel. exp. rank', 'rel. corr. rank', 'rank prod.'}, i, 'A2');

        % write the numeric data
        xlswrite(oF, [[1:length(allGenesNEW)]' zeros(length(allGenesNEW),1) ...
            zeros(length(allGenesNEW),1) zeros(length(allGenesNEW),1)...
            sortedCorrVals' pvalNEW(sortingInd)' ...
            currStrAvgExpNEW(sortingInd)' ...
            ((currStrAvgExpNEW(sortingInd)') ./ currcStrAvgExp_BRAIN_NEW(sortingInd)') ...
            currLocScoreNEW(sortingInd) ...
            nmVoxelsNEW(sortingInd)' expRanks relExpRank relcorrRank rankProd], i, 'A3');

        % write the textual data
        xlswrite(oF, allGenesNEW(sortingInd)', i, 'B3');
        xlswrite(oF, allExpNumbersNEW(sortingInd)', i, 'C3');
        xlswrite(oF, allExpPlanesNEW(sortingInd)', i, 'D3');

        clear sortedCorrVals; clear sortingInd; 
    end
    % create new sheet names (experiment number_experiment plane)
    newSheetNames{i} = [gene_experimentNos{i} '_' gene_experimentPlanes{i}];
end
% rename the file sheets
xlsheets(newSheetNames, oF);


