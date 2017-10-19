%%% 19 Feb. 2013
%%% save topCorrgenes results to an excel sheet
%%% also save a volume showing the structure of interest
%%% 23 April 2013: 'TYPE' is an option to run the correlation based on the
%%% coronal experiments only (C), sagittal only (S), or on both (MIXED)


function txtCorrRes(geneOfInterest, strOfInterest, TYPE, extension, ...
    filesFolder, resultsDir, geneCorr, xlsFiles, K)

% point to the files directory
% filesFolder = 'E:/Ahmed/HP/work/Data/NuclearReceptors/files/';

% load experiments' plane
load([filesFolder 'allExpPlanes.mat']);
% load gene names
load([filesFolder 'allGenes.mat']);
% load experiments' numbers
load([filesFolder 'allExpNumbers.mat']);
% load the average expression, localization and fitting score and keep only
% the structure of interest info
% resultsDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/NuclearReceptors/';
load([resultsDir 'fitScore.mat']);
load([resultsDir 'locScore.mat']);
load([resultsDir 'strAvgExp/allGenes/strAvgExp.mat']);
[num txt] = xlsread([filesFolder 'ARAontology.xls']);
strNames = txt(2:end,5); clear txt; clear num;
strIndex = find(strcmp(strNames, strOfInterest) == 1);

% define TYPE-specific variables
if strcmp(TYPE, 'C')
    resFolder = [geneCorr 'CoronalOnly/' ...
        geneOfInterest '/' geneOfInterest '_' strOfInterest '_CoronalOnly/'];
    if ~exist(resFolder, 'dir')
        display([resFolder '...does not exist']);
        return;
    else
        expDir = [xlsFiles '/CoronalOnly/'];
        if ~exist(expDir, 'dir')
            mkdir(expDir);
        end
        outDir = [expDir geneOfInterest '/'];
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end
        oF = [outDir geneOfInterest '_' strOfInterest '_CoronalOnly.' extension];
        allGenes(4346:end) = [];
        allExpNumbers(4346:end) = [];
        
        currStrAvgExp = strAvgExp(strIndex,1:4345);
        currStrAvgExp_BRAIN = strAvgExp(2,1:4345);
        currFitScore = fitScore(1:4345,strIndex);
        currLocScore = locScore(1:4345,strIndex);
    end
    
elseif strcmp(TYPE, 'All')
    resFolder = [geneCorr 'All/' ...
        geneOfInterest '/' geneOfInterest '_' strOfInterest '/'];
    if ~exist(resFolder, 'dir')
        display([resFolder '...does not exist']);
        return;
    else
%         expDir = ['E:\Ahmed\HP\work\Results\NuclearReceptors\xlsFiles/All/'];
        expDir = [xlsFiles 'All/'];
        if ~exist(expDir, 'dir')
            mkdir(expDir);
        end
        outDir = [expDir geneOfInterest '/'];
        if ~exist(outDir, 'dir')
            mkdir(outDir);
        end
        oF = [outDir geneOfInterest '_' strOfInterest '.' extension];
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

% for K = 1 : length(gene_experimentNos)
    
    % remove all useless values
    vToKeep = find(nmVoxels(K,:) >= 5);
    if ~isempty(vToKeep)
        cNEW = c(K,vToKeep);
        pvalNEW = pval(K,vToKeep);
        nmVoxelsNEW = nmVoxels(K,vToKeep);
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
        
        % make table
        T = table([1:length(allGenesNEW)]', allGenesNEW(sortingInd)',...
            allExpNumbersNEW(sortingInd)', allExpPlanesNEW(sortingInd)',...
            sortedCorrVals', pvalNEW(sortingInd)',...
            currStrAvgExpNEW(sortingInd)',...
            ((currStrAvgExpNEW(sortingInd)') ./ currcStrAvgExp_BRAIN_NEW(sortingInd)'),...
            currLocScoreNEW(sortingInd),...
            nmVoxelsNEW(sortingInd)', expRanks, relExpRank, relcorrRank, rankProd,...
            'VariableNames',{'rank', 'gene_symbol' 'experiment_no' 'plane'...
            'corr' 'p_value' 'avg_exp' 'norm_avg_exp' 'loc_score'...
            'number_of_voxels_used_in_the_computation' 'exp_rank'...
            'rel_exp_rank' 'rel_corr_rank' 'rank_prod'});
        writetable(T, oF);

        clear sortedCorrVals; clear sortingInd; 
    end
    % create new sheet names (experiment number_experiment plane)
%     newSheetNames{i} = [gene_experimentNos{i} '_' gene_experimentPlanes{i}];
% end
% rename the file sheets
% xlsheets(newSheetNames, oF);


