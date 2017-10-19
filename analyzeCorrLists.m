%%% 24 April 2013
%%% Analyze top correlated gene lists and ranking of specific genes
%%% (e.g. coregulators)


%% define the files and results folders
if ispc
    filesDirectory = 'E:\Ahmed\HP\work\Data\NuclearReceptors\files\';
    resultsDirectory = 'E:\Ahmed\HP\work\Results\NuclearReceptors\';
else
    filesDirectory = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Data/NuclearReceptors/files/';
    resultsDirectory = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/NuclearReceptors/';
end
% define the NRs to analyze
% geneOfInterest = {'Ar', 'Nr3c1', 'Esr1', 'Esr2', 'Nr3c2', 'Pgr'};
geneOfInterest = {'Esr1'};
% define the structures to analyze
%%% Main brain structures 
mainStr = {'grey', 'Isocortex', 'OLF', 'HPF' 'CTXsp', 'STR', 'PAL', 'CB', ...
    'TH', 'HY', 'MB', 'P', 'MY'};
% %%% interesting structures
%%% HY substructures
% [num txt] = xlsread([filesDirectory 'HY_substructures.xls']);
% hySubStr = txt(2:end,3);
% hySubStr = [{'HY'} hySubStr'];
% clear num; clear txt;
% structures = hySubStr(2:5);
%%% final set from Onno (March 2014)
% structures = {'BST', 'CEA', 'sAMY', 'NTS', 'VTA', 'SNr', 'SNc', 'DR', 'ACB'};
structures_HIP = {'HIP', 'CA', 'DG', 'CA1', 'CA2', 'CA3'};
structures_PAL = {'PALd', 'PALv', 'PALm', 'PALc' 'BST'};
% structures = [mainStr, structures_HIP, structures_PAL];
structures = {'grey', 'HY'};
% structures = {'VTA', 'SNr', 'SNc'};
% structures = {'grey', 'Isocortex', 'OLF', 'HIP' 'CTXsp', 'STR', 'PAL', 'CB', ...
%     'TH', 'HY', 'MB', 'P', 'MY'};
% structures = mainStr;
% select a gene set ('coregulators' or 'grCoReg', xuCell2012, endo_CORT, CORT);
geneSetName = 'xuCell2012';
% select an experiment type ('C', 'S' or 'All')
expType = 'All';
% fileExtensions = {'controlTotal', 'controlOnly', 'controlTotalUp', 'controlTotalDown'...
%     'stressTotal', 'stressOnly', 'stressTotalUp', 'stressTotalDown', 'shared'};
% fileExtensions = {'controlOnlyUp', 'controlOnlyDown', 'stressOnlyUp', 'stressOnlyDown'...
%     'sharedUp', 'sharedDown', 'all'};
fileExtensions = {'all', 'control_Only', 'CRS_Only', 'shared'};

%% check the ranks of input genes
datasetFile = 'C:\Users\amahfouz\SURFdrive\Projects\Nuclear_Receptors\PNAS\Supplementary_Materials\Dataset_S1.xlsx';
strS = {'CB' , 'CTXsp', 'grey', 'HPF', 'HY', 'Isocortex', 'MB', 'MY', 'OLF', ...
    'P', 'PAL', 'STR', 'TH'};
% Esr1
gName = 'Esr1';
sheetNo = 2;
data = readtable(datasetFile, 'Sheet', sheetNo);
varNames = data.Properties.VariableNames;
gList = data.Var16(4:13);
% Gr
gName = 'Gr';
sheetNo = 4;
data = readtable(datasetFile, 'Sheet', sheetNo);
varNames = data.Properties.VariableNames;
gList = data.Var16(3:12);
% Mr
gName = 'Mr';
sheetNo = 5;
data = readtable(datasetFile, 'Sheet', sheetNo);
varNames = data.Properties.VariableNames;
gList = data.Var16(3:12);
% find the ranks of the genes in each list
for L = 1 : length(strS)
    currGeneList = data.(varNames{2+(7*(L-1))});
    for i = 1 : length(gList)
        tempVar = find(ismember(currGeneList,gList{i})==1,1);
        geneRank(i,L) = tempVar-1;
    end
end
% save
xlswrite([resultsDirectory gName '_top10.xlsx'], geneRank, 1, 'B2');
xlswrite([resultsDirectory gName '_top10.xlsx'], strS, 1, 'B1');
xlswrite([resultsDirectory gName '_top10.xlsx'], gList, 1, 'A2');

%% check overlap with known TFBS
promdist = 10000;
noOfGenes = 2000;
selection = 'top';
encodeFLAG =  0; 
wgrvistaFLAG = 0; 
premodFLAG = 0;
swissregulonFLAG = 0; 
pwmscanFLAG = 0;
ecrbaseFLAG = 1;
for S = 1 : length(structures)
    TFBSanalyze(filesDirectory, resultsDirectory, ...
        geneOfInterest, expType, structures{S}, ...
        promdist, noOfGenes, selection, ...
        encodeFLAG, wgrvistaFLAG, premodFLAG, swissregulonFLAG,pwmscanFLAG, ecrbaseFLAG);
end

%% group results in one excel sheet
for i = 1 : length(geneOfInterest)
    numExp = groupResults(filesDirectory, resultsDirectory, geneOfInterest(i), structures, ...
        geneSetName, expType);
end

%% group results in one excel sheet for GR coreg
for i = 1 : length(fileExtensions)
    numExp = groupResults(filesDirectory, resultsDirectory, geneOfInterest(1), structures, ...
        [geneSetName '_' fileExtensions{i}], expType);
end

%% Analyze CORT results 10 Dec 2014
% for i = 1 : length(fileExtensions)
%     [num txt] = xlsread([resultsDirectory 'CORT_' fileExtensions{i} '_Nr3c1_All.xls'],2); % read sheet 2, GR_sagittal_728
%     corrRank(:,:,i) = num(:,1:3:end);
%     corr(:,:,i) = num(:,2:3:end);
% end
% % boxplot per brain region
% clear Data; Data(:,:) = corr(:,:,1);
% Data(sum(Data') == 0,:) = [];
% figure, boxplot(Data, 'labels', structures)

load([filesDirectory 'allGenes.mat']);
load([filesDirectory 'allExpNumbers.mat']);
load([filesDirectory 'allExpPlanes.mat']);
offset = 0;
for i = 1 : length(fileExtensions)
    geneSetName = ['CORT_' fileExtensions{i}];
    P = geneSetRankEval(filesDirectory, resultsDirectory, geneOfInterest, ...
        structures', geneSetName, expType, allGenes, allExpNumbers, allExpPlanes);
%     R = xlcalcrange('B2', offset, 0, 1, length(P));
%     offset = offset+1;
%     xlswrite([resultsDirectory 'CORT_results_July2015_allSTR.xlsx'], P(2,:), 1, R);
%     xlswrite([resultsDirectory 'CORT_rankSumpValCorrected_extra'], {geneSetName}, 1, 'A2');
    pFinal(i,:) = P(2,:);
end
% xlswrite([resultsDirectory 'CORT_rankSumpValCorrected'], structures, 1, 'B1');
xlswrite([resultsDirectory 'CORT_results_July2015_allSTR.xlsx'], pFinal, 1, 'B2');
save([resultsDirectory 'CORT_results_July2015_allSTR.mat'], 'pFinal')

%% Assess sinificance of co-expression
for g = 1 : length(geneOfInterest)
    if strcmp(expType, 'C')
        resultsFilename = [resultsDirectory geneSetName '_' geneOfInterest{g} '_CoronalOnly.xls'];
        n = 4345;
    elseif strcmp(expType, 'All')
        resultsFilename = [resultsDirectory geneSetName '_' geneOfInterest{g} '_All.xls'];
        n = 26020;
    end
    for E = 1 : numExp
        clear num; clear txt; clear rankSum; clear pVal; 
        [num txt] = xlsread(resultsFilename, E);
        ranks = num(:, 1:3:(3*length(structures)-2));
        rankSum = sum(ranks,2);
        % randomly generate rank sums
        nPerm = 10000;
        for i = 1 : nPerm
            randRankSum(i) = randi(n,1) + randi(n,1) + randi(n,1);
        end
        for i = 1 : length(rankSum)
            pVal(i) = (length(find(randRankSum <= rankSum(i))) + 1) / nPerm;
        end
        Z = find(rankSum == 0); pVal(Z) = 0;
        
%         % plot the null distribution and the siginificant gene for one
%         % experiment
%         figure, hist(randRankSum,100);
%         h = findobj(gca,'Type','patch');
%         set(h,'FaceColor',[80/255 80/255 80/255],'EdgeColor','w')
%         grid on;
%         X = find(pVal < 0.01);
%         X(ismember(X,find(rankSum == 0))) = [];
%         G = txt(3:end,1);
%         for xx = 1 : length(X)
%             line([rankSum(X(xx)) rankSum(X(xx))], [0 50*xx], 'linewidth', 2, 'color', 'black', 'linestyle', '--')
%             text(rankSum(X(xx))+500, 50*xx, G{X(xx)}, 'Color', 'black', 'FontAngle', 'italic', 'FontSize', 15)
%         end
%         xlabel('RankSum', 'FontWeight', 'bold', 'FontSize', 15)
%         ylabel('Count', 'FontWeight', 'bold', 'FontSize', 15)
%         title('RankSum PDF for 3 structures', 'FontWeight', 'bold', 'FontSize', 15);

    %     correctedP = pVal * (length(rankSum)-length(Z));
        R1 = xlcalcrange('A3', 0, 3*length(structures)+1, length(rankSum), 1);
        xlswrite(resultsFilename, rankSum, E, R1);
        R2 = xlcalcrange('A3', 0, 3*length(structures)+2, length(rankSum), 1);
        xlswrite(resultsFilename, pVal', E, R2);
    %     R3 = xlcalcrange('A3', 0, 3*length(structures)+3, length(rankSum), 1);
    %     xlswrite(resultsFilename, correctedP', E, R3);
        R4 = xlcalcrange('A1', 0, 3*length(structures)+1, 1, 2);
        xlswrite(resultsFilename, [{'Rank Sum'}, {'p-value'}], E, R4);
    end
end

%% group results in one excel sheet
load([filesDirectory 'allGenes.mat']);
load([filesDirectory 'allExpNumbers.mat']);
load([filesDirectory 'allExpPlanes.mat']);
singleNrResults(filesDirectory, resultsDirectory, geneOfInterest, structures, ...
    geneSetName, expType, allGenes, allExpNumbers, allExpPlanes);

%% plot ranks and correlations
load([filesDirectory 'allGenes.mat']);
load([filesDirectory 'allExpNumbers.mat']);
load([filesDirectory 'allExpPlanes.mat']);
plotcorrs_ranks(filesDirectory, resultsDirectory, geneOfInterest, structures, ...
    geneSetName, expType, allGenes, allExpNumbers, allExpPlanes);

%% calculate the correlation between a gene of interest and a set of genes 
% based on the expression of the main structures. Plot the average
% expression profiles
load([filesDirectory 'allGenes.mat']);
load([filesDirectory 'allExpNumbers.mat']);
load([filesDirectory 'allExpPlanes.mat']);
load([resultsDirectory 'strAvgExp/allGenes/strAvgExp.mat']);
load([resultsDirectory 'strAvgExp/allGenes/strSize.mat']);
corrWithAvgExp(filesDirectory, resultsDirectory, geneOfInterest, structures, geneSetName,...
    allGenes, allExpNumbers, allExpPlanes, strAvgExp, strSize);

%% calculate a p-value for the rank sum of a gene set
load([filesDirectory 'allGenes.mat']);
load([filesDirectory 'allExpNumbers.mat']);
load([filesDirectory 'allExpPlanes.mat']);
geneSetRankEval(filesDirectory, resultsDirectory, geneOfInterest, ...
    structures', geneSetName, expType, allGenes, allExpNumbers, allExpPlanes);

%% calculate a p-value for the rank sum of xu et al with all NRs
load([filesDirectory 'allGenes.mat']);
load([filesDirectory 'allExpNumbers.mat']);
load([filesDirectory 'allExpPlanes.mat']);
for i = 1 : length(geneOfInterest)
    gene{1} = geneOfInterest{i};
    [P{i}, N{i}] = geneSetRankEval(filesDirectory, resultsDirectory, gene, ...
        structures', geneSetName, expType, allGenes, allExpNumbers, allExpPlanes);
end
save([resultsDirectory 'xuCell2012_allNRs.mat'],'P','N')

%% differential coexpression analysis
load([filesDirectory 'allGenes.mat']);
load([filesDirectory 'allExpNumbers.mat']);
load([filesDirectory 'allExpPlanes.mat']);
%%% structure1 is the structure of interest (HY in case of Esr1) an
%%% structure2 is the "control" structure (Isocortex in case of Esr1)
structure1 = 'CEA';
structure2 = 'BST';
diffCoexp(filesDirectory, resultsDirectory, geneOfInterest, ...
    structure1, structure2, geneSetName, expType, allGenes, allExpNumbers, allExpPlanes);

%% plot the results of the ranksum analysis of the DCW
strs = {'HY', 'CB', 'Isocortex'};
load([filesDirectory 'allGenes.mat']);
load([filesDirectory 'allExpNumbers.mat']);
load([filesDirectory 'allExpPlanes.mat']);
plotDCWresults(filesDirectory, resultsDirectory, geneOfInterest, strs, ...
    geneSetName, expType, allGenes, allExpNumbers, allExpPlanes);

%% ESR1 target genes test set
[num txt] = xlsread('C:\Users\amahfouz\Documents\MATLAB\Results\NuclearReceptors\motifEnrichment\ESR1_testTargets.xlsx');
esr1TargetGenes = txt(2:end);
clear num; clear txt;
upStream = 2000;
motifLength = 7;
alphaSign = 0.05;
outDir = [resultsDirectory 'motifEnrichment/fastaFiles/'];
fileName = ['ESR1_TargetGenes' ...
        '_' num2str(upStream/1000) 'kb' '_' 'query' '.fa'];
species = 'mm10';
motif = 'AGGTCACCRTGRCCY';
UCSCdir = ['C:\Users\amahfouz\Documents\MATLAB\Results\MotifEnrichment\New folder\UCSC_DB\' species '\'];
geneListFasta(esr1TargetGenes, species, upStream, motif, ...
        UCSCdir, outDir, fileName);
[queryHeader, querySeq] = fastaread([outDir fileName]);
outDir = [resultsDirectory 'motifEnrichment/' 'ESR1_TargetGenes/'];
if ~exist(outDir)
    mkdir(outDir);
end
%database information
if ispc
    jaspartransfac.Dir = 'C:\Users\amahfouz\Documents\MATLAB\Scripts\motifenrichment\data\';
else
    jaspartransfac.Dir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Scripts/motifenrichment/data/';
end
jaspartransfac.Ext = '.pwm';
jaspartransfac.Name = {'PWMs_union'; 'JASPAR_vertebrates'; 'TRANSFAC'};
jaspartransfac.Corr = 0.6;
%run the motif enrichment analysis
outMotifsRegion = motifsearch(querySeq, ...
        'motiflength', motifLength, ...
        'alpha', alphaSign, ...
        'multipletest', 2, ...
        'strandboth', 1, ...
        'writedir', outDir, ...
        'filterrepeat', 1, ...
        'fontsize', 8, 'fontname', 'Helvetica', 'figname', 'CIS', ...
        'filename', 'CIS',...
        'jaspartransfacDir', jaspartransfac.Dir, ...
        'jaspartransfacExt', jaspartransfac.Ext, ...
        'jaspartransfacName', jaspartransfac.Name, ...
        'jaspartransfacCorr', jaspartransfac.Corr, ...
        'printloaddata', 0);

%% generate fasta files for motif enrichment analysis
upStream = 2000;
noOfGenes = 200;
load([filesDirectory 'allGenes.mat']);
load([filesDirectory 'allExpNumbers.mat']);
load([filesDirectory 'allExpPlanes.mat']);
addpath(genpath('C:\Users\amahfouz\Documents\MATLAB\Scripts\motifenrichment'));
for S = 1 : length(structures)
    geneMotifEnrichFiles(filesDirectory, resultsDirectory, geneOfInterest, ...
        structures{S}, expType, upStream, noOfGenes, allGenes, allExpNumbers, allExpPlanes);
end

%% do the motif enrichment analysis
upStream = 2000;
noOfGenes = 200;
motifLength = 7;
alphaSign = 0.05;
for S = 1 : length(structures)
    motifEnrichAnalysis(filesDirectory, resultsDirectory, geneOfInterest, ...
        structures{S}, expType, upStream, noOfGenes, motifLength, alphaSign);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% analyze a geneSet
% [~, txt] = xlsread([filesDirectory 'gene sets/' geneSetName '.xlsx']);
% geneSet = txt(2:end,1);
% clear txt;
% 
% [num txt] = xlsread([filesDirectory 'GR-coregulatorsIP.xlsx']);
% pList = txt(2:end,1);
% pScores = num(:,3);
% clear txt; clear num;
% for i = 1 : length(pList)
%     coRegTemp = pList{i};
%     coReg{i} = coRegTemp(1:strfind(coRegTemp, '_')-1);
%     pos = find(strcmpi(coReg{i}, geneSet) == 1);
%     if ~isempty(pos)
%         coRegScore(pos) = pScores(i);
%     else
%         coRegScore(pos) = 0;
%     end
% end
% % coRegList = unique(coReg);
% 

%% [] prepare files for the GSEA in R
for i = 1 : 2%length(genesOfInterest)
    geneSetRank(geneOfInterest{i}, geneSetName, structures, expType);
end

%% combine corr/rank matrices of each geneSet in one matlab matrix
for i = 1 : length(geneOfInterest)
    [gS gC(:,:,i) gR(:,:,i)] = geneSetAnalysis(filesDirectory, resultsDirectory, ...
        geneOfInterest{i}, structures, geneSetName, expType);
end
save(['results/' geneSetName '/' geneSetName '_' expType '_corr_allSTR.mat'], 'gC');
save(['results/' geneSetName '/' geneSetName '_' expType '_rank_allSTR.mat'], 'gR');
save(['results/' geneSetName '/' geneSetName '_' expType '_genes_allSTR.mat'], 'gS');

%% cluster correlation/rank patterns of geneSet with NRs
load(['results/' geneSetName '/' geneSetName '_' expType '_corr_smallSTR.mat']);
load(['results/' geneSetName '/' geneSetName '_' expType '_rank_smallSTR.mat']);
load(['results/' geneSetName '/' geneSetName '_' expType '_genes_smallSTR.mat']);
for i = 1 : size(gS,1)
    clear dataToPlot;
    dataToPlot(:,:) = -log10(gR(i,:,:)./4345);
%     dataToPlot = repmat(4345, size(gR,2), size(gR,3)) - dataToPlot;
    f = figure;
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    bar(dataToPlot, 'histc'); grid on; colormap jet
    set(gca, 'XLim', [0 size(dataToPlot,1)+1]);
    set(gca, 'XTickLabel', structures, 'FontWeight', 'bold', 'XTick', 1:numel(structures));
    rotateXLabels(gca, 45)
    if strcmpi(expType, 'C')
        N = 4345;   
    elseif strcmpi(expType, 'S')
        N = 21677;
    elseif strcmpi(expType, 'All')
        N = 26022;
    end
    set(gca, 'YLim', [0 -log10(1/N)]);
    line([0 size(dataToPlot,1)+1], [-log10(10/N) -log10(10/N)], ...
        'linewidth', 2, 'color', 'r');
    text(0.1, -log10(10/N)+0.05, 'Top 10', 'FontWeight', 'bold', ...
        'color', 'r');
    line([0 size(dataToPlot,1)+1], [-log10(100/N) -log10(100/N)], ...
        'linewidth', 2, 'color', 'r');
    text(0.1, -log10(100/N)+0.05, 'Top 100', 'FontWeight', 'bold', ...
        'color', 'r');
    title(gS(i,1), 'FontSize', 15, 'FontWeight', 'bold');
%     xlabel('Brain Structures', 'FontSize', 15, 'FontWeight', 'bold');
    ylabel('-log_1_0(Rank / Total Number of Genes)', 'FontSize', 15, 'FontWeight', 'bold');
    legend(geneOfInterest, 'FontWeight', 'bold');
    saveas(f, ['results/' geneSetName '/' geneSetName '_' expType '_' gS{i} '_relR.fig']);
    saveas(f, ['results/' geneSetName '/' geneSetName '_' expType '_' gS{i} '_relR.jpg']);
    close(f);
end

%% Analyze Nr3c1/Gr interaction prediction data
geneSetAnalysis(filesDirectory, resultsDirectory, geneOfInterest{1}, structures, ...
        geneSetName, expType);
    
[num txt] = xlsread([filesDirectory 'GR-coregulatorsIP.xlsx']);
origAvgDMSO = num(:,1);
origAvgDEX = num(:,2);
origIPScore = num(:,3); %abs(DEX-DMSO)
originalCoregList = txt(2:end, 2);
clear num; clear txt;
figure, bar([1:length(origIPScore)], origIPScore);
grid on
set(gca, 'XLim', [1 length(origIPScore)+1]);

%% Plot correlation results (3 July 2013)
load(['results/' geneSetName '/' geneSetName '_' expType '_corr_allSTR.mat']);
load(['results/' geneSetName '/' geneSetName '_' expType '_rank_allSTR.mat']);
load(['results/' geneSetName '/' geneSetName '_' expType '_genes_allSTR.mat']);
clear dataToPlot; clear tempData;
dataToPlot(:,:) = -log10(gR(:,:,1)./4345);
for i = 2 : size(gR,3)
    clear tempData;
    tempData(:,:) = -log10(gR(:,:,i)./4345);
    dataToPlot = [dataToPlot; tempData];
end
f = figure;
imagesc(dataToPlot); %colormap('redbluecmap')
set(gca, 'YTickLabel', repmat(gS, 5, 1), 'FontWeight', 'bold', 'FontSize', 7, ...
    'YTick', 1:numel(gS)*length(geneOfInterest));
set(gca, 'XTickLabel', structures, 'FontWeight', 'bold', 'XTick', 1:numel(structures));
rotateXLabels(gca, 45)
title({[geneSetName ' - ' expType], '-log_1_0(Rank / Total Number of Genes)'}, 'FontSize', 15, 'FontWeight', 'bold');
% lines to separate the NRs
line([0 size(dataToPlot,1)+1], [length(gS)+0.5 length(gS)+0.5], 'linewidth', 2, 'color', 'white');
line([0 size(dataToPlot,1)+1], [length(gS)*2+0.5 length(gS)*2+0.5], 'linewidth', 2, 'color', 'white');
line([0 size(dataToPlot,1)+1], [length(gS)*3+0.5 length(gS)*3+0.5], 'linewidth', 2, 'color', 'white');
line([0 size(dataToPlot,1)+1], [length(gS)*4+0.5 length(gS)*4+0.5], 'linewidth', 2, 'color', 'white');

% saveas(f, ['results/' geneSetName '/' geneSetName '_' expType '_' gS{i} '_relR.fig']);
% saveas(f, ['results/' geneSetName '/' geneSetName '_' expType '_' gS{i} '_relR.jpg']);







