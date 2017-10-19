%%% 27 Sep 2013
%%% calculate the ranksum p-value for a gene set within all genes 

function [pValCorr_corrected, newSheetNames] = geneSetRankEval(filesDirectory, resultsDirectory, geneOfInterest, ...
    structures, geneSetName, expType, allGenes, allExpNumbers, allExpPlanes)

[~, txt] = xlsread([filesDirectory 'gene sets/' geneSetName '.xlsx']);
geneSet = txt(2:end,1);
clear txt;
geneSetInd = find(ismember(allGenes, geneSet) == 1);
geneSetNames = allGenes(geneSetInd);
geneSetExpNo = allExpNumbers(geneSetInd);
geneSetExpPlanes = allExpPlanes(geneSetInd);
%%% create the outfile for this analysis
outDir = [resultsDirectory geneSetName];
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
if strcmp(expType, 'C')
    outFile = [outDir  '/' geneSetName '_' geneOfInterest{1} 'rank-pVal_C.xls'];
    xlsDIR = [resultsDirectory 'xlsFiles/CoronalOnly/'];
    gene_index = find(strcmpi(allGenes(1:4345), geneOfInterest) == 1);
    extension = '_CoronalOnly';
    nGenes = 4345;
elseif strcmp(expType, 'All')
    outFile = [outDir  '/' geneSetName '_' geneOfInterest{1} 'rank-pVal_All.xls'];
    xlsDIR = [resultsDirectory 'xlsFiles/All/'];
    gene_index = find(strcmpi(allGenes, geneOfInterest) == 1);
    extension = '';
    nGenes = 26022;
else
    display('expType not defined correctly');
end
%%% find all experiments of the geneOfInterest
gene_experimentNos = allExpNumbers(gene_index);
gene_experimentPlanes = allExpPlanes(gene_index);
%%% loop on all the gene experiments
for experiment = 1 : length(gene_index)
    newSheetNames{experiment} = [allGenes{gene_index(experiment)} '_' gene_experimentNos{experiment} '_' gene_experimentPlanes{experiment}];
    %%% loop on all the structures in the structures list
    for s = 1 : length(structures)
        %%% load the file corresponding to the current structure and
        %%% gene experiment 
        [num txt] = xlsread([xlsDIR geneOfInterest{1} '/' geneOfInterest{1} '_' structures{s} extension '.xls'], experiment);
        rankedGenes = txt(3:end, 2);
        rankedExpNo = num(1:end,3);
        rankedExpPlane = txt(3:end, 4);
        rankedCorrVals = num(1:end,5);
        avgExpVals = num(1:end,7);
        locScores = num(1:end,8);
        % generate sorted gene lists based on the average expression
        % and localization scores for the current structure
        [sortedAvgExpVals sortingInd_avgExpVals] = sort(avgExpVals, 'descend');%%% watch out for tied ranks!!!
        [sortedCLocScores sortingInd_locScores] = sort(locScores, 'descend');
        clear txt; clear num;
        % return the corr, avg. epression and the localization score
        % for all the experiments of each gene
        indecies = 0;
        for g = 1 : length(geneSetNames)
            % find the index of the gene within the ranked list
            tempR = find(ismember(rankedGenes, geneSetNames{g}) == 1);
            if ~isempty(tempR)
                tempRE = find(rankedExpNo(tempR) == str2double(geneSetExpNo{g}));
                if ~isempty(tempRE)
                    currGeneSetInd(g) = tempR(tempRE);
                    geneCorr(g) = rankedCorrVals(tempR(tempRE));
                    geneCorrRank(g) = tempR(tempRE);
                    geneAvgExp(g) = avgExpVals(tempR(tempRE));
                    geneExpRank(g) = sortingInd_avgExpVals(tempR(tempRE));
                    geneLocScore(g) = locScores(tempR(tempRE));
                    geneLocRank(g) = sortingInd_locScores(tempR(tempRE));
                else
                    currGeneSetInd(g) = 0;
                    geneCorr(g) = 0;
                    geneCorrRank(g) = 0;
                    geneAvgExp(g) = 0;
                    geneExpRank(g) = 0;
                    geneLocScore(g) = 0;
                    geneLocRank(g) = 0;
                end                
            else
                currGeneSetInd(g) = 0;
                geneCorr(g) = 0;
                geneCorrRank(g) = 0;
                geneAvgExp(g) = 0;
                geneExpRank(g) = 0;
                geneLocScore(g) = 0;
                geneLocRank(g) = 0;
            end
            clear tempR;
        end
        clear g;
        %%% clean the data!!!!!!!!!!!!
        remRows = find(currGeneSetInd == 0);
        currGeneSetInd(remRows) = [];
        geneCorr(remRows) = [];
        geneCorrRank(remRows) = [];
        geneAvgExp(remRows) = [];
        geneExpRank(remRows) = [];
        geneLocScore(remRows) = [];
        geneLocRank(remRows) = [];
        %%% (1) calculate the ranksum pval based on correlation
        tempCorrArr = rankedCorrVals;
        tempCorrArr(currGeneSetInd) = [];
        pValCorr(experiment,s) = ranksum(tempCorrArr, geneCorr);
        %%% (1') prepare data for boxplot 
        if s == 1
            labels = zeros(1,length(geneCorr))+s;
            geneCorrPlot = geneCorr;
            geneRankPlot = geneCorrRank;
        else
            labels = [labels zeros(1,length(geneCorr))+s];
            geneCorrPlot = [geneCorrPlot geneCorr];
            geneRankPlot = [geneRankPlot geneCorrRank];
        end
        %%% (2) calculate the sum of the ranks of the current gene set
        corrSumOfRanks = sum(-log10(geneCorrRank/nGenes));
%         %%% (3) plot a distribution of rank sums for random sets of genes
%         %%% of the same size as the gene set of interest
%         noRandSets = 10000;
%         for rs = 1 : noRandSets
%             clear inds; clear randCorrArr; clear tempTempCorrArr; 
%             inds = randperm(length(tempCorrArr), length(currGeneSetInd));
%             randCorrArr = tempCorrArr(inds);
%             tempTempCorrArr = tempCorrArr; tempTempCorrArr(inds) = [];
%             randPValCorr(rs) = ranksum(tempTempCorrArr, randCorrArr);
%             randCorrSumOfRanks(rs) = sum(-log10(inds/ nGenes));
%         end
%         f = figure('Visible', 'off'); 
%         set(gcf,'units','normalized','outerposition',[0 0 1 1]);
%         hold on
%         subplot(1,2,1); hist(randCorrSumOfRanks, 1000); grid on
%         line([corrSumOfRanks corrSumOfRanks], [0 50], 'linewidth', 3, 'color', 'red', 'linestyle', '--')
%         text(corrSumOfRanks+2, 20, num2str(corrSumOfRanks), 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 15, 'Rotation', 90)
%         set(gca, 'XLim', [0 100]);
%         set(gca, 'YLim', [0 60]);
%         xlabel('Sum of the ranks', 'FontWeight', 'bold', 'FontSize', 15)
%         ylabel('Frequency', 'FontWeight', 'bold', 'FontSize', 15)
%         title({[allGenes{gene_index(experiment)} ' ' gene_experimentNos{experiment} ' ' gene_experimentPlanes{experiment}], ...
%             structures{s}}, 'FontWeight', 'bold', 'FontSize', 15);
%         subplot(1,2,2); hist(-log10(randPValCorr), 1000); grid on
%         line([-log10(pValCorr(experiment,s)) -log10(pValCorr(experiment,s))], [0 80], 'linewidth', 3, 'color', 'red', 'linestyle', '--')
%         text(-log10(pValCorr(experiment,s))+1, 40, num2str(pValCorr(experiment,s)), 'Color', 'r', 'FontWeight', 'bold', 'FontSize', 15, 'Rotation', 90)
%         set(gca, 'XLim', [0 20]);
%         set(gca, 'YLim', [0 120]);
%         xlabel('-log_1_0(ranksum p-value)', 'FontWeight', 'bold', 'FontSize', 15)
%         ylabel('Frequency', 'FontWeight', 'bold', 'FontSize', 15)
%         title({[allGenes{gene_index(experiment)} ' ' gene_experimentNos{experiment} ' ' gene_experimentPlanes{experiment}], ...
%             structures{s}}, 'FontWeight', 'bold', 'FontSize', 15);
%         hold off
%         saveas(f, [outDir '/' newSheetNames{experiment} '_' expType '_' structures{s} '_randomSet.fig']);
%         saveas(f, [outDir '/' newSheetNames{experiment} '_' expType '_' structures{s} '_randomSet.jpg']);
    end
    clear offset; clear Workbook; clear sheet1; clear eActivesheetRange;
    clear ri; clear r11; clear r2; clear r22; clear r3;
    %%% correct for multiple testing
    pValCorr_corrected(experiment,:) = multtest(pValCorr(experiment,:),'method','BH');
    %%% plot boxplots of the correlations in different structures together
    %%% with the corresponding p-value
%     f = figure; 
%     subplot(2,1,1), notBoxPlot(geneCorrPlot, labels); grid on;
%     ylabel('Correlation with Esr1', 'FontWeight', 'bold', 'FontSize', 15)
%     set(gca, 'XTickLabel', structures, 'FontWeight', 'bold', 'xlim', [0 numel(structures)+1]);
%     rotateXLabels(gca(), 45)
% %     title([allGenes{gene_index(experiment)} ' ' gene_experimentNos{experiment} ' ' gene_experimentPlanes{experiment} ' - ' expType],...
% %         'FontWeight', 'bold', 'FontSize', 15);
%     subplot(2,1,2), bar(-1*log10(pValCorr_corrected)); grid on;
%     ylabel('- log_1_0 (p-value)', 'FontWeight', 'bold', 'FontSize', 15)
%     set(gca, 'XTickLabel', structures, 'FontWeight', 'bold', ...
%         'XTick', 1:numel(structures), 'xlim', [0 numel(structures)+1], ...
%         'ylim', [0 8]);
% %     rotateXLabels(gca(), 45)
% %     title('Mann-Whitney U-test', 'FontWeight', 'bold', 'FontSize', 15);
%     line([0 numel(structures)+1], [-log10(0.05) -log10(0.05)], 'LineWidth', 3, 'color', 'r')
%     hold off
%     saveas(f, [outDir '/' newSheetNames{experiment} '_' expType '_CorrelationAnalysis.fig']);
%     saveas(f, [outDir '/' newSheetNames{experiment} '_' expType '_CorrelationAnalysis.jpg']);
end
% %%% write the results to excel
% xlswrite(outFile, structures, 1, 'B1');
% xlswrite(outFile, pValCorr, 1, 'B2');
% xlswrite(outFile, newSheetNames', 1, 'A2');






