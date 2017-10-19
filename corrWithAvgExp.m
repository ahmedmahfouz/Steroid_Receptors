%%% 24 Sep 2013
%%% calculate the correlation between a gene of interest and a set of genes 
%%% based on the expression of the main structures. Plot the average
%%% expression profiles

function corrWithAvgExp(filesDirectory, resultsDirectory, geneOfInterest, structures, ...
    geneSetName, allGenes, allExpNumbers, allExpPlanes, strAvgExp, strSize)

% load the set of genes
[~, txt] = xlsread([filesDirectory 'gene sets/' geneSetName '.xlsx']);
geneSet = txt(2:end,1);
clear txt;
%%% create the outdirectory for this analysis
outDir = [resultsDirectory geneSetName '/'];
if ~exist(outDir, 'dir')
    mkdir(outDir);
end
% find the indicies for the structures of interest and extract avg
% expression data for these structures
[num txt] = xlsread([filesDirectory 'ARAontology.xls']);
strNames = txt(2:end,5); 
clear num; clear txt;
strInd = find(ismember(strNames, structures) == 1);
%%% plot the structures' sizes
currStrSize = strSize(strInd(2:end));
currStrNames = strNames(strInd(2:end));
f = figure, bar(currStrSize), grid on
set(gca, 'XLim', [0 size(currStrSize,2)+1]);
set(gca, 'XTickLabel', currStrNames, 'FontWeight', 'bold', 'XTick', 1:numel(currStrNames));
rotateXLabels(gca, 45)
ylabel('structure size (voxels)', 'FontSize', 15, 'FontWeight', 'bold');
saveas(f, [outDir 'structureSizes_BigStructures.fig']);
saveas(f, [outDir 'structureSizes_BigStructures.jpg']);
%%% find all experiments of the geneOfInterest
gene_index = find(strcmpi(allGenes, geneOfInterest) == 1);
gene_experimentNos = allExpNumbers(gene_index);
gene_experimentPlanes = allExpPlanes(gene_index);
% prepare figure legend and xls file row names
for i = 1 : length(gene_index)
    rowNames{i} = [allGenes{gene_index(i)} ' ' gene_experimentNos{i} ' ' gene_experimentPlanes{i}];
end
%%% retreive the average expression profiles for the gene of interest
avgExpArr = strAvgExp(strInd,gene_index);
% [num txt] = xlsread([resultsDirectory 'strAvgExp/subStrAvgExp_' geneOfInterest{1} '_grey.xls']);
% strNames = txt(2:end,1);
% avgExpArr = num(2:end,3:2+length(gene_index));
%%% Plot the expression levels of the gene of ineterst
f = figure('Visible', 'off');
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
bar(avgExpArr); grid on
set(gca, 'XLim', [0 size(avgExpArr,1)+1]);
set(gca, 'XTickLabel', structures, 'FontWeight', 'bold', 'XTick', 1:numel(structures));
rotateXLabels(gca, 45)
ylabel('Average Expression Energy', 'FontSize', 15, 'FontWeight', 'bold');
legend(rowNames)
title(geneOfInterest, 'FontSize', 15, 'FontWeight', 'bold');
saveas(f, [outDir geneOfInterest{1} '_avgExp.fig']);
saveas(f, [outDir geneOfInterest{1} '_avgExp.jpg']);

%%% retreive the average expression profiles for the structures of interest
%%% for each gene
indecies = 0;
for i = 1 : length(geneSet)
    currGene_index = find(strcmpi(allGenes, geneSet{i}) == 1);
    currGene_experimentNos = allExpNumbers(currGene_index);
    currGene_experimentPlanes = allExpPlanes(currGene_index);
    indecies = indecies(end)+1 : indecies(end)+length(currGene_index);
    geneNames(indecies) = allGenes(currGene_index);
    geneExpNos(indecies) = currGene_experimentNos;
    geneExpPlanes(indecies) = currGene_experimentPlanes;
    avgExpMat(:,indecies) = strAvgExp(strInd,currGene_index);
%     [num txt] = xlsread([resultsDirectory 'strAvgExp/subStrAvgExp_' geneSet{i} '_grey.xls']);
%     avgExpMat(:,indecies) = num(2:end,3:2+length(currGene_index));
    %%% prepare figure legend and xls file col names    
    for j = 1 : length(currGene_index);
        colNames{indecies(j)} = [geneNames{indecies(j)} ' ' geneExpNos{indecies(j)} ' ' geneExpPlanes{indecies(j)}];
    end
    %%% Plot the expression levels of the gene of ineterst
    f = figure('Visible', 'off');
    set(gcf,'units','normalized','outerposition',[0 0 1 1]);
    bar(strAvgExp(strInd,currGene_index)); grid on
    set(gca, 'XLim', [0 size(avgExpArr,1)+1]);
    set(gca, 'XTickLabel', structures, 'FontWeight', 'bold', 'XTick', 1:numel(structures));
    rotateXLabels(gca, 45)
    ylabel('Average Expression Energy', 'FontSize', 15, 'FontWeight', 'bold');
    legend(colNames{indecies})
    title(geneSet{i}, 'FontSize', 15, 'FontWeight', 'bold');
    saveas(f, [outDir geneSet{i} '_avgExp.fig']);
    saveas(f, [outDir geneSet{i} '_avgExp.jpg']);
end
%%% plot the average expressions 
avgExpMat = [avgExpArr avgExpMat];
% avgExpMat = avgExpMat ./ repmat(max(avgExpMat), size(avgExpMat,1), 1);
allGeneNames = [rowNames colNames];
f = figure;
imagesc(avgExpMat'); colorbar
set(gcf,'units','normalized','outerposition',[0 0 1 1]);
set(gca, 'XTickLabel', structures, 'FontWeight', 'bold', 'XTick', 1:numel(structures));
rotateXLabels(gca, 45)
set(gca, 'YTickLabel', allGeneNames, 'FontWeight', 'bold', 'YTick', 1:numel(allGeneNames));
title('Average Expression - HY Substructures', 'FontSize', 15, 'FontWeight', 'bold');
saveas(f, [outDir 'avgExp_HYsubstructures.fig']);
saveas(f, [outDir 'avgExp_HYsubstructures.jpg']);
%%% calculate the correlation matrix
corrMat = corr(avgExpArr, avgExpMat);
save([outDir 'corrAvgExp_Esr1_' geneSetName '.mat'], 'corrMat');
outFile = [outDir 'corrAvgExp_Esr1_' geneSetName '.xlsx'];
%%% save to xls
xlswrite(outFile, rowNames', 1, 'A2');
xlswrite(outFile, colNames, 1, 'B1');
xlswrite(outFile, corrMat, 1, 'B2');









