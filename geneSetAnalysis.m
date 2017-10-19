%%% 24 May 2013
%%% 

function [geneSet geneCorr geneRank] = geneSetAnalysis(filesDirectory, resultsDirectory, geneOfInterest, structures, ...
    geneSetName, expType);

[~, txt] = xlsread([filesDirectory 'gene sets/' geneSetName '.xlsx']);
geneSet = txt(2:end,1);
clear txt;

if strcmp(expType, 'C')
    %%% loop on all the structures in the structures list
    for s = 1 : length(structures)
        [num txt] = xlsread([resultsDirectory 'xlsFiles/CoronalOnly/' geneOfInterest '/' geneOfInterest '_' structures{s} '_CoronalOnly.xls']);
        rankedGenes = txt(3:end, 2);
        rankedCorrVals = num(1:end,3);
        rankedExpNo = num(1:end,5);
        clear txt; clear num;
        % return only the highest rank for each gene (if the gene has more than 
        % one experiment)
        for g = 1 : length(geneSet)
            tempR = find(strcmpi(rankedGenes, geneSet{g}) == 1);
            if ~isempty(tempR)
                geneRank(g,s) = min(tempR);
                geneCorr(g,s) = rankedCorrVals(min(tempR));
            else
                geneRank(g,s) = 0;
                geneCorr(g,s) = 0.5;
            end
            clear tempR;
        end
        clear g;
    end
    
elseif strcmp(expType, 'S')
    %%% loop on all the structures in the structures list
    for s = 1 : length(structures)
        [num txt] = xlsread([resultsDirectory 'xlsFiles/SagittalOnly/' geneOfInterest '/' geneOfInterest{NR} '_' structures{s} '_SagittalOnly.xls']);
        rankedGenes = txt(3:end, 2);
        rankedCorrVals = num(1:end,3);
        rankedExpNo = num(1:end,5);
        clear txt; clear num;
        % return only the highest rank for each gene (if the gene has more than 
        % one experiment)
        for g = 1 : length(geneSet)
            tempR = find(strcmpi(rankedGenes, geneSet{g}) == 1);
            if ~isempty(tempR)
                geneRank(g,s) = min(tempR);
                geneCorr(g,s) = rankedCorrVals(min(tempR));
            else
                geneRank(g,s) = 0;
                geneCorr(g,s) = 0.5;
            end
            clear tempR;
        end
        clear g;
    end
    
elseif strcmp(expType, 'All')
    %%% determine how many probes are there per NR
    load([filesDirectory 'allGenes.mat']);
    gene_index = find(strcmpi(allGenes, geneOfInterest) == 1);
    load([filesDirectory 'allExpNumbers.mat']);
    load([filesDirectory 'allExpPlanes.mat']);
    nrExpNo = allExpNumbers(gene_index);
    nrExpP = allExpPlanes(gene_index);
    clear allGenes; clear allExpNumbers; clear allExpPlanes;
    for i = 1 : length(gene_index)
        %%% loop on all the structures in the structures list
        for s = 1 : length(structures)
            [num txt] = xlsread([resultsDirectory 'xlsFiles/All/' geneOfInterest '/' geneOfInterest '_' structures{s} '.xls'], i);
            if ~isempty(txt)
                rankedGenes = txt(3:end, 2);
                rankedCorrVals = num(1:end,3);
                rankedExpNo = num(1:end,5);
                clear txt; clear num;
                % return only the highest rank for each gene (if the gene has more than 
                % one experiment)
                for g = 1 : length(geneSet)
                    tempR = find(strcmpi(rankedGenes, geneSet{g}) == 1);
                    if ~isempty(tempR)
                        geneRank(g,s,i) = min(tempR);
                        geneCorr(g,s,i) = rankedCorrVals(min(tempR));
                    else
                        geneRank(g,s,i) = 0;
                        geneCorr(g,s,i) = 0.5;
                    end
                    clear tempR;
                end
                clear g;
            end
        end
    end
    
else
    display('expType not defined correctly');
end

% %%% load the GR-coregulators file to plot the interaction-prediction scores
% [num txt] = xlsread([filesDirectory 'GR-coregulatorsIP.xlsx']);
% origAvgDMSO = num(:,1);
% origAvgDEX = num(:,2);
% origIPScore = num(:,3); %abs(DEX-DMSO)
% originalCoregList = txt(2:end, 2);
% clear num; clear txt;

%%% return the size of the geneCorr/geneRank matrix
SIZE = size(geneCorr);
if length(SIZE) > 2 % a gene has more than one experiment
    for R = 1 : SIZE(3)
        currGeneRank = geneRank(:,:,R);
        currGeneCorr = geneCorr(:,:,R);
        currGeneSet = geneSet;
        % remove rows/genes with no ranking/correlations (missing experiments)
        gToR = find(sum(currGeneRank, 2) == 0);
        currGeneRank(gToR,:) = [];
        currGeneCorr(gToR,:) = [];
        currGeneSet(gToR) = [];

%         corrType = 'correlation';
%         linkType = 'average'; % linkage type
%         C1 = clustergram(currGeneCorr, 'Standardize', 'none', 'Cluster', 1, 'RowPDist', ...
%             corrType, 'Linkage', linkType, 'Colormap', 'redbluecmap', ...
%             'RowLabels', currGeneSet, 'ColumnLabels',structures);
%         h1 = plot(C1);
%         saveas(h1, ['results\plots\' geneOfInterest '_' nrExpNo{R} '_' nrExpP{R} '_' geneSetName '_' expType '_corr.fig']);
    %     saveas(h1, ['results\plots\' geneOfInterest '_' nrExpNo{R} '_' nrExpP{R} '_' geneSetName '_' expType '_corr.jpg']);
        
%         [~, permrows] = ismember(get(C1, 'RowLabels'), currGeneSet);
%         currIPS = origIPScore;
%         currIPS(gToR) = [];
%         f = figure; barh([1:length(currIPS)], currIPS(permrows)); grid on
%         set(gca, 'YLim', [1 length(currIPS)+1]);
%         set(gca, 'YTickLabel', currGeneSet(permrows), 'YTick', 1:numel(currGeneSet));
%         saveas(f, ['results\plots\' geneOfInterest '_' nrExpNo{R} '_' nrExpP{R} '_' geneSetName '_' expType '_ips.fig']);
        
%         oF = ['results\plots\' geneOfInterest '_' nrExpNo{R} '_' nrExpP{R} '_' geneSetName '_' expType '_corr.xls'];
%         xlswrite(oF, currGeneSet(permrows), 1, 'A2');
%         xlswrite(oF, structures, 1, 'B1');
%         xlswrite(oF, currGeneCorr(permrows,:), 1, 'B2');
% %         xlswrite(oF, {'|DEX-DMSO|'}, 1, 'AL1');
% %         xlswrite(oF, currIPS(permrows), 1, 'AL2');
%         
%         xlswrite(oF, currGeneSet(permrows), 2, 'A2');
%         xlswrite(oF, structures, 2, 'B1');
%         xlswrite(oF, currGeneRank(permrows,:), 2, 'B2');
% %         xlswrite(oF, {'|DEX-DMSO|'}, 2, 'AL1');
% %         xlswrite(oF, currIPS(permrows), 2, 'AL2');
        
        
%         C2 = clustergram(currGeneRank, 'Standardize', 'none', 'Cluster', 1, 'RowPDist', ...
%             corrType, 'Linkage', linkType, 'Colormap', 'redbluecmap', ...
%             'RowLabels', currGeneSet, 'ColumnLabels',structures);
%         h2 = plot(C2);
%         saveas(h2, ['results\plots\' geneOfInterest '_' nrExpNo{R} '_' nrExpP{R} '_' geneSetName '_' expType '_rank.fig']);
%     %     saveas(h2, ['results\plots\' geneOfInterest '_' nrExpNo{R} '_' nrExpP{R} '_' geneSetName '_' expType '_rank.jpg']);
    end
    
else
    % remove rows/genes with no ranking/correlations (missing experiments)
    gToR = find(sum(geneRank, 2) == 0);
    geneRank(gToR,:) = [];
    geneCorr(gToR,:) = [];
    geneSet(gToR) = [];
    
    % show a heat map of the corr & rank of the geneSet
    % figure, imagesc(geneRank);
    % figure, imagesc(geneCorr);

%     corrType = 'correlation';
%     linkType = 'average'; % linkage type
%     C1 = clustergram(geneCorr, 'Standardize', 'none', 'Cluster', 1, 'RowPDist', ...
%         corrType, 'Linkage', linkType, 'Colormap', 'redbluecmap', ...
%         'RowLabels', geneSet, 'ColumnLabels',structures);
%     C2 = clustergram(geneRank, 'Standardize', 'none', 'Cluster', 1, 'RowPDist', ...
%         corrType, 'Linkage', linkType, 'Colormap', 'redbluecmap', ...
%         'RowLabels', geneSet, 'ColumnLabels',structures);
% 
%     h1 = plot(C1);
%     saveas(h1, ['results\plots\' geneOfInterest '_' geneSetName '_' expType '_corr.fig']);
% %     saveas(h1, ['results\plots\' geneOfInterest '_' geneSetName '_' expType '_corr.jpg']);
% 
%     h2 = plot(C2);
%     saveas(h2, ['results\plots\' geneOfInterest '_' geneSetName '_' expType '_rank.fig']);
% %     saveas(h2, ['results\plots\' geneOfInterest '_' geneSetName '_' expType '_rank.jpg']);
end








