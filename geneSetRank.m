%%% 5 March 2013
%%% Assign a rank-measure for an input set of genes, based on the
%%% correlation to a seed gene in different structures.
%%% a p-value is assigned to each rank-measure based on a permutation test
%%% (repeating the test 10000 times for a random set of genes of the same
%%% size as the input set

function geneSetRank(seedGene, geneSetName, listOfStructures, TYPE)

filesDirectory = 'files/';
resultsDirectory = 'results/';
resFolder = [resultsDirectory 'xlsFiles/'];
outDir = [resultsDirectory 'R/' geneSetName '_' TYPE '/'];
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% how many probes does the seed gene have?
load([filesDirectory 'allGenes.mat']);
load([filesDirectory 'allExpNumbers.mat']);
gene_index = find(strcmp(allGenes, seedGene) == 1);
clear allGenes;

% read the gene set 
[~, txt] = xlsread([filesDirectory 'gene sets/' geneSetName '.xlsx']);
geneSet = txt(2:end,1);
clear txt;

for p = 1 : length(gene_index)
%     clear rankedGenes_brain; clear rankedCorrVals_brain; clear rankedExpNo_brain;
%     clear rankedPlanes_brain; clear geneRank_brain; 
%     % calculate the gene set rank measure for the correlation list
%     % corresponding to the whole brain
%     [num txt] = xlsread([resFolder seedGene '_grey.xls'], p);
%     rankedGenes_brain = txt(4:end, 2);
%     rankedCorrVals_brain = num(2:end,3);
%     rankedExpNo_brain = num(2:end,4);
%     rankedPlanes_brain = txt(4:end, 5);
%     clear txt; clear num;
%     % return only the highest rank for each gene (if the gene has more than 
%     % one experiment)
%     for g = 1 : length(geneSet)
%         tempR = find(strcmp(rankedGenes_brain, geneSet{g}) == 1);
%         geneRank_brain(g) = min(tempR);
%         clear tempR;
%     end
%     clear g;
    % calculate the gene set rank measure for each correlation list (in
    % different structures)
    clear rankedGenes; clear rankedCorrVals; clear rankedExpNo;
    clear rankedPlanes; clear geneRank; 
    for s = 1 : length(listOfStructures)
        [num txt] = xlsread([resFolder seedGene '/' seedGene '_' listOfStructures{s} '.xls'], p);
        rankedGenes = txt(3:end, 2);
        rankedCorrVals = num(1:end,3);
        rankedExpNo = num(1:end,5);
        rankedPlanes = txt(3:end,6);
        clear txt; clear num;
        % return only the highest rank for each gene (if the gene has more than 
        % one experiment)
        for g = 1 : length(geneSet)
            tempR = find(strcmp(rankedGenes, geneSet{g}) == 1);
            if ~isempty(tempR)
                geneRank(g,s) = min(tempR);
            else
                geneRank(g,s) = 0;
            end
            clear tempR;
        end
        clear g;
    %     % Measure#1: Wilcoxon signed rank test
    %     pVals(s) = signrank(geneRank_brain, geneRank(:,s));
    %     % Measure#2: Rank Product
    %     rankProduct(s) = prod(geneRank(:,s)) / length(geneSet);
        % Measure#3: Mean-Rang Gene Test (prepare for the test in R)
        xlswrite([outDir seedGene '_'  allExpNumbers{gene_index(p)} '_' geneSetName '.xls'], {listOfStructures{s}}, s, 'A1');
        xlswrite([outDir seedGene '_'  allExpNumbers{gene_index(p)} '_' geneSetName '.xls'], {'gene ranks'}, s, 'A2');
        xlswrite([outDir seedGene '_'  allExpNumbers{gene_index(p)} '_' geneSetName '.xls'], geneRank(:,s), s, 'A3');
        xlswrite([outDir seedGene '_'  allExpNumbers{gene_index(p)} '_' geneSetName '.xls'], {'corr'}, s, 'B2');
        xlswrite([outDir seedGene '_'  allExpNumbers{gene_index(p)} '_' geneSetName '.xls'], rankedCorrVals, s, 'B3');
    end
end

    




    
    

