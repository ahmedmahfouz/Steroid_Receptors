%%% a function to read an excel file with a correlation list (output of
%%% xlsCorrRes.m) and return data in a struct

function output = readCorrFile(xlsDIR, geneOfInterest, structure, ...
        noOfGenes, extension, experiment);
    
% read the correlation file
[num txt] = xlsread([xlsDIR geneOfInterest{1} '/' geneOfInterest{1} '_' structure extension '.xls'], experiment);
if ~isempty(num)
    output.seedGene = txt{1,1};
    output.strOfInterest = txt{1,2};
    output.rank = num(:,1);
    output.rankedGenes = txt(3:end, 2);
    output.expNum = num(:,3);
    output.expPlane = txt(3:end,4);
    output.corrVals = num(:,5);
    output.pVals = num(:,6);
    output.avgExp = num(:,7);
    output.norAvgExp = num(:,8);
    output.noVox = num(:,10);
    output.expRank = num(:,11);
    output.relExpRank = num(:,12);
    output.relCorrRank = num(:,13);
    output.rankProd = num(:,14);
    % select the top n genes
    output.topNgenesInd = 1:noOfGenes;
    output.topNgenes = output.rankedGenes(1:noOfGenes);
    % select the bottom n genes
    lowestCorrInd = find(output.corrVals > 0, 1, 'last');
    ind1 = lowestCorrInd-noOfGenes+1;
    if ind1 <= 0
        ind1 = 1;
    end
    output.bottomNGenesInd = ind1:lowestCorrInd;
    output.bottomNgenes = output.rankedGenes(ind1:lowestCorrInd);  
else
    output = [];
end




