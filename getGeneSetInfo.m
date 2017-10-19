%%% 2 October 2015
%%% get gene set info from a correlation list

%% read lists of sorted co-expressed genes
geneName = 'Esr1';
structures = {'grey', 'HY'};
geneSetName = 'xuCell2012';
sheetNo = 1; % sheet corresponding to Esr1 coronal experiment #79591677
if ispc
    dataDir = 'E:\Ahmed\HP\work\Results\NuclearReceptors\xlsFiles\All\Esr1\';
    filesDirectory = 'E:\Ahmed\HP\work\Data\NuclearReceptors\files\';
    resultsDir = 'E:\Ahmed\HP\work\Results\NuclearReceptors\results_esr1_xuetal_2Oct2015\';
end

% read gene set
load([filesDirectory 'allGenes.mat']);
load([filesDirectory 'allExpNumbers.mat']);
load([filesDirectory 'allExpPlanes.mat']);
[~, txt] = xlsread([filesDirectory 'gene sets/' geneSetName '.xlsx']);
geneSet = txt(2:end,1);
clear txt;
geneSetInd = find(ismember(allGenes, geneSet) == 1);
geneSetNames = geneSet';
% geneSetNames = allGenes(geneSetInd);
% geneSetExpNo = allExpNumbers(geneSetInd);
% geneSetExpPlanes = allExpPlanes(geneSetInd);

for i = 1 : length(structures)
    T = readtable([dataDir geneName '_' structures{i} '.xls'],...
        'Sheet',sheetNo,'Range','A2:N26020');
    expRank = tiedrank(T.avg_Exp_);
    expRank = max(expRank) - expRank + 1;
    T.exp_Rank = max(T.exp_Rank) - T.exp_Rank;
    for G = 1 : length(geneSetNames)
        temp = find(strcmpi(T.geneSymbol,geneSetNames{G})==1,1);
        if ~isempty(temp)
            I1(G) = temp;
            I2(G) = T.corr_(temp);
            I_exp(G) = T.norm_Avg_Exp_(temp);
            I_exp_rank(G) = T.exp_Rank(temp);
            I_avgExp(G) = T.avg_Exp_(temp);
            I_avgExp_rank(G) = expRank(temp);
        else
            I1(G) = NaN;
            I2(G) = NaN;
            I_exp(G) = NaN;
            I_exp_rank(G) = NaN;
            I_avgExp(G) = NaN;
            I_avgExp_rank(G) = NaN;
        end
    end
    ranks{i,:} = I1;
    corr{i,:} = I2;
    exp{i,:} = I_exp;
    exp_rank{i,:} = I_exp_rank;
    avg_exp{i,:} = I_avgExp;
    avg_exp_rank{i,:} = I_avgExp_rank;
    clear I1; clear I2; clear I_exp; clear I_exp_rank; clear I_avgExp; clear I_avgExp_rank
end

%% save the results
for i = 1 : length(structures)
    T = table(geneSetNames', ranks{i,:}', corr{i,:}', exp{i,:}', exp_rank{i,:}',...
        avg_exp{i,:}',avg_exp_rank{i,:}',...
        'VariableNames',{'gene_symbol', 'correlation_rank','correlation',...
        'normaalized_average_expression','normalized_expression_rank','average_expression','average_expression_rank'});
    writetable(T, [resultsDir geneName 'coexpression.xlsx'],...
        'Sheet',i);
end



