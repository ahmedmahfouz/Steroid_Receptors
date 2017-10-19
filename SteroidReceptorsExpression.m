%%% 11 March 2015
%%% Steroid Receptors expression analysis 

%% define the files and results folders
if ispc
    filesDirectory = 'E:\Ahmed\HP\work\Data\NuclearReceptors\';
    resultsDirectory = 'E:\Ahmed\HP\work\Results\NuclearReceptors\';
else
    filesDirectory = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Data/NuclearReceptors/files/';
    resultsDirectory = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/NuclearReceptors/';
end
% define the NRs to analyze
geneOfInterest = {'Esr1', 'Esr2', 'Ar', 'Nr3c1', 'Nr3c2', 'Pgr'};
%%% Main brain structures 
mainStr = {'grey', 'Isocortex', 'OLF', 'HPF' 'CTXsp', 'STR', 'PAL', 'CB', ...
    'TH', 'HY', 'MB', 'P', 'MY'};
HIP = {'CA1', 'CA2', 'CA3', 'DG'};
structures = [mainStr HIP];
% structuresFullName = {'Cortex', 'Olfactory areas', 'Hippocampal formation' 'Cortical supplate',...
%     'Striatum', 'Pallidum', 'Cerebellum', 'Thalamus', 'Hypothalamus', 'Midbrain',...
%     'Pons', 'Medulla'};
% structuresFullName = {'Ventral tegmental area', 'Substantia nigra_reticular part', 'Substantia nigra_compact part'};
% type of experiments to use (All or C)
expType = 'All';

load('geneOfInterest.mat')

%% load the expression data and calculate the average expressions
if strcmp(expType, 'C')
    % load gene names
    load([filesDirectory 'listOfGenes.mat']);
    % load experiments' numbers
    load([filesDirectory 'listOfExperimentNumbers.mat']);
    % load experiments' plane [useless in case of 'coronal only']
    load([filesDirectory 'allExpPlanes.mat']);
    % load voxels annotations
    load([filesDirectory 'voxelAnnotation_CoronalOnly.mat']);
    % load brain voxels 
    load([filesDirectory 'brainInd_CoronalOnly.mat']);
    % load expression matrix
    load([filesDirectory 'expressionMatrix_CoronalOnly.mat']);
    expressionMatrix_CoronalOnly(expressionMatrix_CoronalOnly == -1) = NaN;
%     expressionMatrix_CoronalOnly = (expressionMatrix_CoronalOnly - repmat(nanmean(expressionMatrix_CoronalOnly,2),1,size(expressionMatrix_CoronalOnly,2)))./repmat(nanstd(expressionMatrix_CoronalOnly,[],2),1,size(expressionMatrix_CoronalOnly,2));
  for S = 1 : length(structures)
        str_voxels = indicateSubStr(structures{S}, filesDirectory);
        s_idx = find(ismember(voxelAnnotation_CoronalOnly,str_voxels)==1);
        for i = 1 : length(geneOfInterest)
            % check wether this gene has a coronal experiment or not
            geneExpCheck = find(strcmpi(listOfGenes, geneOfInterest{i}) == 1);
            if isempty(geneExpCheck)
    %                 display(['"' geneOfInterest{i} '" has no coronal experiment']);
                avgExp(i,S) = NaN;
            else
                expMatPerStr{i,S} = expressionMatrix_CoronalOnly(s_idx,geneExpCheck)';
                avgExp(i,S) = nanmean(expressionMatrix_CoronalOnly(s_idx,geneExpCheck));
            end
        end
    end
%     save([resultsDirectory 'structures_expression_NoNormalization_C.mat'],'avgExp','expMatPerStr');
    save([resultsDirectory 'structures_expression_NoNormalization_C_dopa_coreg.mat'],'avgExp','expMatPerStr');
elseif strcmp(expType, 'All')
    % load gene names
    load([filesDirectory 'allGenes.mat']);
    % load experiments' numbers
    load([filesDirectory 'allExpNumbers.mat']);
    % load voxels annotations
    load([filesDirectory 'voxelAnnotation.mat']);
    % load brain voxels 
    load([filesDirectory 'brainInd.mat']);
    % load expression matrix
    load([filesDirectory 'fullExpressionMatrix.mat']);
    fullExpressionMatrix(fullExpressionMatrix == -1) = NaN;
    % check wether this gene has a coronal experiment or not
    geneExpCheck = find(ismember(allGenes(4346:end), geneOfInterest) == 1)+4345;
    genes = allGenes(geneExpCheck);
    for S = 1 : length(structures)
        str_voxels = indicateSubStr(structures{S}, filesDirectory);
        s_idx = find(ismember(voxelAnnotation,str_voxels)==1);
        for i = 1 : length(geneExpCheck)
            expMatPerStr{i,S} = fullExpressionMatrix(s_idx,geneExpCheck(i))';
            avgExp(i,S) = nanmean(fullExpressionMatrix(s_idx,geneExpCheck(i)));
        end
    end
%     save([resultsDirectory 'structures_expression_NoNormalization_All.mat'],'avgExp','expMatPerStr','genes');
    save([resultsDirectory 'structures_expression_NoNormalization_All_dopa_coreg.mat'],'avgExp','expMatPerStr','genes');
end

%%
load([filesDirectory 'structures_expression_NoNormalization.mat'])

figure, 
B = bar(avgExp(1:5,:),'EdgeColor','w')
C = jet(13); colormap(C);
set(gca, 'XTick', 1:5, 'XTickLabel', geneOfInterest(1:5),'FontWeight', 'bold', 'FontSize', 15)
ylabel('Average Expression Energy','FontWeight', 'bold', 'FontSize', 15)
legend(structures)

% figure, imagesc(log2(avgExp(1:5,:)+1)), colormap('redbluecmap')
% set(gca, 'XTick', 1:length(structures), 'XTickLabel', structures)
% set(gca, 'YTick', 1:length(geneOfInterest)-1, 'YTickLabel', geneOfInterest(1:5))
% legend on

%%
load([filesDirectory 'structures_expression_NoNormalization_All.mat'])

figure, 
B = bar(avgExp,'EdgeColor','w');
C = jet(13); colormap(C);
set(gca, 'XTickLabel', genes,'FontWeight', 'bold', 'FontSize', 15)
ylabel('Average Expression Energy','FontWeight', 'bold', 'FontSize', 15)
legend(structures)

%% data
geneOfInterest = {'Ar', 'Pgr', 'Esr1', 'Esr2', 'Nr3c1', 'Nr3c2'};
% AR, Esr1 and Esr2: coronal
load([filesDirectory 'structures_expression_NoNormalization.mat'])
dat([1,3,4],:) = avgExp([1,3,4],:);
% GR, MR, and Pgr: sagittal
load([filesDirectory 'structures_expression_NoNormalization_All.mat'])
dat([6,5,2],:) = avgExp([4,6,8],:);
dat = dat ./ repmat(dat(:,1),1,size(dat,2));

figure; imagesc(log2(dat(:,2:end))), colormap('redbluecmap')
colorbar

set(gca, 'YTickLabel', geneOfInterest, 'YTick', 1:size(dat,1), 'FontWeight', 'bold', 'FontSize', 15)
set(gca, 'XTickLabel', structuresFullName, 'XTick', 1:size(dat,2)-1, 'FontWeight', 'bold', 'FontSize', 15)

rotateXLabels(gca(),45)
xticklabel_rotate(1:size(dat,2)-1,45,structuresFullName,'FontWeight', 'bold', 'FontSize', 15)

figure, 
B = bar(dat,'EdgeColor','w');
C = jet(13); colormap(C);
set(gca, 'XTickLabel', geneOfInterest,'FontWeight', 'bold', 'FontSize', 15)
ylabel('Average Expression Energy','FontWeight', 'bold', 'FontSize', 15)
legend(structures)
grid on

%% VTA_SN
geneOfInterest = {'Ar', 'Pgr', 'Esr1', 'Esr2', 'Nr3c1', 'Nr3c2'};
% AR, Esr1 and Esr2: coronal
load([filesDirectory 'structures_expression_NoNormalization_C_dopa_receptors.mat'])
dat([1,3,4],:) = avgExp([1,3,4],:);
% GR, MR, and Pgr: sagittal
load([filesDirectory 'structures_expression_NoNormalization_All_dopa_receptors.mat'])
dat([6,5,2],:) = avgExp([4,6,8],:);
% dat = dat ./ repmat(dat(:,1),1,size(dat,2));

figure; imagesc(log2(dat)), colormap('redbluecmap')
colorbar

set(gca, 'YTickLabel', geneOfInterest, 'YTick', 1:size(dat,1), 'FontWeight', 'bold', 'FontSize', 15)
set(gca, 'XTickLabel', structuresFullName, 'XTick', 1:size(dat,2)-1, 'FontWeight', 'bold', 'FontSize', 15)

rotateXLabels(gca(),45)
xticklabel_rotate(1:size(dat,2)-1,45,structuresFullName,'FontWeight', 'bold', 'FontSize', 15)

figure, 
B = bar(dat,'EdgeColor','w');
C = jet(13); colormap(C);
set(gca, 'XTickLabel', geneOfInterest,'FontWeight', 'bold', 'FontSize', 15)
ylabel('Average Expression Energy','FontWeight', 'bold', 'FontSize', 15)
legend(structures)
grid on

%%
load([filesDirectory 'structures_expression_NoNormalization_All_coreg.mat'])


for i = 1 : length(geneOfInterest)
    idx = find(strcmpi(genes,geneOfInterest{i})==1);
    if length(idx) == 1
        dat(i,:) = avgExp(idx,:);
    else
        dat(i,:) = nanmean(avgExp(idx,:));
    end
end

dat = dat ./ repmat(dat(:,1),1,size(dat,2));

figure; imagesc(log2(dat(:,2:end))), colormap('redbluecmap')
colorbar
% set(gca, 'XTick', 1:5, 'XTickLabel', structuresFullName,'FontWeight', 'bold', 'FontSize', 15)
set(gca, 'YTick', 1:length(geneOfInterest), 'YTickLabel', geneOfInterest)

% figure, 
% B = bar(avgExp(1:5,:),'EdgeColor','w')
% C = jet(13); colormap(C);
% set(gca, 'XTick', 1:5, 'XTickLabel', geneOfInterest(1:5),'FontWeight', 'bold', 'FontSize', 15)
% ylabel('Average Expression Energy','FontWeight', 'bold', 'FontSize', 15)
% legend(structures)
