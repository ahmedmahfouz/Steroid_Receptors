%%%  Dec 2014
%%% PNAS Figures

%% Figure 1A
filesDirectory = 'E:\Ahmed\HP\work\Data\NuclearReceptors\';
geneOfInterest = {'Ar', 'Pgr', 'Esr1', 'Esr2', 'Nr3c1', 'Nr3c2'};
% AR, Esr1 and Esr2: coronal
load([filesDirectory 'structures_expression_NoNormalization.mat'])
dat([1,3,4],:) = avgExp([1,3,4],:);
% GR, MR, and Pgr: sagittal
load([filesDirectory 'structures_expression_NoNormalization_All.mat'])
dat([6,5,2],:) = avgExp([4,6,8],:);
dat = dat ./ repmat(dat(:,1),1,size(dat,2));
figure; imagesc(log2(dat(:,2:end))), colormap('redbluecmap')
X = dat(:,2:end);
hold on;
for i = 1:size(X,1)
   plot([.5,size(X,2)+.5],[i-.5,i-.5],'k-');
end
for i = 1:size(X,2)
    plot([i-.5,i-.5],[.5,size(X,1)+.5],'k-');
end
hold off
axis image %// equal scale on both axes

%% Figure 2A
resultsDirectory = 'E:\Ahmed\HP\work\Results\NuclearReceptors\';
load([resultsDirectory 'CORT_results_July2015.mat']);
fileExtensions = {'GC-resposive genes', 'GC-resposive genes in control rats',...
    'GC-resposive genes in CRS rats', 'common GC-resposive genes between control & CRS rats'};
structures = {'Whole Brain', 'Isocortex', 'HY', 'Hippocampus', 'CA', 'DG', 'CA1', 'CA2', 'CA3'};
% remove HY, CTX and CA
pFinal(:,[2,3,5]) = [];
structures([2,3,5]) = [];

f = figure,
set(f,'Position',[200, 200, 1024, 512])
b = bar(-log10(pFinal));
b(1).FaceColor = [0.4 0.4 0.4]; b(1).EdgeColor = [0.4 0.4 0.4];
b(2).FaceColor = [0 0 153/255]; b(2).EdgeColor = [0 0 153/255];
b(3).FaceColor = [1 69/255 0]; b(3).EdgeColor = [1 69/255 0];
b(4).FaceColor = [0 100/255 0]; b(4).EdgeColor = [0 100/255 0];
b(5).FaceColor = [50/255 205/255 50/255]; b(5).EdgeColor = [50/255 205/255 50/255];
b(6).FaceColor = [60/255 169/255 113/255]; b(6).EdgeColor = [60/255 169/255 113/255];

grid on
grid minor
legend(structures)
line([0.5 4.5],[-log10(0.05) -log10(0.05)],'Color','black','LineWidth',2,'LineStyle','--');

%% Figure 2B
resultsDirectory = 'E:\Ahmed\HP\work\Results\NuclearReceptors\';
load([resultsDirectory 'CORT_results_July2015_allSTR.mat']);
fileExtensions = {'GC-resposive genes (n = 497)', 'GC-resposive genes in control rats (n = 147)',...
    'GC-resposive genes in CRS rats (n = 182)', 'Common GC-resposive genes between control & CRS rats (n = 168)'};

f = figure,
set(f,'Position',[200, 200, 1024, 512])
b = bar(-log10(pFinal)');
b(1).FaceColor = [178/255 34/255 34/255]; b(1).EdgeColor = [178/255 34/255 34/255];
b(2).FaceColor = [70/255 130/255 180/255]; b(2).EdgeColor = [70/255 130/255 180/255];
b(3).FaceColor = [1 215/255 0]; b(3).EdgeColor = [1 215/255 0];
b(4).FaceColor = [50/255 205/255 50/255]; b(4).EdgeColor = [50/255 205/255 50/255];

grid on
grid minor
legend(fileExtensions)
line([0 14],[-log10(0.05) -log10(0.05)],'Color','black','LineWidth',2,'LineStyle','--');

%% Figure 3A
resultsDirectory = 'E:\Ahmed\HP\work\Results\NuclearReceptors\xuCell2012\';
pValues = xlsread([resultsDirectory 'xuCell2012_Esr1rank-pVal_All.xls'],1,'B2:N2');
HEX = {'70FF71', '7ED04B', '9AD2BD', '8ADA87', '98D6F9', '8599CC', ...
    'F0F080', 'FF7080', 'E64438', 'FF64FF', 'FF9B88', 'FF9BCD'};

f = figure, hold on
set(f,'Position',[200, 200, 1024, 512])
b = bar(1,-log10(pValues(1)),'FaceColor',[0.4 0.4 0.4],'EdgeColor',[0.4 0.4 0.4]);
for i = 2 : length(pValues)
    C = hex2rgb(HEX{i-1})/255;
    b = bar(i,-log10(pValues(i)),'FaceColor', C,'EdgeColor',C);
end
grid on
grid minor
line([0 14],[-log10(0.05) -log10(0.05)],'Color','black','LineWidth',2,'LineStyle','--');
hold off

%% Figure 3B
resultsDirectory = 'E:\Ahmed\HP\work\Results\NuclearReceptors\';
load([resultsDirectory 'xuCell2012_allNRs.mat'])
% select which experiment to display for each receptor
selectionInd = [1,2,1,1,2,2];
for i = 1 : numel(selectionInd)
    newP(i,:) = P{i}(selectionInd(i),:); expName{i} = N{i}{selectionInd(i)};
end

f = figure,
set(f,'Position',[200, 200, 1024, 512])
b = bar(-log10(newP));
b(1).FaceColor = [0.4 0.4 0.4]; b(1).EdgeColor = [0.4 0.4 0.4];
b(2).FaceColor = hex2rgb('E64438')/255; b(2).EdgeColor = hex2rgb('E64438')/255;
grid on
grid minor
legend({'Whole Brain', 'Hypothalamus'})
line([0 7],[-log10(0.05) -log10(0.05)],'Color','black','LineWidth',2,'LineStyle','--');

%% Figure 4A
filesDirectory = 'E:\Ahmed\HP\work\Data\NuclearReceptors\';
load([filesDirectory 'structures_expression_NoNormalization_All_coreg.mat']);
load('geneOfInterest.mat')
structures = {'Isocortex', 'OLF', 'HPF' 'CTXsp', 'STR', 'PAL', 'CB', ...
    'TH', 'HY', 'MB', 'P', 'MY'};
for i = 1 : length(geneOfInterest)
    idx = find(strcmpi(genes,geneOfInterest{i})==1);
    if length(idx) == 1
        dat(i,:) = avgExp(idx,:);
    else
        dat(i,:) = nanmean(avgExp(idx,:));
    end
end
dat = dat ./ repmat(dat(:,1),1,size(dat,2));

figure; imagesc(log2(dat(:,2:end)))
set(gca, 'YTick', 1:length(genes), 'YTickLabel', genes, 'FontSize', 8)
set(gca, 'XTick', 1:length(structures), 'XTickLabel', structures, 'FontSize', 8)
colormap('redbluecmap')
X = log2(dat(:,2:end));
hold on;
for i = 1:size(X,1)
   plot([.5,size(X,2)+.5],[i-.5,i-.5],'k-');
end
for i = 1:size(X,2)
    plot([i-.5,i-.5],[.5,size(X,1)+.5],'k-');
end
hold off
axis image %// equal scale on both axes
rotateXLabels(gca(), 45)

%% Figure 4B-D
[num txt] = xlsread('C:\Users\amahfouz\Dropbox\projects\Steroid receptors in the ABA\PNAS\Figures\Book1.xlsx',1);
corrVals = num;
genes = txt(3:end,1);
NRnames = txt(1,2:12:end);
structures = {'Isocortex', 'OLF', 'HPF' 'CTXsp', 'STR', 'PAL', 'CB', ...
    'TH', 'HY', 'MB', 'P', 'MY'};
transCorrVals = 2.^(-corrVals/100);
% HEX = {'70FF71', '7ED04B', '9AD2BD', '8ADA87', '98D6F9', '8599CC', ...
%     'F0F080', 'FF7080', 'E64438', 'FF64FF', 'FF9B88', 'FF9BCD'};
textC = lines(length(genes));

custom_map = csvread('reds_256_rgb.txt');
custom_map = custom_map / 255;

figure, hold on
S = 0;
subplot(1,3,1), imagesc(transCorrVals(:,1+(12*S):12+(12*S)))
set(gca, 'YTick', 1:length(genes), 'YTickLabel', genes, 'FontSize', 8)
caxis([0 0.5])
colormap(custom_map)
hold on;
X = transCorrVals(:,1+(12*S):12+(12*S));
for i = 1:size(X,1)
   plot([.5,size(X,2)+.5],[i-.5,i-.5],'k-');
end
for i = 1:size(X,2)
    plot([i-.5,i-.5],[.5,size(X,1)+.5],'k-');
end
hold off
axis image %// equal scale on both axes
S = 3;
subplot(1,3,2), imagesc(transCorrVals(:,1+(12*S):12+(12*S)))
% set(gca, 'YTick', 1:length(genes), 'YTickLabel', genes)
caxis([0 0.5])
colormap(custom_map)
hold on;
X = transCorrVals(:,1+(12*S):12+(12*S));
for i = 1:size(X,1)
   plot([.5,size(X,2)+.5],[i-.5,i-.5],'k-');
end
for i = 1:size(X,2)
    plot([i-.5,i-.5],[.5,size(X,1)+.5],'k-');
end
hold off
axis image %// equal scale on both axes
S = 4;
subplot(1,3,3), imagesc(transCorrVals(:,1+(12*S):12+(12*S)))
% set(gca, 'YTick', 1:length(genes), 'YTickLabel', genes)
caxis([0 0.5])
colormap(custom_map)
hold on;
X = transCorrVals(:,1+(12*S):12+(12*S));
for i = 1:size(X,1)
   plot([.5,size(X,2)+.5],[i-.5,i-.5],'k-');
end
for i = 1:size(X,2)
    plot([i-.5,i-.5],[.5,size(X,1)+.5],'k-');
end
hold off
axis image %// equal scale on both axes
hold off

%% Figure 4E
[num txt] = xlsread('C:\Users\amahfouz\Dropbox\projects\Steroid receptors in the ABA\PNAS\Figures\Figure 4_5\receptors_coregulators_VTA&SN_pValues.xlsx',1);
rankSumVals = num(:,1:2:end-1);
genes = txt(3:end,1);
NRnames = txt(1,2:2:end);
for N = 1 : length(NRnames)
    NRs{N} = NRnames{N}(1:strfind(NRnames{N},'_')-1);
end
transrankSumVals = 2.^(-rankSumVals/1000);
textC = lines(length(genes));

custom_map = csvread('reds_256_rgb.txt');
custom_map = custom_map / 255;

figure, 
imagesc(transrankSumVals)
set(gca, 'YTick', 1:length(genes), 'YTickLabel', genes)
set(gca, 'XTick', 1:length(NRs))
colormap(custom_map), colorbar
hold on;
X = transrankSumVals;
for i = 1:size(X,1)
   plot([.5,size(X,2)+.5],[i-.5,i-.5],'k-');
end
for i = 1:size(X,2)
    plot([i-.5,i-.5],[.5,size(X,1)+.5],'k-');
end
hold off
axis image %// equal scale on both axes

%% Supplementary Figure 3A
filesDirectory = 'E:\Ahmed\HP\work\Data\NuclearReceptors\';
geneOfInterest = {'Ar', 'Pgr', 'Esr1', 'Esr2', 'Nr3c1', 'Nr3c2'};
% load the whole brain data for normalization
% AR, Esr1 and Esr2: coronal
load([filesDirectory 'structures_expression_NoNormalization.mat'])
dat_brain([1,3,4],:) = avgExp([1,3,4],1);
% GR, MR, and Pgr: sagittal
load([filesDirectory 'structures_expression_NoNormalization_All.mat'])
dat_brain([6,5,2],:) = avgExp([4,6,8],1);
% load the VTA & SN data
% AR, Esr1 and Esr2: coronal
load([filesDirectory 'structures_expression_NoNormalization_C_dopa_receptors.mat'])
dat([1,3,4],:) = avgExp([1,3,4],:);
% GR, MR, and Pgr: sagittal
load([filesDirectory 'structures_expression_NoNormalization_All_dopa_receptors.mat'])
dat([6,5,2],:) = avgExp([4,6,8],:);
dat = dat ./ repmat(dat_brain,1,size(dat,2));
figure; imagesc(log2(dat)), colormap('redbluecmap')
X = dat;
hold on;
for i = 1:size(X,1)
   plot([.5,size(X,2)+.5],[i-.5,i-.5],'k-');
end
for i = 1:size(X,2)
    plot([i-.5,i-.5],[.5,size(X,1)+.5],'k-');
end
hold off
axis image %// equal scale on both axes

%% Supplementary Figure 3B
filesDirectory = 'E:\Ahmed\HP\work\Data\NuclearReceptors\';
load([filesDirectory 'structures_expression_NoNormalization_All_coreg.mat']);
load('geneOfInterest.mat')
for i = 1 : length(geneOfInterest)
    idx = find(strcmpi(genes,geneOfInterest{i})==1);
    if length(idx) == 1
        dat_brain(i,:) = avgExp(idx,1);
    else
        dat_brain(i,:) = nanmean(avgExp(idx,1));
    end
end
% load the VTA & SN data
filesDirectory = 'E:\Ahmed\HP\work\Data\NuclearReceptors\';
load([filesDirectory 'structures_expression_NoNormalization_All_dopa_coreg.mat']);
load('geneOfInterest.mat')
structures = {'VTA', 'SNr', 'SNc'};
for i = 1 : length(geneOfInterest)
    idx = find(strcmpi(genes,geneOfInterest{i})==1);
    if length(idx) == 1
        dat(i,:) = avgExp(idx,:);
    else
        dat(i,:) = nanmean(avgExp(idx,:));
    end
end
dat = dat ./ repmat(dat_brain,1,size(dat,2));

figure; imagesc(log2(dat))
set(gca, 'YTick', 1:length(genes), 'YTickLabel', genes, 'FontSize', 8)
set(gca, 'XTick', 1:length(structures), 'XTickLabel', structures, 'FontSize', 8)
colormap('redbluecmap')
X = log2(dat);
hold on;
for i = 1:size(X,1)
   plot([.5,size(X,2)+.5],[i-.5,i-.5],'k-');
end
for i = 1:size(X,2)
    plot([i-.5,i-.5],[.5,size(X,1)+.5],'k-');
end
hold off
axis image %// equal scale on both axes
rotateXLabels(gca(), 45)

%% Figure 1 - boxplot of correlaion with Esr1_coronal across 13 big structures
[num txt] = xlsread('C:\Users\amahfouz\Documents\MATLAB\Results\NuclearReceptors\xuCell2012_Esr1_All.xls',1);
corrVals = num(:,2:3:end);
structures = {'Whole Brain', 'Isocortex', 'OLF', 'HPF' 'CTXsp', 'STR', 'PAL', 'CB', ...
    'TH', 'HY', 'MB', 'P', 'MY'};
genes = txt(3:end,1);
ScatterBoxPlot(corrVals, 80, 'jet', structures)
set(gca, 'ylim', [-0.5 1])

%% Figure 3
[num txt] = xlsread('C:\Users\amahfouz\Dropbox\projects\Steroid receptors in the ABA\PNAS\Figures\Figure 4_5\receptors_coregulators_bigStructures.xlsx',1);
corrVals = num;
genes = txt(3:end,1);
NRnames = txt(1,2:12:end);
structures = {'Isocortex', 'OLF', 'HPF' 'CTXsp', 'STR', 'PAL', 'CB', ...
    'TH', 'HY', 'MB', 'P', 'MY'};
transCorrVals = 2.^(-corrVals/100);
% HEX = {'70FF71', '7ED04B', '9AD2BD', '8ADA87', '98D6F9', '8599CC', ...
%     'F0F080', 'FF7080', 'E64438', 'FF64FF', 'FF9B88', 'FF9BCD'};
textC = lines(length(genes));
% Select receptor
S = 0;
for y = 100 : 100 : floor(max(max(transCorrVals(:,1+(12*S):12+(12*S))))*10)*100
    X(1,y/100) = 2^(-(600-y)/100);
end
f = figure, hold on
for i = 1 : length(structures)
    scatter(ones(size(transCorrVals(:,i)))*i,transCorrVals(:,i+(12*S)), 50, textC, ...%hex2rgb(HEX{i})./255, ...
        'o', 'filled');
    colormap('lines')

    for j = 1 : size(transCorrVals(:,i),1)
        if transCorrVals(j,i+(12*S)) > 2^(-500/100)
            text(i+0.1,transCorrVals(j,i+(12*S)),genes(j), 'color', textC(j,:),...
                'FontSize', 15, 'FontWeight', 'bold');
        end
    end
end
hold off
line([0 numel(structures)+1], [2^(-500/100) 2^(-500/100)], ...
    'LineStyle', '--', 'LineWidth', 2, 'color', [0.5 0.5 0.5])
set(gca, 'xlim', [0 length(structures)+1], 'XTick', 1:numel(structures), ...
    'XTickLabel', structures, 'FontSize', 15)
rotateXLabels(gca(), 45)
set(gca, 'ylim', [0 (floor(max(max(transCorrVals(:,1+(12*S):12+(12*S))))*10)/10)+0.1], ...
    'YTick', X, 'YTickLabel', round(-100*log2(X)))
ylabel('Correlation Rank', 'FontSize', 20)
grid on
set(gcf,'units','normalized','outerposition',[0 0 1 1])
saveas(f,[NRnames{S+1} '.fig']);
saveas(f,[NRnames{S+1} '.png']);

%% Figure 4
[num txt] = xlsread('C:\Users\amahfouz\Dropbox\projects\Steroid receptors in the ABA\PNAS\Figures\Figure 4_5\receptors_coregulators_VTA&SN_pValues.xlsx',1);
rankSumVals = num(:,1:2:end-1);
genes = txt(3:end,1);
NRnames = txt(1,2:2:end);
for N = 1 : length(NRnames)
    NRs{N} = NRnames{N}(1:strfind(NRnames{N},'_')-1);
end
transrankSumVals = 2.^(-rankSumVals/1000);
textC = lines(length(genes));


figure, 
imagesc(transrankSumVals)
set(gca, 'YTick', 1:length(genes), 'YTickLabel', genes)
set(gca, 'XTick', 1:length(NRs))
colormap('redbluecmap'), colorbar

for y = 1000 : 1000 : floor(max(max(transrankSumVals(:)))*10)*1000
    X(1,y/1000) = 2^(-(6000-y)/1000);
end
f = figure, hold on
for i = 1 : length(NRs)
    scatter(ones(size(transrankSumVals(:,i)))*i,transrankSumVals(:,i), 50, textC, ...%hex2rgb(HEX{i})./255, ...
        'o', 'filled');
    colormap('lines')

    for j = 1 : size(transrankSumVals(:,i),1)
        if transrankSumVals(j,i) > 2^(-500/100)
            text(i+0.1,transrankSumVals(j,i),genes(j), 'color', textC(j,:),...
                'FontSize', 15, 'FontWeight', 'bold');
        end
    end
end
hold off
line([0 numel(NRs)+1], [2^(-5000/1000) 2^(-5000/1000)], ...
    'LineStyle', '--', 'LineWidth', 2, 'color', [0.5 0.5 0.5])
set(gca, 'xlim', [0 length(NRs)+1], 'XTick', 1:numel(NRs), ...
    'XTickLabel', NRs, 'FontSize', 15)
rotateXLabels(gca(), 45)
set(gca, 'ylim', [0 (floor(max(max(transrankSumVals(:)))*10)/10)+0.1], ...
    'YTick', X, 'YTickLabel', round(-1000*log2(X)))
ylabel('Sum of Correlation Ranks', 'FontSize', 20)
grid on
set(gcf,'units','normalized','outerposition',[0 0 1 1])

%% Figure 3
[num txt] = xlsread('C:\Users\amahfouz\Dropbox\projects\Steroid receptors in the ABA\PNAS\Figures\Book1.xlsx',1);
corrVals = num;
genes = txt(3:end,1);
NRnames = txt(1,2:12:end);
structures = {'Isocortex', 'OLF', 'HPF' 'CTXsp', 'STR', 'PAL', 'CB', ...
    'TH', 'HY', 'MB', 'P', 'MY'};
transCorrVals = 2.^(-corrVals/100);
% HEX = {'70FF71', '7ED04B', '9AD2BD', '8ADA87', '98D6F9', '8599CC', ...
%     'F0F080', 'FF7080', 'E64438', 'FF64FF', 'FF9B88', 'FF9BCD'};
textC = lines(length(genes));
% Select receptor
figure, hold on
S = 0;
subplot(1,3,1), imagesc(transCorrVals(:,1+(12*S):12+(12*S)))
set(gca, 'YTick', 1:length(genes), 'YTickLabel', genes)
caxis([0 0.5])
colormap('redbluecmap')
S = 3;
subplot(1,3,2), imagesc(transCorrVals(:,1+(12*S):12+(12*S)))
set(gca, 'YTick', 1:length(genes), 'YTickLabel', genes)
caxis([0 0.5])
colormap('redbluecmap')
S = 4;
subplot(1,3,3), imagesc(transCorrVals(:,1+(12*S):12+(12*S)))
set(gca, 'YTick', 1:length(genes), 'YTickLabel', genes)
caxis([0 0.5])
colormap('redbluecmap')
hold off

% S = 4;
% for y = 100 : 100 : floor(max(max(transCorrVals(:,1+(12*S):12+(12*S))))*10)*100
%     X(1,y/100) = 2^(-(600-y)/100);
% end
% f = figure, hold on
% for i = 1 : length(structures)
%     scatter(ones(size(transCorrVals(:,i)))*i,transCorrVals(:,i+(12*S)), 50, textC, ...%hex2rgb(HEX{i})./255, ...
%         'o', 'filled');
%     colormap('lines')
% 
%     for j = 1 : size(transCorrVals(:,i),1)
%         if transCorrVals(j,i+(12*S)) > 2^(-500/100)
%             text(i+0.1,transCorrVals(j,i+(12*S)),genes(j), 'color', textC(j,:),...
%                 'FontSize', 15, 'FontWeight', 'bold');
%         end
%     end
% end
% hold off
% line([0 numel(structures)+1], [2^(-500/100) 2^(-500/100)], ...
%     'LineStyle', '--', 'LineWidth', 2, 'color', [0.5 0.5 0.5])
% set(gca, 'xlim', [0 length(structures)+1], 'XTick', 1:numel(structures), ...
%     'XTickLabel', structures, 'FontSize', 15)
% rotateXLabels(gca(), 45)
% set(gca, 'ylim', [0 (floor(max(max(transCorrVals(:,1+(12*S):12+(12*S))))*10)/10)+0.1], ...
%     'YTick', X, 'YTickLabel', round(-100*log2(X)))
% ylabel('Correlation Rank', 'FontSize', 20)
% grid on
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
% saveas(f,[NRnames{S+1} '.fig']);
% saveas(f,[NRnames{S+1} '.png']);

