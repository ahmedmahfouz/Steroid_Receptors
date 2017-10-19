%%% 1 March 2013
%%% Run a batch of experiments to retrieve top correlated genes with
%%% different NRs in multiple brain regions

if ispc
    filesDirectory = 'E:/Ahmed/HP/work/Data/NuclearReceptors/files/';
    dataDir = 'C:\Users\amahfouz\SURFdrive\Projects\Lisa_Koorneef\Data\';
else
    filesDirectory = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Data/NuclearReceptors/files/';
    dataDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/Lisa_Koorneef/';
end

% genes = {'Ar', 'Nr3c1', 'Esr1', 'Esr2', 'Nr3c2', 'Pgr'};
% genes = {'Neurod1', 'Neurod2'};
% genes = {'Crh', 'Crhr1','Crhr2'};
genes = {'Avp', 'Oxt'};
% structures = {'HIP','CA','DG','CA1','CA2','CA3','DG-mo','DG-po','DG-sg','STR','TH','Isocortex'};
% structures = {'DG-mo','DG-po','DG-sg'};
% structures = {'HIP', 'CA', 'DG', 'CA1', 'CA2', 'CA3'};
% structures = {'PALd', 'PALv', 'PALm', 'PALc' 'BST', 'BAC'};
% structures = [structures {'BST', 'CEA', 'sAMY', 'NTS', 'VTA', 'SNr', 'SNc', 'DR', 'ACB'}];
% structures = {'BST', 'CEA', 'CTXsp', 'HIP', 'PAL', 'sAMY'};
%%% HY substructures
% [num txt] = xlsread([filesDirectory 'HY_substructures.xls']);
% hySubStr = txt(2:end,3);
% structures = {'CA1slm' 'CA1so', 'CA1sp', 'CA1sr', 'CA3slm', 'CA3slu', 'CA3so', 'CA3sp', 'CA3sr'};
% structures = {'VTA', 'SNr', 'SNc', 'DR'};
% structures = {'HIP','CA','DG','CA1','CA2','CA3','DG-mo','DG-po','DG-sg','STR','TH','Isocortex',...
%     'CA1slm' 'CA1so', 'CA1sp', 'CA1sr', 'CA3slu', 'CA3so', 'CA3sp', 'CA3sr'}


% mainStr = {'grey', 'Isocortex', 'OLF', 'HPF' 'CTXsp', 'STR', 'PAL', 'CB', ...
%     'TH', 'HY', 'MB', 'P', 'MY'};
% % read the regions of interest
% T = readtable([dataDir 'structures.csv']);
% structure.acronym = T.structure_acronym;
% clear T
% structures = [mainStr structure.acronym'];
% structures(31) = [];

structures = {'Hy'};

experimentType = 'All';
extension = 'csv';

if ~ispc
    filesDirectory = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Data/NuclearReceptors/files/';
    %%%(1) calculate the correlations
    for i = 1 : length(genes)
        for j = 1 : length(structures)
            disp([genes{i} '_' structures{j}])
            topCorrGenes(genes{i}, structures{j}, experimentType);
            disp('..... DONE')
        end
    end
% else
%     filesDirectory = 'E:/Ahmed/HP/work/Data/NuclearReceptors/files/';
%     resDirectory = 'E:\Ahmed\HP\work\Results\NuclearReceptors\';
%     geneCorr = 'E:\Ahmed\HP\work\Results\NuclearReceptors\geneCorr/';
%     xlsFiles = 'E:\Ahmed\HP\work\Results\NuclearReceptors\xlsFiles/';
    resDirectory = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/NuclearReceptors/';
    geneCorr = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/NuclearReceptors/geneCorr/';
    xlsFiles = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/NuclearReceptors/xlsFiles/';
    %%%(2) save results to xls files
%     idx = [2,2,3]; %crh, crhr1, crhr2
%     idx = [1,1, 4]; %Slc6a3, Slc6a4, Dbh
    for i = 1 : length(genes)
        for j = 1 : length(structures)
            disp([genes{i} '_' structures{j}])
%             xlsCorrRes(genes{i}, structures{j}, experimentType);
%             K = idx(i);
            tic
            txtCorrRes(genes{i}, structures{j}, experimentType, extension,...
                filesDirectory, resDirectory, geneCorr, xlsFiles, K);
            disp('..... DONE')
            toc
        end
    end
end





