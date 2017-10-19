%%% 15 Feb. 2013
%%% calculate the correlation between two genes within a
%%% specific structure (defined by the ABA acronym)
%%%to run on the server
%%% 22 April 2013: 'TYPE' is an option to run the correlation based on the
%%% coronal experiments only (C), sagittal only (S), or on both (MIXED)

function topCorrGenes(geneOfInterest, strOfInterest, TYPE)

% define the files directory
if ispc
    filesDirectory = 'C:/Users/amahfouz/Documents/MATLAB/Data/NuclearReceptors/files/';
    coronalDir = 'C:/Users/amahfouz/Documents/MATLAB/Results/NuclearReceptors/geneCorr/CoronalOnly/';
    sagittalDir = 'C:/Users/amahfouz/Documents/MATLAB/Results/NuclearReceptors/geneCorr/SagittalOnly/';
    allDir = 'C:/Users/amahfouz/Documents/MATLAB/Results/NuclearReceptors/geneCorr/All/';
else
    filesDirectory = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Data/NuclearReceptors/files/';
    coronalDir = '/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/NuclearReceptors/geneCorr/CoronalOnly/';
    sagittalDir ='/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/NuclearReceptors/geneCorr/SagittalOnly/';
    allDir ='/tudelft.net/staff-bulk/ewi/insy/DBL/amahfouz/MATLAB/Results/NuclearReceptors/geneCorr/All/';
end

if strcmp(TYPE, 'C')
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
    
    % check wether this gene has a coronal experiment or not
    geneExpCheck = find(strcmpi(listOfGenes, geneOfInterest) == 1);
    if isempty(geneExpCheck)
        display('seed gene does not have a coronal experiment');
    else
        % take the natural logarithm of the expression energy values
        % NOTE: voxels with no expression of a gene have a value of "-1", we adda a
        % value of "+2" to all data points so the -1 -> 0 after taking the log
        load([filesDirectory 'expressionMatrix_CoronalOnly.mat']);
        % expMat_logTransformed = log(fullExpressionMatrix + 2);
        expMat_logTransformed = expressionMatrix_CoronalOnly;
        clear expressionMatrix_CoronalOnly;   

        % define the outputfile
        if ~exist(coronalDir, 'dir')
            mkdir(coronalDir);
        end
        geneDirectory = [coronalDir geneOfInterest '/'];
        if ~exist(geneDirectory, 'dir')
            mkdir(geneDirectory);
        end
        resultsDirectory = [geneDirectory geneOfInterest '_' strOfInterest '_CoronalOnly/'];
        if ~exist(resultsDirectory, 'dir')
            mkdir(resultsDirectory);
        end

        % calculate the correlation and return the correlation vector and the
        % list of voxels used for the calculation
        corrWithGene(geneOfInterest, strOfInterest, listOfGenes, listOfExperimentNumbers, allExpPlanes, ...
            voxelAnnotation_CoronalOnly, brainInd_CoronalOnly, filesDirectory, expMat_logTransformed, resultsDirectory);
    end
    
elseif strcmp(TYPE, 'S')
    % load gene names
    load([filesDirectory 'allGenes.mat']);
    % load experiments' numbers
    load([filesDirectory 'allExpNumbers.mat']);
    % load experiments' plane
    load([filesDirectory 'allExpPlanes.mat']);
    % load voxels annotations
    load([filesDirectory 'voxelAnnotation.mat']);
    % load brain voxels 
    load([filesDirectory 'brainInd.mat']);
    
    %%% select sagittal experiments only
    sagInd = 4346 : length(allExpNumbers);
    
    % check wether this gene has a sagittal experiment or not
    geneExpCheck = find(strcmpi(allGenes(sagInd), gene) == 1);
    if isempty(geneExpCheck)
        display('seed gene does not have a sagittal experiment');
    else
        % take the natural logarithm of the expression energy values
        % NOTE: voxels with no expression of a gene have a value of "-1", we adda a
        % value of "+2" to all data points so the -1 -> 0 after taking the log
        load([filesDirectory 'fullExpressionMatrix.mat']);
        % expMat_logTransformed = log(fullExpressionMatrix + 2);
        expMat_logTransformed = fullExpressionMatrix;
        clear fullExpressionMatrix;
        
        % define the outputfile
        if ~exist(sagittalDir, 'dir')
            mkdir(sagittalDir);
        end
        geneDirectory = [sagittalDir geneOfInterest '/'];
        if ~exist(geneDirectory, 'dir')
            mkdir(geneDirectory);
        end
        resultsDirectory = [geneDirectory geneOfInterest '_' strOfInterest '_SaigittalOnly/'];
        if ~exist(resultsDirectory, 'dir')
            mkdir(resultsDirectory);
        end

        % calculate the correlation and return the correlation vector and the
        % list of voxels used for the calculation
        corrWithGene(geneOfInterest, strOfInterest, allGenes(sagInd), allExpNumbers(sagInd), allExpPlanes, ...
            voxelAnnotation, brainInd, filesDirectory, expMat_logTransformed(sagInd,:), resultsDirectory);

    end
    
elseif strcmp(TYPE, 'All')
    % take the natural logarithm of the expression energy values
    % NOTE: voxels with no expression of a gene have a value of "-1", we adda a
    % value of "+2" to all data points so the -1 -> 0 after taking the log
    load([filesDirectory 'fullExpressionMatrix.mat']);
    % expMat_logTransformed = log(fullExpressionMatrix + 2);
    expMat_logTransformed = fullExpressionMatrix;
    clear fullExpressionMatrix;

    % load gene names
    load([filesDirectory 'allGenes.mat']);
    % load experiments' numbers
    load([filesDirectory 'allExpNumbers.mat']);
    % load experiments' plane
    load([filesDirectory 'allExpPlanes.mat']);
    % load voxels annotations
    load([filesDirectory 'voxelAnnotation.mat']);
    % load brain voxels 
    load([filesDirectory 'brainInd.mat']);

    % define the outputfile
    if ~exist(allDir, 'dir')
        mkdir(allDir);
    end
    geneDirectory = [allDir geneOfInterest '/'];
    if ~exist(geneDirectory, 'dir')
        mkdir(geneDirectory);
    end
    resultsDirectory = [geneDirectory geneOfInterest '_' strOfInterest '/'];
    if ~exist(resultsDirectory, 'dir')
        mkdir(resultsDirectory);
    end

    % calculate the correlation and return the correlation vector and the
    % list of voxels used for the calculation
    corrWithGene(geneOfInterest, strOfInterest, allGenes, allExpNumbers, allExpPlanes, ...
        voxelAnnotation, brainInd, filesDirectory, expMat_logTransformed, resultsDirectory);
end







