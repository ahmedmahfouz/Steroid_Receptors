%%% 18 Feb. 2013
%%% A function that returns a sorted list of genes based on the correlation
%%% with a gene of interest based on the expression within a specific
%%% structure

function corrWithGene(gene, str, geneslist, expNums, expPlanes, ...
    vAnn, bIndx, fDir, eMat, resDir)

% find the index of the gene within the list of genes
gene_index = find(strcmpi(geneslist, gene) == 1);
if isempty(gene_index)
    display('seed gene not found');
else
    gene_experimentNos = expNums(gene_index);
    gene_experimentPlanes = expPlanes(gene_index);

    % get the list of voxels covering the structure of interest
    inclVoxAnnotations = indicateSubStr(str, fDir);
    inclVoxelsInd = find(ismember(vAnn, inclVoxAnnotations) == 1);

    % brain indecies corresponding to str
    strInds = bIndx(inclVoxelsInd);

    % prepare data for the correlation function
    geneOfInterestExp = eMat(inclVoxelsInd, gene_index);
    dataMat = eMat(inclVoxelsInd, :);

    % calculate the correlation between all probe-pairs of the gene of interest 
    % and all the other genes
    % NOTE: calculations are done pair-wise, because at each step we want to
    % ignore voxels with missing values
    for i = 1 : size(geneOfInterestExp,2)
        [notMissing1 ~] = find(geneOfInterestExp(:,i) ~= -1);
        for j = 1 : size(dataMat,2)
            [notMissing2 ~] = find(dataMat(:,j) ~= -1);
            nm = intersect(notMissing1, notMissing2);
            if length(nm) >= 5
                [c(i,j) pval(i,j)] = corr(geneOfInterestExp(nm,i), dataMat(nm,j), 'type', 'Pearson');
                nmVoxels(i,j) = length(nm);
            else
                c(i,j) = -2;
                pval(i,j) = -1;
                nmVoxels(i,j) = 0;
            end
            clear notMissing2; clear nm;
        end
        clear missing1;
    end

    % save the results
    save([resDir gene '_' str '.mat'], 'c');
    save([resDir gene '_' str '_PVAL.mat'], 'pval');
    save([resDir str '_voxels.mat'], 'inclVoxelsInd');
    save([resDir str '_strInds.mat'], 'strInds');
    save([resDir str '_nmVoxels.mat'], 'nmVoxels');
    save([resDir gene '_expNumbers.mat'], 'gene_experimentNos');
    save([resDir gene '_expPlane.mat'], 'gene_experimentPlanes');
end




