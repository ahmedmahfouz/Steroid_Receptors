%%% 11 March 2013
%%% give the number of children structures and the number of voxels spanned
%%% by each of the input structures

filesDirectory = 'files/';
load([filesDirectory 'voxelAnnotation.mat']);

[num txt] = xlsread([filesDirectory 'structuresOfInterest.xlsx']);
strList = txt(22:end,3);
clear num; clear txt;

for i = 1 : length(strList)
    strChildren = indicateSubStr(strList(i), filesDirectory);
    numOfChildren(i) = length(strChildren);
    inclVoxelsInd = find(ismember(voxelAnnotation, strChildren) == 1);
    numOfVoxels(i) = length(inclVoxelsInd);
    clear strChildren; clear inclVoxelsInd;
end

outFile = [filesDirectory 'HY_children_voxels.xls'];
xlswrite(outFile, {'Structure'}, 1, 'B1');
xlswrite(outFile, strList, 1, 'B2');
xlswrite(outFile, {'Number of children'}, 1, 'C1');
xlswrite(outFile, numOfChildren', 1, 'C2');
xlswrite(outFile, {'Number of voxels'}, 1, 'D1');
xlswrite(outFile, numOfVoxels', 1, 'D2');
