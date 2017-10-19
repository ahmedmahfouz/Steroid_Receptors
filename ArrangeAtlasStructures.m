%%% 14 Feb. 2013
%%% Arrange Atlas structures

clear all;
addpath('jsonlab');

dataDir = 'E:\Ahmed\HP\work\Data\NuclearReceptors\files';
ontData = loadjson([ dataDir '/Ontology.json']);
listOfStructures(1,1) = ontData.msg;

if length(ontData.msg.children) ~= 0;
    listOfStructures = extractChildrenInfo(ontData.msg, listOfStructures);
end

for i = 1 : size(listOfStructures,1);
    
    if isempty(listOfStructures(i,1).atlas_id)
        structureAtlasID(i) = 0;
    else
        structureAtlasID(i) =  listOfStructures(i,1).atlas_id;
    end
    
    if i == 1
        structureParentID(i) =  0;
    else
        structureParentID(i) =  listOfStructures(i,1).parent_structure_id;
    end
    
    structureID(i) =  listOfStructures(i,1).id;
    structureAcronym{i} = listOfStructures(i,1).acronym;
    structureColor{i} = listOfStructures(i,1).color_hex_triplet;
    structureName{i} = listOfStructures(i,1).name;
end

outFile = [dataDir '/ARAontology.xls'];
xlswrite(outFile, [{'Atlas ID'} {'Structure ID'} {'Parent ID'} {'Acronym'} {'Color_HEX'} {'Name'}], 1, 'B1');
xlswrite(outFile, structureAtlasID', 1, 'B2');
xlswrite(outFile, structureID', 1, 'C2');
xlswrite(outFile, structureParentID', 1, 'D2');
xlswrite(outFile, structureAcronym', 1, 'E2');
xlswrite(outFile, structureColor', 1, 'F2');
xlswrite(outFile, structureName', 1, 'G2');





