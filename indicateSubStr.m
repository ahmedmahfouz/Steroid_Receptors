%%% 15 Feb. 2013
%%% indicate substructures of a structure (defined by the acronym)
%%% the first ID in the output list is the id of the input structure
%%% NOTE:  THE FUNCTION RETURNS STRUCTURE ID (THE REFERENCE VOLUME INCLUDES THESE IDs)

function children = indicateSubStr(inputStr, filesDir);

% load the ontology data
ontologyFile = [filesDir 'ARAontology.xls'];
[num txt] = xlsread(ontologyFile);
strAcronyms = txt(2:end,5);
clear num; clear txt;

% find the index of the input structure
inputIndex = find(strcmpi(strAcronyms, inputStr));

if isempty(inputIndex)
    display('EROOR: structure not found')
else
    % load the list of structures
    load([filesDir 'listOfStructures.mat']);
    listOfSubStr = listOfStructures(inputIndex,1);
    if length(listOfStructures(inputIndex,1).children) ~= 0
        listOfSubStr = extractChildrenInfo(listOfStructures(inputIndex,1), listOfSubStr);
        for i = 1 : size(listOfSubStr,1);
            if ~isempty(listOfSubStr(i,1).id)
                children(i) =  listOfSubStr(i,1).id;
            end
        end
    else
        children = listOfSubStr.id;
        disp([inputStr ' has no children'])
    end   
end







