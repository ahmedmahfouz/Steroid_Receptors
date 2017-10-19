%%% 14 Feb. 2013
%%% a function to extract children's info

function list = extractChildrenInfo(root, list)

numberOFChildren = length(root.children);
if numberOFChildren ~= 0;
    for i = 1 : numberOFChildren
        listSize = size(list,1);
        list(listSize+1,1) = root.children(1,i);
        
        if length(root.children(1,i).children) ~= 0;
            list = extractChildrenInfo(root.children(1,i), list);
        end
    end
end

