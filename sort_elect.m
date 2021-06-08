% With array containing coordinates of electrodes obtained, this function
% sorts electrodes according to their types and outputs cell array with 
% each cell storing data of those of the same type. e.g. given a set of 
% electrodes with names {LPI1, LPI2, ROF1, ROF2, LOF1, LOF2}, the types of
% electrodes identified are LPI, ROF, and LOF. Hence, the output array will
% contain one cell for each type identified, with all info. such as
% coordinates and name tags stored within itself.

% ELECT_SORTED = SORT_ELECT(ELECT_AR)
% ELECT_SORTED is the output cell array with electrodes sorted
% ELECT_AR is the input cell array with all electrodes known

% REMARKS:
% input electrode array has to follow the format outlined below, 
% row 1 contains headers of columns
% column 2 contains name tags of eletrodes, 
% column 3:5 contain their x, y, and z-coordinates, respectively

%---------------------------------------------------------------------------

function [elect_sorted] = sort_elect(elect_ar)

%---------------------------------------------------------------------------
% PART (I): clean electrode input array

% file containing info. about electrodes may contains rows with value of 0
% assigned to missing name tags and coordinates,
% remove such rows and return resultant array to electar_cleaned

% initialize empty array named electar_cleaned
electar_cleaned = {};

% use loop to go through each row of electar array , if name tag (col. 2) 
% is NOT missing, assign info. of that row to electar_cleaned array
rownum = 1;   % initialize row number as 1
for i = 1:size(elect_ar, 1)   % for each row in electar
    % if name tag is found or if coordinates are not missing
    % ismissing finds missing entries in array, have to separtely evaluate
    % col. 2 and col. 3:5, as they are of different data types, 
    % any means if any of the entries is non-zero or logical 1
    if any(~ismissing(elect_ar{i, 2})) && any(~ismissing([elect_ar{i, 3:5}]))
        electar_cleaned(rownum, :) = elect_ar(i,:);   % assign relevant info. to array
        rownum = rownum + 1;   % increment row number by 1
    end
end

% END PART (I): clean electrode input array
%---------------------------------------------------------------------------
% PART (II): sort cleaned electrode array

% initialize cell array named ele_sorted with first entry having info. about
% the first electrode in the cleaned electrode array (electar_cleaned)
ele_sorted{1} = electar_cleaned(2, 2:5);

% initialize electrode count
ele_count = 1;

% each electrode name contains a compo. of string and a compo. of number
% e.g. 'RPF1' has string compo. of 'RPF' and numeric compo. of '1'
% use for loop to go through each entry in the name column, 
% find max. index of string in name tag and number associated with name,
% compare name in current name tag with the next one, 
% if names match and number of the next one is larger than the current
% one, store the info. within current electrode type, 
% else create another cell with new name tag

for i = 2:size(electar_cleaned, 1)-1
    ele_name = electar_cleaned{i, 2};   % name of current electrode 
    ele_name1 = electar_cleaned{i+1, 2};   % name of next electrode 
    
    if length(ele_name) <= length(ele_name1)
        maxind = max(regexp(ele_name, '[A-Z]'));   % find max. index of string in name
        if strcmp(ele_name(1:maxind), ele_name1(1:maxind)) && ...
                str2double(ele_name(maxind+1:end)) < str2double(ele_name1(maxind+1:end))   % if names match and num. agree
            ele_sorted{ele_count} = [ele_sorted{ele_count}; electar_cleaned(i+1,2:5)];   % add info. of next one to current type
        else
            ele_count = ele_count + 1;   % increment electrode count by 1
            ele_sorted{ele_count} = electar_cleaned(i+1,2:5);   % add info. to new cell
        end
    else
        ele_count = ele_count + 1;   % increment electrode count by 1
        ele_sorted{ele_count} = electar_cleaned(i+1,2:5);   % add info. to new cell
    end
end

% print msg. after sorting
sprintf('electrodes sorted, %d of them found', ele_count)

% assign sorted electrode array to output variable, elect_sorted
elect_sorted = ele_sorted;

% END PART (II): sort cleaned electrode array
%---------------------------------------------------------------------------
 


