% With format of filename known, this function first locates files
% interested in directory specified and output the full paths of those files

% [NUMFILE, FILE_PATH] = GET_PATH(DIRECTNAME, FILENAME_FORMAT)
% NUMFILE is number of files interested located
% FILE_PATH is full path of files interested
% FILE_NAME is names of files interested
% DIRECTNAME is name of directory where files interested are located at
% FILENAME_FORMAT is format of filename interested
%---------------------------------------------------------------------------

function [numfile, file_path, file_name] = get_path(directname, filename_format)

% with format of filename entered, use dir to search and 
% get full path of file(s) in directory specified, store results 
% in file_struct
file_struct = dir(fullfile(directname, filename_format));  

% if no file is found, an empty struct is formed,
% display warning message if no file is found
if isempty(file_struct)   % check if struct. is empty
    disp('warning! when using get_path, no file is found.')
end

% each file found is given a dimension in struct, get number of
% files located and store their corresponding path to cell
% array named file_path
% similarly, store filename found in struct named file_name
numfile = length(file_struct);   % number of file(s) found
file_path = cell(numfile, 1);   % initialize array for storing paths of files interested
file_name = cell(numfile, 1);   % initialize array for storing names of files interested
for i = 1:numfile
    file_path{i, 1} = fullfile(file_struct(i).folder, file_struct(i).name);
    file_name{i, 1} = file_struct(i).name;
end


end



