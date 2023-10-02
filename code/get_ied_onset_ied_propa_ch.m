% function [inside, q, t, n] = get_ied_onset_ied_propa_ch(subnum)
% the function requires,
% subnum = subject number (str) as input 
% and generates the following arrays with ieds registered during fmri scan,
% ch_list_onset, ch_list_propa, and ch_list_all
% Their format is given below, 
% 1st col. = subject number
% 2nd col. = event type
% 3rd col. = name of ied onset/propagation/all channel, depending on the
% array

function [ch_list_onset, ch_list_propa, ch_list_all] = get_ied_onset_ied_propa_ch(subnum)

% initialize output as empty array
ch_list_onset = {};   % ied onset channels
ch_list_propa = {};   % ied propagation channels
ch_list_all = {};   % all channels with ieds registered during fmri scan

ch_list = {};

% state path of .txt file listing ied onset and ied propagation channels
% for each event type of every subject
ied_ch_file_path = '/work/levan_lab/mtang/fmri_project/ied_onset_ied_propa_ch.txt';
% ied_ch_file_path = 'C:\Users\siumichael.tang\Downloads\fmri_project\ied_onset_ied_propa_ch.txt';
file_id = fopen(ied_ch_file_path, 'r');

% % read header of input file by appyling the format spec. (%s) four times, 
% % then scan the remainder of file using format spec. (%s %s %s %s),
% which denotes the four columns of data of string type,
% with empty space separating each field in format spec., but delimiter of
% comma is used for extracting data in each line.
data_header = textscan(file_id, '%s', 4, 'Delimiter', ',');
data_cell = textscan(file_id, '%s %s %s %s', 'Delimiter', ',');   % comma-delimited file

% close file using the file identifier
fclose(file_id);

% data_cell has the following format, 
% col. 1 = subject number
% col. 2 = event type
% col. 3 = names of ied onset channels, each separated by empty space
% col. 4 = names of ied propagation channels, each separated by empty space

% use strcmp (string compare) to find rows that have subject number equal
% to the one that is being examined
row_req = find(strcmp(data_cell{1}, subnum));

% check if row_req is empty (i.e. if subject number interested in included
% in the list of channels , if so, proceed to execute lines below, otherwise,
% output empty array and return control to line calling this function
if isempty(row_req)
    return
end

% scan across each line in data_cell, then scan through each location
% name recorded (col. 3 for ied onset channels and col. 4 for ied propagation
% channels), use regexp to search for alphabetic and numeric patterns, 
% e.g. given LA1-5, further decompose the number sequence to intergal
% increments, so that we have LA1-2, LA2-3, LA3-4, LA4-5

% initialize array with the following format, 
% col. 1 = subject number
% col. 2 = event type
% col. 3 = ied onset channel names, one name per cell
% col. 4 = ied propagation channel names, one name per cell
ch_list = {};
row_num = 1;   % initialize row number for ch_list
for i = 1:numel(row_req)   % for each row with matching sub. num. in data_cell
    ch_list{row_num, 1} = data_cell{1}{row_req(i)};   % assign sub. num.
    ch_list{row_num, 2} = data_cell{2}{row_req(i)};   % assign event type
    ch_list{row_num, 3} = {};   % initialize empty cell for ied onset ch. name
    ch_list{row_num, 4} = {};   % initialize empty cell for ied propa. ch. name
    
    % only proceed if ied onset ch. names in curr. row of data_cell is not
    % empty (i.e. name is provided)
    if ~isempty(data_cell{3}{row_req(i)})
        
        % separate channel names into cells
        ied_onset_locations = strsplit(data_cell{3}{row_req(i)}, ' ');
    
        % initialize ied onset ch. name array with its count
        pair_count = 1;
        ied_onset_ch_name = {};
    
        % for each channel location, use regexp to obtain it alphabetic 
        % and numeric components, then expand its numeric sequence, e.g.
        % if given LA1-5, expand it to LA1-2, LA2-3, LA3-4, etc.
        for item = 1:numel(ied_onset_locations)
            ch_type = regexp(ied_onset_locations{item}, '[A-Z]*', 'match');   % alphabetic compo.
            ele_con_num = regexp(ied_onset_locations{item}, '[\d]*', 'match');   % numeric comp.
  
            num_i = str2num(ele_con_num{1});   % intitial number of sequence
            num_f = str2num(ele_con_num{2});   % final number of sequence
            num_s = num_i;   % assign start number
    
            while num_s < num_f   % as long as start number is less than final number
                num_pair = num_s:num_s+1;   % form number pair
                ied_onset_ch_name{pair_count} = [ch_type{:}, num2str(num_pair(1)), '-', ...
                    ch_type{:}, num2str(num_pair(2))];   % construct channel name
                num_s = num_s + 1;   % increment start number by 1 for next entry
                pair_count = pair_count + 1;   % increment pair count by 1 for next enty
            end 
        end

        % after expanding all ied onset location names, take transpose and
        % assign results to col. 3 of ch_list
        ch_list{row_num, 3} = transpose(ied_onset_ch_name);
    
    end   % end if ~isempty(data_cell{3}{row_req(i)})
    
    % only proceed if ied propagation ch. names in curr. row of data_cell is not
    % empty (i.e. name is provided)
    if ~isempty(data_cell{4}{row_req(i)})

        % separate channel names in cells
        ied_propa_locations = strsplit(data_cell{4}{row_req(i)}, ' ');
    
        % initialize ied propagation ch. name array with its count
        pair_count = 1;
        ied_propa_ch_name = {}; 
        
        % for each channel location, use regexp to obtain it alphabetic 
        % and numeric components, then expand its numeric sequence, e.g.
        % if given LA1-5, expand it to LA1-2, LA2-3, LA3-4, etc.
        for item = 1:numel(ied_propa_locations)
            ch_type = regexp(ied_propa_locations{item}, '[A-Z]*', 'match');   % alphabetic compo.
            ele_con_num = regexp(ied_propa_locations{item}, '[\d]*', 'match');   % numeric comp.
    
            num_i = str2num(ele_con_num{1});   % intitial number of sequence
            num_f = str2num(ele_con_num{2});   % final number of sequence
            num_s = num_i;   % assign start number
    
            while num_s < num_f   % as long as start number is less than final number
                num_pair = num_s:num_s+1;   % form number pair
                ied_propa_ch_name{pair_count} = [ch_type{:}, num2str(num_pair(1)), '-', ...
                    ch_type{:}, num2str(num_pair(2))];   % construct channel name
                num_s = num_s + 1;   % increment start number by 1 for next entry
                pair_count = pair_count + 1;   % increment pair count by 1 for next enty
            end  
        end
        
        % after expanding all ied propagtion location names, take transpose and
        % assign results to col. 4 of ch_list
        ch_list{row_num, 4} = transpose(ied_propa_ch_name);
    
    end
    
    % increment row number by 1 for next row in data_cell with matching
    % subject number
    row_num = row_num + 1;
end   % end for i = 1:numel(row_req)   % for each row with matching sub. num. in data_cell

% convert cells into single layer array for each type of channels (onset,
% propagation, and all)
ch_list_onset = {};   % array for ied onset channels
ch_list_propa = {};   % array for ied propagation channels
ch_list_all = {};   % array for all channels with ieds registered during scan

row_num = 1;   % initialize row number as 1 
for i = 1:size(ch_list, 1)   % for each row in ch_list
    for j = 1:numel(ch_list{i, 3})   % for each item in cell of ied onset ch. name 
        ch_list_onset(row_num, 1:2) = ch_list(i, 1:2);
        ch_list_onset{row_num, 3} = ch_list{i, 3}{j};   % each item in 3rd col.
        row_num = row_num + 1;   % increment row number by 1 for next channel
    end
end

row_num = 1;   % initialize row number as 1 
for i = 1:size(ch_list, 1)   % for each row in ch_list
    for k = 1:numel(ch_list{i, 4})   % for each item in cell of ied onset ch. name
        ch_list_propa(row_num, 1:2) = ch_list(i, 1:2);
        ch_list_propa{row_num, 3} = ch_list{i, 4}{k};   % each item in 3rd col.
        row_num = row_num + 1;   % increment row number by 1 for next channel
    end
end

% concatenate ied onset and propa arrays vertically to obtain array for all
% channels with ieds registered during scan
ch_list_all = [ch_list_onset; ch_list_propa];

end   % end function ch_list = get_ied_onset_ied_propa_ch(subnum)