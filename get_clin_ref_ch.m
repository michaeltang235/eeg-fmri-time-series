% This function generates an array of channels based on input file which
% contains the names of a range of electrodes determined clinically.
% e.g. if a record, 'RAI1-5', is shown in the input file, this function
% generates a series of channels corresponding to this record, namely
% RAI1-2, RAI2-3, RAI3-4, and RAI4-5. 
% Then, recall that the input signal array, which is constructed in
% rho_spikes.m, lists the ids, names, and signals of all channels
% available.
% This function will then match each entry in the series of channels formed
% (A) with entries in the signal array (B), locating if a channel in (A) also 
% appears in (B). If so, that particular channel is saved in the ouput 
% variable, otherwise, it is discarded.

% CH_NAME_CLINICAL = GET_CLIN_REF_CH(SUBNUM, CLIN_REF_CH_PATH, SIG_BOX_ALL)
% SUBNUM is the subject number interested (int),
% CLIN_REF_CH_FILE_PATH is the full path of file containing info. about
% ref. channels recorded clinically by doctors
% SIG_BOX_ALL is the array which contains ids, names, and signals of all
% channels recorded in the subject
% CH_NAME_CLINICAL is an output array with 1st layer listing the types 
% of ref. channels determined clinically, and 2nd layer listing
% the ids and names of those reference channels

function ch_name_clinical = get_clin_ref_ch(subnum, clin_ref_ch_file_path, sig_box_all)

% initialize output var. ch_name_clinical, whose values will change
% depending on if string of channel names is empty in input file
ch_name_clinical = {};

% open input file, assign read permission, and obtain file identifier 
file_id = fopen(clin_ref_ch_file_path, 'r');

% read header of input file by appyling the format spec. (%s) two times, 
% then scan the remainder of file using format spec. (%d %s), which denotes
% the first column as integer and second column as string, with delimiter
% of empty space separating each field in format spec.
data_header = textscan(file_id, '%s', 2, 'Delimiter', ' ');
data_cell = textscan(file_id, '%d %s', 'Delimiter', ' ');

% close file using the file identifier
fclose(file_id);

% find row index corresponding to subject number interested
% subnum = 14;
row_number = find(data_cell{1} == subnum);

% initialize exec_num = 0 (not execute lines below, 1 = execute lines)
exec_num = 0;

% check if cell is empty, if so, print message to console and output empty
% array, otherwise, separate strings from cell by ','
if isempty(data_cell{2}{row_number})
    disp('warning, no channel names recorded in clinical file')
else
    ele_array = strsplit(data_cell{2}{row_number}, ',');
    exec_num = 1;   % set num. to 1 so lines below will be executed
end

% if exec_num == 1 (cell containing channel names is not empty in input
% file), then execute the following lines
if exec_num == 1
    
% initialize array for storing channel names 
    ch_name = cell(numel(ele_array), 1);

% for each electrode string detected from input, 
% extract alphabetical and numeric components, 
% then use numbers obtained to construct sequence, e.g. 'LOR1-5' is
% separated into 'LOR' and '1-5' components, then forming a sequence 1-2,
% 2-3, 3-4, 4-5. Afterwards, channel names are LOR1-2, LOR2-3, etc.
    for j = 1:numel(ele_array)
        ch_type = regexp(ele_array{j}, '[A-Z]*', 'match');   % alphabetic compo.
        ele_con_num = regexp(ele_array{j}, '[\d]*', 'match');   % numeric comp.
        num_i = str2num(ele_con_num{1});   % intitial number of sequence
        num_f = str2num(ele_con_num{2});   % final number of sequence
        num_s = num_i;   % assign start number
        pair_count = 1;   % assign count of number pairs
        while num_s < num_f   % as long as start number is less than final number
            num_pair = num_s:num_s+1;   % form number pair
            ch_name{j}{pair_count} = [ch_type{:}, num2str(num_pair(1)), '-', ...
            ch_type{:}, num2str(num_pair(2))];   % construct channel name
            num_s = num_s + 1;   % increment start number by 1 for next entry
            pair_count = pair_count + 1;   % increment pair count by 1 for next entry
        end
    end   % end for j = 1:numel(ele_array)
       
% with channel names from input file established, scan through each entry
% in array and match it with entries in sig_box_all, which contains signals
% of all channels recorded in functional images. 
% this step is done because sometimes channel names entered in clinical
% report are missing in electrode file

% initialize array for storing matched channel names, 1st layer denotes
% electrode type and 2nd layer denotes ch. names of current type, with ids
% in col. 1 and names in col. 2
    ch_name_clinical = cell(numel(ch_name), 1);  
    row_ind = 1;   % initialize row index as 1 for current electrode type
    for entry = 1:numel(ch_name)   % for each type of electrode est. in ch_name
        for i = 1:length(ch_name{entry})   % for each channel name est. in current type
            for j = 1:size(sig_box_all, 1)   % for each row in sig_box_all
                if strcmp(ch_name{entry}{i}, sig_box_all{j, 2})   % if names in both arryas matched 
                    ch_name_clinical{entry}{row_ind, 1} = sig_box_all{j, 1};   % id of channel
                    ch_name_clinical{entry}{row_ind, 2} = ch_name{entry}{i};   % name of channel
                    row_ind = row_ind + 1;   % increment row index by 1 for next match
                    break   
                end
            end
        end
        row_ind = 1;   % increment row index by 1 for next electrode type
    end   % end if entry = 1:numel(ch_name)
    
end   % end if exec_num == 1 (if cell in containing channel names is not empty)
    

end   % end function ch_name_clinical = get_clin_ref_ch(ref_ch_file_path, sig_box_all)
