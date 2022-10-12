ied_ch_file_path = '/work/levan_lab/mtang/fmri_project/ied_onset_ied_propa_ch.txt';
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

subnum = '41';

% use strcmp (string compare) to find rows that have subject number equal
% to the one that is being examined, then use find to get non-zero entries
% from output of strcmp, assign results to row_req
row_req = find(strcmp(data_cell{1}, subnum));

% check if row_req is empty (i.e. if subject number interested in included
% in the list of channels , if so, proceed to the lines below, otherwise,
% return control to line calling this function
if isempty(row_req)
    return
end

ch_list = {};
row_num = 1;
for i = 1:numel(row_req)
    ch_list{row_num, 1} = data_cell{1}{row_req(i)};
    ch_list{row_num, 2} = data_cell{2}{row_req(i)};
    ch_list{row_num, 3} = {};
    ch_list{row_num, 4} = {};
    
    if ~isempty(data_cell{3}{row_req(i)})
        
    ied_onset_locations = strsplit(data_cell{3}{row_req(i)}, ' ');
    
    pair_count = 1;
    ied_onset_ch_name = {};
    
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
    
    ch_list{row_num, 3} = transpose(ied_onset_ch_name);
    
    end   % end if ~isempty(data_cell{3}{row_req(i)})
    
    if ~isempty(data_cell{4}{row_req(i)})
        ied_propa_locations = strsplit(data_cell{4}{row_req(i)}, ' ');
    
        pair_count = 1;
        ied_propa_ch_name = {}; 
        
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
        
        ch_list{row_num, 4} = transpose(ied_propa_ch_name);
    
    end
    
    row_num = row_num + 1;
end