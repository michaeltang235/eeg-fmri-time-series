close all
clear all

% This script reads event onset times .txt file(s) and outputs struct. that
% stores the times in seconds for all runs interested in a .mat file.

% Remarks:
% (i) use this script only if event times in each run are already reported
% in separate .txt file, otherwise, use onset_times_events.m
% (ii) filename of event times file must follow the format of 'Run%d_EV%d.txt',
% where the first '%d' denotes run number and the second '%d' denotes event
% types interested
% (iii) each input file represent one event type. If a run has 2 event files, 
% then there are 2 event types in that run.
% (iv) event onset times must be stored in first column of input file

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter subject number (str)
subnum = '14';

% enter path to directory where all input files are located
directname = ['C:\Users\siumichael.tang\Downloads\fmri_project\', 'sub', subnum];
% directname = ['/Users/michaeltang/Downloads/fmri_project/', 'sub', subnum, '_imthres0_exmask'];

% enter format of filename of input file(s)
filename_format = 'Run*EV*.txt';   % wildcard char is added for searching

% enter event types interested
% e.g. {[1, 11], [12], [13]} denotes event 1 is composed of types 1 and 11, 
% event 2 is composed of type 12, while event 3 is composed of type 13
events = {[1], [2]}; 

% enter output directory and filename 
directname_output = [directname, filesep, 'matrices'];
filename_output = ['onset_times_events_extract_sub', subnum, '.mat'];

% enter if output should be written to path (1 = yes, 0 = no)
op_results = 0;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% get full path of event times file(s) using get_path
[numfile, file_path, file_name_event] = get_path(directname, filename_format);

% get run number from input filenames by first finding (i) the end index of
% 'Run\d' in filename, then (ii) end index of last numeric digit, 
% \d+ denotes matching consecutive numeric digits as many as possible
sind = regexp(file_name_event, 'Run\d', 'end');   % start index of run number in filename
eind = regexp(file_name_event, 'Run\d+', 'end');   % end index of run number in filename

% after indices required are stored in cell array, loop through each entry
% and extract run number, store values to runnum array
% i^th entry in runnum array denotes the run number of the i^th entry in
% file_path
runnum = cell(numel(file_name_event), 1);   % initialize cell array
for i = 1:numel(file_name_event)   % for each entry in file_name_event
    % extract run number using corresponding indices obtained
    runnum{i} = file_name_event{i}(sind{i}:eind{i});
end

% similarly, event types represented by input files are shown in the
% filenames, e.g. in 'Run1_EV2.txt', onset times of event type 2 are stored
% get indices required in filename, \d denotes numeric digit, \d+ denotes
% matching consecutive numeric digits as many as possible, for cases like _EV25.txt
sind_event = regexp(file_name_event, '_EV\d', 'end');   % start index of event type
eind_event = regexp(file_name_event, 'EV\d+', 'end');   % end index of event type

% assign event defintion to array eventdef (numeric) and eventdefstr (str)
eventdef = cell(1, numel(runnum));   % initialize eventdef array (numeric)
eventdefstr = cell(1, numel(runnum));   % initialize eventdefstr array (str)
for i = 1:numel(runnum)   % for each input file
    % extract event type in filename, convert from str to numeric
    eventdef{i} = str2num(file_name_event{i}(sind_event{i}:eind_event{i}));  
    % assemble event def. string, get i^th eventdef, and corresponding
    % run number 
    eventdefstr{i} = ['type_', num2str(eventdef{i}), '_run_', runnum{i}];
end

% check if eventdef agrees with events entered
for i = 1:numel(eventdef)
    if ~(ismember(eventdef{i}, cell2mat(events)))
        sprintf('event types identified do not match with events entered')
    end
end

% create onset times array (ontimesar) to store event onset times reported
% in each input file, whose 1st column stores onset times interested
ontimesar = cell(1, numel(runnum));   % initialize cell array ontimesar
% i^th entry in ontimesar represents onset times for the 1^th input file
for i = 1:numel(runnum)   % for each input file
    inputfile_ar = readcell(file_path{i});   % read entire i^th file
    % assign 1st column of file to ontimesar, convert from cell to double
    ontimesar{i} = cell2mat(inputfile_ar(:,1));
end

% create struct. named terms to output event onset times and related info
terms = struct;

% initialize entry index (entind) as 1 for assigning cell to their
% corresponding field (run)
entind = 1;   % initialize entry index as 1, 

% the i^th entries in eventdef, eventdefstr, and ontimesar cell arrays
% agree with the i^th input file provided, whose run number is also
% reported in the i^th entry in runnum array

% assemble terms struct using the variables calcu. above
for i = 1:numel(runnum)   % for i^th entry in runnum array, aka the i^th input file
    fieldname = sprintf('run%s', runnum{i});   % create fieldname, i^th runnum
    if ~(isfield(terms, fieldname))   % if field doesn't exist
        % create new field with current run number
        % assign the i^th cell in eventdefstr to current run number,
        % do the same for other variables 
        terms.(fieldname).results.num_cond = 1;   % initialize number of conditions estaliblished in each run
        terms.(fieldname).results.event_def_str = eventdefstr(i);   % string describing def. of each event for spec. run
        terms.(fieldname).results.event_def = eventdef(i);   % defition of each event for spec. run
        terms.(fieldname).results.onsettimes = ontimesar{i};   % onset times recorded in dataset provided
        terms.(fieldname).results.evti = ontimesar{i};   % initial time (seconds) of each event in each run, the same as ontimesar
        entind = i;   % update current entry index 
    else   % if fieldname already exists
        % assign entind:i^th entries in eventdefstr array to var
        terms.(fieldname).results.event_def_str = eventdefstr(entind:i);
        terms.(fieldname).results.event_def = eventdef(entind:i);
        terms.(fieldname).results.onsettimes = ontimesar(entind:i);
        terms.(fieldname).results.evti = ontimesar(entind:i);    
        terms.(fieldname).results.num_cond = numel(ontimesar(entind:i));   
    end
end

% to struct. format consistent with terms struct. created using
% onset_times_events.m, another field 'data' is created to store all var.
% calculated. When using onset_times_events.m, onset times of certain event type
% might be missing in input files (.ev2), the corresponding entry in onset
% times array is then given a value of 0. Under the field 'data', it
% contains all var. calculated without omitting any 0 entries found.
% In this script, onset times are extracted from input files directly, but
% to keep format of struct. consistent, assign field 'data' with var. in
% field 'results'
for i = 1:length(fieldnames(terms))   % for each field (run number)
    run_name = fieldnames(terms);   % get fieldname
    terms.(run_name{i}).data = terms.(run_name{i}).results;   % assign var. in 'results' to 'data'
end

%---------------------------------------------------------------------------
% output struct. array to directory
if op_results == 1
    save(fullfile(directname_output,filename_output),'terms');
end



    






