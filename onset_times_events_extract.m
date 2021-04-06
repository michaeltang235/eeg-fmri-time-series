close all
clear all

% This script reads event onset times .txt file(s) and outputs struct. 
% named terms that stores the times in seconds for all runs interested
% in a .mat file.

% ** use this script only if event times in each run are already reported
% in separate .txt file, otherwise, use onset_times_events.m **
% The only diff. in output struct. generated between the two scripts is
% that a field named "durations" are included in the output struct.
% produced by this script (onset_times_event_extract.m), while that is 
% excluded when using the other script (onset_times_event_extract.m),
% becasue event durations reported in input files applicable to this script
% are sometimes non-zero.

% Remarks:
% (i) filename of event times file must follow the format of 'Run%d_EV%d.txt',
% where the first '%d' denotes run number and the second '%d' denotes event
% types interested
% (ii) each input file represent one event type. If a run has 2 event files, 
% then there are 2 event types in that run.
% (iii) event onset times must be stored in first column of input file
% (iv) times reported in input files are relative to fmri trigger
% (v) event durations are assumed to be zero (close to instantaneous),
% unless otherwise specified in col. 3 of input files

% Updates: April 6, 2021, allowed non-zero event durations to be included in
% output struct. If col. 3 of input files have value of 1, event durations
% are assumed to be zero, otherwise, values reported in that col. are
% extracted and saved as the event durations required.
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

% get unique run numbers identified in input files
uni_runnum = unique(runnum);

% create numevents array to store number of events reported in each run
numevents = cell(1, length(uni_runnum));   % initialize cell array 

% scan through each entry in uni_runnum array (B)
% use ismember to return logical indices if entries in runnum array (A)
% are in B. If entries in A are in B, value of 1 for those entries are
% returned.
% then, use nnz to count number of non-zero entries returned to get number
% of events in each run
for entry = 1:length(uni_runnum)   % for each unique run number
    numevents{entry} = nnz(ismember(runnum, uni_runnum{entry}));
end

% similarly, event types represented by input files are shown in the
% filenames, e.g. in 'Run1_EV2.txt', onset times of event type 2 are stored
% get indices required in filename, \d denotes numeric digit, \w* denotes
% matching any group (more than or equal to 1 char) comprised of
% alphabetic, numeric, or underscore character, all the way to the char before
% '.txt', the file extension, for cases like _EV2_5.txt,
%
sind_event = regexp(file_name_event, '_EV\d', 'end');   % start index of event type
eind_event = regexp(file_name_event, 'EV\w*', 'end');   % end index of event type
% eind_event = regexp(file_name_event, 'EV\d+', 'end');   % end index of event type

% assign event defintion to array eventdef (numeric) and eventdefstr (str)
eventdef = cell(1, numel(uni_runnum));   % initialize eventdef array (numeric)
eventdefstr = cell(1, numel(uni_runnum));   % initialize eventdefstr array (str)

rowind = 1;   % initialize rowind as 1
for i = 1:numel(uni_runnum)   % for each unqiue event
    for j = 1:numevents{i}   % for each event type identified in each unique run
        % use rowind to access the corresponding entry in file_name_event,
        % extract event type in filename
        strext = file_name_event{rowind}(sind_event{rowind}:eind_event{rowind}); 
        % replace '_' with ' ' (if there is such pattern) 
        % and convert from str to numeric
        strextrep = strrep(strext, '_', ' ');   % str replaced
        eventdef{i}{j} = str2num(strextrep);   % convert from str to numeric
        
        % get place holder by number of digits found in current event def
        placestr = repmat('%d_', 1, length(eventdef{i}{j}));   
        % format event def. into string, get rowind^th run number,
        eventdefstr{i}{j} = ['type_', sprintf(placestr, eventdef{i}{j}),...
            'run_', runnum{rowind}];
        rowind = rowind + 1;   % increment rowind by 1
    end
end

% check if entries in eventdef agree with events entered
for i = 1:numel(eventdef)
    for j = 1:numel(eventdef{i})
        if ~(ismember(eventdef{i}{j}, cell2mat(events)))
            sprintf('event types identified do not match with events entered')
        end
    end
end

% scan through each entry in uni_runnum array (B)
% use ismember to return logical indices if entries in runnum array (A)
% are in B. If entries in A are in B, value of 1 for those entries are
% returned.
% then, use nnz to count number of non-zero entries returned to get number
% of events in each run
for entry = 1:length(uni_runnum)   % for each unique run number
    numevents{entry} = nnz(ismember(runnum, uni_runnum{entry}));
end

% create onset times array (ontimesar) and duration array (duratar) to store
% event onset times (col. 1) and event durations (col. 2) reported, repsectively,
% in each input file
ontimesar = cell(1, numel(uni_runnum));   % initialize cell array ontimesar
duratar = cell(1, numel(uni_runnum));   % initialize cell array duratar

% initialize rownumber as 1 for accessing entry in file_path cell array
rownum = 1;  

for i = 1:numel(uni_runnum)   % for each unique run number
    for j = 1:numevents{i}   % for each event identified in current run
        % select the rownum^th path from file_path array,
        % and read the entire file 
        inputfile_ar = readcell(file_path{rownum});  
        % assign 1st column of file to ontimesar, convert from cell to double
        ontimesar{i}{j} = cell2mat(inputfile_ar(:,1));
        
        % if entries in col. 2 aren't equal to 1, extract corresponding
        % values reported and save them as double,
        % otherwise, zero event duration is assumed
        if inputfile_ar{1,2} ~= 1   % if entries in col. 2 aren't equal to 1
            duratar{i}{j} = cell2mat(inputfile_ar(:,2));   % extract values
        else   % if col. 2 contains entries of value of only 1
            duratar{i}{j} = zeros(numel(inputfile_ar(:,2)), 1);   % assume zero durations
        end
       
        rownum = rownum + 1;   % increment rownum by 1
    end
end

% create struct. named terms to output event onset times and related info
terms = struct;

% assemble terms struct using the variables calcu. above
for i = 1:numel(uni_runnum)   % for each unique run number
    fieldname = sprintf('run%s', uni_runnum{i});   % create fieldname, i^th uni_runnum
    terms.(fieldname).results.num_cond = numevents{i};   % number of conditions estaliblished in each run
    terms.(fieldname).results.event_def_str = eventdefstr{i};   % string describing def. of each event for spec. run
    terms.(fieldname).results.event_def = eventdef{i};   % defition of each event for spec. run
    terms.(fieldname).results.onsettimes = ontimesar{i};   % onset times recorded in dataset provided
    terms.(fieldname).results.evti = ontimesar{i};   % initial time (seconds) of each event in each run, the same as ontimesar
    terms.(fieldname).results.durations = duratar{i};  % event durations (seconds)
end

% to keep struct. format consistent with terms struct. created using
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



    






