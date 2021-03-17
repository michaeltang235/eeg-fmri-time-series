close all 
clear all
%---------------------------------------------------------------------------
% This script takes file(s) the following inputs:
% (i) path to input file(s)
% (ii) format the filenames (.ev2) containing details of events recorded 
% (iii) def. of events (neurologist markings, e.g. 21 and 41, ...),
% (iv) eeg sampling frequency 
% (v) fmri trigger in seconds (if known, optional)
% and outputs struct. array named terms which contains 
% the following:
% (i): number of conditions established in each run
% (ii): string of definition of each event in each run
% e.g. type_1_22_run_3
% (iii): definition of each event (i.e. types each event contains) (array)
% (iv): onsettimes array that (col. 1) lists the times provided in input 
% dataset, (col. 2) onset times in unit of seconds with respec to the first 
% fmri image acquistion
% (iv): onset times of each event in unit of (s) from col. 2 of onsettimes array
% (v): frequency of occurence of each event in each run

% Remarks: 

% *** assume format of input file(s) as follows: 
% 2nd column denotes event types, 
% 3rd column denotes response, and
% 6th column denotes event times in (seconds), but sometimes they are
% expressed as integer offsets instead, which correspond to the number of 
% EEG samples from the beginning of the recording.
% ***

% - directory specified only stores the necessary number of ev.2 file(s)
% - string "runX_ev" must be contained in input filename

% - diff. between onset_times_events_1.m and onset_times_event.m are
% that 
% (i): when onset time (ti) of an event is not found in a certain run, instead of
% setting ti of that regressor as 0 s, this new version simply doesn't even
% create a regressor for that event type in that run. Thus, total number of
% conditions (regressors) in all runs depends on the number of event types
% sucessfully identified in input file(s).
% (ii): less manual labor work is required for running, in that only
% entering number of runs and format of input filename, and the rest will be
% handled.
% (iii): when no event is identified in a run, an empty cell is assigned to 
% the onset times array, a new field "results" is constructued to store
% non-empty entries of all relevant arrays, while keeping the orginal
% (untrimmed) data in separate field "data".

% Updates: 
% (i) March 4, 2021, instead of using eventmat(1,6) (1st row, 6th column)
% as the time when the first func. image was aquired, find the first
% occurrence of such event (type 0 and response 5), and only consider
% event types that happened after such event was found.

% (ii) March 16, 2021, parenthese in event_def_str is removed to avoid
% confusion. Also, allow the user to enter fmri trigger when *.ev2 file(s) 
% contains no such information (lacking type 0 and response 5).

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter subject number (str)
subnum = '27';

% enter directory and filename of input file(s)
directname = ['/Users/michaeltang/Downloads/fmri_project/', 'sub', subnum, '_imthres0_exmask'];

% % enter format of input file(s)
evfilename_format = ['Sub', subnum, '*.ev2'];   % %wildcard char is added for searching matched files in directname

% enter event types interested
% e.g. {[1, 11], [12], [13]} denotes event 1 is composed of types 1 and 11, 
% event 2 is composed of type 12, while event 3 is composed of type 13
events = {[1 11], [12], [13]}; 

% enter output directory and filename 
directnameoutput = [directname, '/matrices'];
filenameoutput = ['onset_times_events_sub', subnum, '.mat'];

% enter eeg sampling frequency
sampfreq = 2000;   % Hz

% enter time in seconds of first fmri frame relative to beginning of eeg 
% recording (if known, otherwise leave it empty ([]) and it will be determined
% in the script)
fmri_frame_ti = [];

% enter if output should be written to path (1 = yes, 0 = no)
op_results = 0;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

%----------------

% with format of *.ev2 filename entered, use dir to search and 
% get full path of file(s)
evfile_path = dir(fullfile(directname , evfilename_format));   

% each .ev2 file found is given a dim. in struct. evfile_path to store 
% its info., get number of .ev2 files located in directory 
input_num = size(evfile_path, 1);

% print message to terminal listing number of .ev2 files found for this
% subject
sprintf('subject %s, number of .ev2 files found = %d', subnum, input_num)

% access name of each .ev2 file from struct. evfile_path
evfilename = cell(input_num, 1);
for i = 1:input_num
    evfilename{i} = evfile_path(i).name;
end

% initialize array for outputting onset times and relevant info.
terms = {};
for i = 1:input_num
    fn_run = sprintf('run%d', i);
    terms.(fn_run) = onset_times(directname, evfilename{i}, sampfreq, events, fmri_frame_ti);
end

% output struct. array to directory
if op_results == 1
    save(fullfile(directnameoutput,filenameoutput),'terms');
end

%--------------------------------

% define function named onset_times that requires the following inputs:
% (i) directory of event details file, 
% (ii) filename of event details (*.ev2 files), 
% (iii) sampling freq., 
% (iv) array(s) of event types interested, and 
% (v) time in seconds when the first fmri frame was acquired (optional)
% then output the following: 
% (i) number of conditions established in each run
% (ii) strings of event definitions
% (iii) onset times matrix (in terms of integer offset and time)
% (iv) events' intitial times (unit of seconds)
% (v) freq. of each event
% remarks: input arguments of function must be placed in the order defined

function run = onset_times(directname, evfilename, sampfreq, events, fmri_trigger)

    % import data from files to struct
    eventdata = importdata(fullfile(directname, evfilename));   % event data for this run
    
    % get event matrix of each run from struct.
    eventmat = eventdata.data;  
    
    % get number of events provided 
    numevents = length(events);
    
    % create eventdef array, use for loop to assign the i^th variable input
    % as the i^th event (e.g. define types 21 and 41 as event 1, types 12
    % as event 2, and so on ...)
    eventdef = cell(1, numevents);
    for i = 1:numevents
        eventdef{i} = events{i};
    end

    % create empty cell array for storing onset times (both recorded time (col. 1) 
    % and time (s) relative to the first fmri image (col. 2))
    ontimesar = {};
    
    % initialize array evfre for storing freq. of each event 
    evfreq = cell(1, numevents);
    for i = 1:numevents
        evfreq{i} = 0;
    end
    
    % check format of onset times reported in dataset 
    % i.e. are times reported in (seconds) or (integer offset)
    % formatnum = 0 means onset times in dataset are in (seconds)
    % formatnum = 1 means onset times in dataset are in (offset)
    % note: output onset times are all in unit of seconds relative to the 
    % start of the beginning of fmri acquisition
    formatnum = logical(eventmat(1,6)/sampfreq > 1);
   
    % determine fmri_ti, the time at which first fmri frame was aqcuired
    % relative to the start of eeg recordings, by checking if the 5th input
    % arg. is an empty matrix. If so, locate the first occurrence of type 0 
    % and response 5, otherwise, assign the 5^th input arg. to the variable.
    if isempty(fmri_trigger)   % if 5th input arg. is empty
        % get time and row number in eventmat at which first func. image was 
        % acquired (fmri_ti), which corresponds to neurological marking 
        % of event type 0 and respone 5
        for i = 1:size(eventmat, 1)   % each row
            if eventmat(i, 2) == 0 && eventmat(i, 3) == 5   % col. 2 = type and col. 3 = response
                fmri_ti = eventmat(i, 6);   % get offset from col. 6 
                fmri_ti_rownum = i;   % get correpsonding row number
                break   % break loop as only first occurrence is interested
            end
        end
        % print fmri_ti and the associated row number in current input file to terminal
        sprintf('In %s, row %d, first fmri frame was aquired at %.3f (offset or seconds)', ...
            evfilename, fmri_ti_rownum, fmri_ti)
    end
    
    % if the 5th input arg. is not empty, assign the 5^th input arg. 
    % to fmri_ti and set assoicated row number as 1, so all rows including
    % the 1st one in input file (.ev2) are scanned in the loop defined below
    if ~(isempty(fmri_trigger)) % if 5th input arg. is not empty
        fmri_ti = fmri_trigger;   % assign the 5^th input arg. to var.
        fmri_ti_rownum = 1;   % set row number
        sprintf('In %s, first fmri frame was set as the 5th input arg. (%.3f)', evfilename, fmri_trigger)
    end
   
    % use for loop to go through each event, for each value in the def. 
    % of the i^th event, check if it is equal to any row in col. 2 of eventmat.
    % Only consider rows after fmri_ti (happened after fmri acquisition
    % had begun)
    % (e.g. def. of event{i} = [1, 11], check which row in col. 2 of
    % eventmat has the same value (1 or 11))
    % if so, extract the corresponding onset time from column 6, 
    % check unit of times reported, and convert the values to unit of
    % seconds if necessary, assign resultant value to ontimesar
    % e.g. ontimesar{3}(5,2) contains the 5^th onset time in (s) of event 3
    
    for i = 1:numevents   % for each event
        for row = fmri_ti_rownum:size(eventmat, 1)   % for each row in eventmat after fmri aquisition began
            for j = 1:length(eventdef{i})   % for each value in def. of event
                if eventmat(row, 2) == eventdef{i}(j)   % if entry in eventmat matches with value in def.
                    evfreq{i} = evfreq{i} + 1;   % increment ev. freq. by 1
                    ontimesar{i}(evfreq{i}, 1) = eventmat(row, 6);   % get onset time from col. 6 in eventmat
                    if formatnum == 0   % check unit of onset times
                        ontimesar{i}(evfreq{i}, 2) = (ontimesar{i}(evfreq{i}, 1) - fmri_ti);   % calcu. onset times when unit is in seconds
                    end
                    if formatnum == 1
                        ontimesar{i}(evfreq{i}, 2) = (ontimesar{i}(evfreq{i}, 1) - fmri_ti)/sampfreq;   % calcu. onset times when unit is in offset
                    end
                end
            end
        end
        if evfreq{i} == 0   % if no rows in eventmat matched with current def. of events
            ontimesar{i} = [];   % assign empty matrix
        end
    end
    
    % extract initial times (ti) in unit of seconds of each event type from ontimesar
    evti = cell(1, numevents);
    for i = 1:numevents
        if isempty(ontimesar{i})   % if current cell is empty (no marking found)
            evti{i} = [];
        else   % if current cell is not empty
            evti{i} = ontimesar{i}(:,2);   % from col. 2 in each cell 
        end
    end
    
    % create string describing event details in each run 
    sind = regexp(lower(evfilename), 'run');   % start index of 'r' in event filename
    eind = strfind(lower(evfilename), '_ev');   % index of '_ev' in event filename
    runnum = evfilename(sind+3:eind-1);   % string of run number contained in event filename
    
    % assemble required strings
    eventdefstr = cell(1, numevents);
    for i = 1:numevents
        % create place holder by num. of types in each event
        placestr = repmat('%d_', 1, length(eventdef{i})); 
        % assign formatted string to eventdefstr array
        eventdefstr{i} = sprintf(append('type_', placestr, 'run_%s'), eventdef{i}, runnum);
    end
 
    % empty matrix is assigned to ontimesar if no event type is found in
    % eventmat, only select non-empty cells  
    
    % create matrix named nonemptyind for storing indices of non-zero
    % entries in evfreq
    nonemptyind = [];
    for j = 1:length(evfreq)
        if evfreq{j} ~= 0
            nonemptyind = [nonemptyind j];
        end
    end
    
    % number of conditions established in each run is equal to the number 
    % of sets of non-empty event types identified, get length of nonemptyid 
    % to assign value to var. num_cond
    num_cond = length(nonemptyind);
     
    % empty matrix is assigned to ontimesar if no event type is found in eventmat
    % only select non-empty cells (i.e. event with non-zero freq. recorded)
    % assign NON-EMPTY values to field "results" under struct. "run"   
    run.results.num_cond = num_cond;   % number of conditions estaliblished in each run
    run.results.event_def_str = eventdefstr(nonemptyind);   % string describing def. of each event for spec. run
    run.results.event_def = eventdef(nonemptyind);    % defition of each event for spec. run
    run.results.onsettimes = ontimesar(nonemptyind);   % onset times recorded in dataset provided
    run.results.evti = evti(nonemptyind);   % initial time (seconds) of each event in each run
    
     % output quantities calcu. (including EMPTY values) to field "data" 
     % under struct. "run"
    run.data.event_def_str = eventdefstr;   % string describing def. of each event for spec. run
    run.data.event_def = eventdef;   % defition of each event for spec. run
    run.data.onsettimes = ontimesar;   % onset times recorded in dataset provided
    run.data.evti = evti;   % initial time (seconds) of each event in (s)
    run.data.evfreq = evfreq;   % freq. of each event 
        
    
end   % end function defined
%--------------------------------

%---------------------------------------------------------------------------


