globals().clear()

import re
#import pandas as pd
from pandas import read_csv
import numpy as np
from prelim_setup import get_path

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# BEGIN USER INPUT

# enter subject number (str)
subnum = '14'

# enter directory and filename of input file(s)
directname = '/work/levan_lab/mtang/fmri_project/' + 'sub' + subnum
# directname = ['/Users/michaeltang/Downloads/fmri_project/', 'sub', subnum, '_imthres0_exmask'];

# enter formats of input file(s)
evfilename_format = ['Sub' + subnum + '.*\.ev2', '^Run\d+_.*\.txt']

# enter output directory and filename 
directname_output = directname + '/py_matrices/'
filename_output = 'onset_times_events_sub' + subnum + '.npy'

# enter eeg sampling frequency in Hz
sampfreq = 2000

# enter time in seconds of first fmri frame relative to beginning of eeg 
# recording (if known, otherwise leave it empty ([]) and it will be determined 
# in the script)
if subnum == '35':
    fmri_frame_ti = 34.3647
else:
    fmri_frame_ti = [];


# END USER INPUT
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# format full path of event types file
ev_type_path = '/work/levan_lab/mtang/fmri_project/event_types.txt'

# load data from event types file into dataframe
ev_type_df = read_csv(ev_type_path)

# find row with matching subject number in ev_type_df
row_matched = ev_type_df.iloc[:, 0] == int(subnum)

# obtain subset of dataframe using row number found, convert subset to list
ev_type_list = ev_type_df[row_matched].iloc[:, -1].tolist()

# remove whie space and split string by ' ' 
ev_type_str = ev_type_list[0].strip().split(' ')

# event type with more than one def. contains the pattern '-', which is used 
# separate each def. check if an entry contains '-', if so, further split string
# into components by '-', append data to ev_type_read list
ev_type_read = []   # intiailize list
for entry in ev_type_str:   # for every entry in list
    if '-' in entry:   # if '-' exists in current entry
        ev_type_read.append(entry.split('-'))   # split string by '-', append to list
    else:   # if '-' doesn't exist in curr. entry, append entry to list
        ev_type_read.append(entry)
        
# assign ev_type_read to events
events = ev_type_read
        
# search for matching filenames in directory specified, use get_path to otain
# full path of event files and its format number
for i in range(len(evfilename_format)):   # every format input
    if get_path(directname, evfilename_format[i]):   # if matching filename is found
        ev_file_path = get_path(directname, evfilename_format[i])   # obtain full path
        ev_file_format_num = i   # assign ev file format number
        break
    
# format of event files (.ev2 or .txt) determines the way event onsets are obtained:

#----------------------------------------------------------
# obtain event onsets by the following if input files are in .ev2 format (ev_file_format_num == 0) 
if ev_file_format_num == 0:
    
    # initialize output dict for curren run
    op_dict = {}
    
    for run_ind in range(len(ev_file_path)):
        
        # get indices of current run number embedded in path
        s_ind = re.search('Run\d.*_', ev_file_path[run_ind], re.IGNORECASE).span()
        
        # get current run number
        runnum = ev_file_path[run_ind][s_ind[0]+3:s_ind[1]-1]
    
        # load data in .ev2 file into numpy array
        event_array = np.loadtxt(ev_file_path[run_ind])


    # check format of onset times reported in dataset 
    # i.e. are times reported in (seconds) or (integer offset)
    # format_bool = TRUE means onset times in dataset are in (offset)
    # format_bool = FALSE means onset times in dataset are in (seconds)
    # note: output onset times are all in unit of seconds relative to the 
    # start of the beginning of fmri acquisition
        format_bool = event_array[0, -1]/sampfreq > 1

    # determine fmri_ti, the time at which first fmri frame was aqcuired
    # relative to the start of eeg recordings by locating the first occurrence of 
    # type 0 (2nd col.) and response 5 (3rd col.)
        if not bool(fmri_frame_ti):   # if no fmri trigger is entered
            for i in range(event_array.shape[0]):   # every row of event array read
               # identify row of first occurence of event type 0 and response 5
               if event_array[i, 1] == 0 and event_array[i, 2] == 5:
                   # obtain the corresponding time from rightmost col.
                   fmri_ti = event_array[i, -1]   
                   break
        else:   # if fmri trigger is provided by user, assign value to fmri_ti
            fmri_ti = fmri_frame_ti
    
# use for loop to go through each event, for each value in the def. 

# (e.g. def. of event{i} = [1, 11], find matching rows in event array,
# extract the corresponding onset time from rightmost col.
# check unit of times reported, and convert the values to unit of
# seconds if necessary, assign resultant value to evti nested list
# e.g.evti[i^th] contains onset times in seconds for event[i^th]

# initialize ind_total list with format given below,
# in i^th entry contains the matching row indices in event array for the i^th 
# event type. If no matching row indices are found, empty list is assigned for
# that event type
        ind_total = []    
        for i in range(len(events)):   # each def. of event type
        #for i in range(1):
            ev = events[i]   # obtain current event type
            ind_found = []
            # event type could contain more than one def., e.g. [1, 11]
            # find row index in event array with matching event type (2nd col.) for each def.
            for j in range(len(ev)):   # every def. of current event type
            # find matching row index for current event def.
                cur_ind_matched = np.where(event_array[:, 1] == int(ev[j]))[0]   # numpy array
                ind_found.extend(cur_ind_matched.tolist())
                #ind_found.append(cur_ind_matched)
                #ind_found = ind_found + list(cur_ind_matched)   # concatenate lists
            # after looping through all def. of current event type,
            # append all row indices found for current event type to ind_total nested list
            ind_total.append(ind_found)   

# with row indices for every event type obtained, get the corresponding onset
# times in seconds
# initialize evti list format given below
# i^th entry contains onset times in seconds for the i^th event type in events list
# emtpy list is assigned when no onset time is found for the event type
        evti = []
        for item in ind_total:   # every entry in ind_total nested list
            # for cases where there are more than one def. for an event type, 
            # convert indices from list to set, then get the union of the sets to obtain
            # unique row indices, and then sort those indices by ascending order
            ind_req = sorted((set(item)))       
            # check unit of times reported, and calculate onset times with respect to the 
            # start of fmri acquisiton, conert numpy array to list, and append list to 
            # evti nested list
            if ~format_bool:
                evti.append((event_array[ind_req, -1] - fmri_ti).tolist())
            else:
                evti.append(((event_array[ind_req, -1] - fmri_ti)/sampfreq).tolist())
        
            # get frequency of occurence for every event type
            evfreq = [len(item) for item in evti]
            
        # create key for current run
        key_name_cur = 'run' + runnum
        
    # assign event def., event onsets in seconds (evti), and their freq. (evfreq) 
    # to output dict. under current key 'results'
        op_dict[key_name_cur] = {}
        op_dict[key_name_cur]['results'] = {}
        op_dict[key_name_cur]['results']['event_def'] = events
        op_dict[key_name_cur]['results']['evti'] = evti
        op_dict[key_name_cur]['results']['evfreq'] = evfreq
        
#----------------------------------------------------------  
        
# obtain event onsets by the following if input files are in .txt format (ev_file_format_num == 1)     
    
if ev_file_format_num == 1:
    
    # initialize lists with format given below,
    # i^th entry in evti denotes onset times for the run number in i^th entry
    # of runnum, and event def in the i^th entry of ev_def
    # e.g. evti[0] gives onset times for the event type defined in ev_def[i]
    # for the run number defined in runnum[i]
    runnum = []
    ev_def = []
    evti = []
    
    # ger current filename
    for i in range(len(ev_file_path)):
        file_name_cur = re.search('Run\d+_.*\.txt', ev_file_path[i]).group() 
    
        # get run number from cur. file name by first finding the start and end 
        # indices of the search pattern in in cur. file name
        eind = re.search('Run\d+', file_name_cur).span()

        # eind is a tuple with two entries, the first and last index,
        # use the last index to get current run numer
        runnum.append(file_name_cur[3:eind[-1]])
    
        # similarly, get event type contained in cur. file name
        eind_event = re.search('EV\d+', file_name_cur).span()
    
        # use indices obtained above to get current event type
        ev_def.append(file_name_cur[eind_event[0]+2:eind_event[-1]])
        
        # load data from .txt file into numpy arrray
        event_array = np.loadtxt(ev_file_path[i])
        
        # get event onsets from the leftmost column in numpy array loaded
        evti.append(event_array[:, 0].tolist())
        
    # after obtaining all necessary info. from all runs, convert runnum lis to
    # set to get unique entries, and sort the set and convert it backt to list
    unique_runnum = list(sorted(set(runnum)))   # get unique run number
    
    # initialize output dict
    op_dict = {}
    
    # for every unique run number, select corresponding indices in runnum
    # then select the associated event def. and event onsets (evti),
    for item in range(len(unique_runnum)):
        ind_req = [i for i, j in enumerate(runnum) if j == unique_runnum[item]]
        ev_def_cur_run = np.array(ev_def)[ind_req].tolist() 
        evti_cur_run = np.array(evti, dtype = object)[ind_req].tolist()
        evfreq_cur_run = [len(entry) for entry in evti_cur_run]
        
        # after obtaining data for current run, save them in output dict
        fd_name = 'run' + unique_runnum[item]   # name key by current run number
        op_dict[fd_name] = {}   # initialize key
        op_dict[fd_name]['results'] = {}
        op_dict[fd_name]['results']['event_def'] = ev_def_cur_run
        op_dict[fd_name]['results']['evti'] = evti_cur_run
        op_dict[fd_name]['results']['evfreq'] = evfreq_cur_run
        
    
#----------------------------------------------------------  


#sub_list = ['14', '15', '18', '19', '20', '27', '30', '31', '32', '33', '34', \
           #'35', '41']
           
output_path = directname_output + filename_output
    
np.save(output_path, op_dict)            
    
        
        
    
    


