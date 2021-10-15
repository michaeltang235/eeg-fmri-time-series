
import os 
import re

import pandas as pd

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# GET_PATH FUNCTION:

# With format of filename known, this function first locates files
# interested in directory specified and output the full paths of those files

# GET_PATH(DIRECTNAME, FILENAME_FORMAT)
# inputs:
# DIRECTNAME: name of directory where files interested are located at
# FILENAME_FORMAT: format of filename interested
# and outputs 
# FILENAME_MATCHED: full paths of matched files
#---------------------------------------------------------------------------

def get_path(directname, filename_format):
    
    # use os.listdir to get a list of filenames in directory provided
    file_list = os.listdir(directname)
    
    # loop through every filename in the list and check if it matches with
    # the search pattern (filename_format), if so, obtain full path of the file
    # with os.path.join and append info to the output list
    filename_matched = []   # initialize output list 
    for item in range(len(file_list)):   # for every filename in curr. dir.
        # obain match object 
        cur_filename = re.search(filename_format, file_list[item])
        # if matched object is found, obtain its full path
        if cur_filename:
            full_path = os.path.join(directname, cur_filename.string)
            filename_matched.append(full_path)   # append info to output list
    
    # format return statement
    return filename_matched

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

############################################################################

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

# SORT_ELECT FUNCTION:

# SORT_ELECT(ELECT_DF)
# input:
# ELECT_DF is the input dataframe with all electrodes known
# outut:
# ELE_SORTED_LIST is the output nested list with electrodes sorted

# With array containing coordinates of electrodes obtained, this function
# sorts electrodes according to their types and outputs nested list with 
# each entry storing data of those of the same type. e.g. given a set of 
# electrodes with names [LPI1, LPI2, ROF1, ROF2, LOF1, LOF2], the types of
# electrodes identified are LPI, ROF, and LOF. Hence, the output nested list
#  will contain one entry for each type identified, with all info. such as
# coordinates and name tags stored within itself.

# REMARKS:
# input electrode dataframe HAS to follow the format outlined below, 
# row 0 contains headers of columns
# column 1 contains name tags of eletrodes, 
# column 2:4 contain their x, y, and z-coordinates, respectively

def sort_elect(elect_df):
    
#---------------------------------------------------------------------------
# PART (I): clean electrode input array

# file containing info. about electrodes may contains rows with value of 0
# assigned to missing name tags and coordinates,
# remove such rows and return resultant array to electar_cleaned

    # initialize output variable
    elect_cleaned_list = []
    
    # use loop to go through each row of electar df , if name tag (col. 2) 
    # and coordinates (col. 3 to 5) are NOT missing, 
    # assign info. of that row to electar_cleaned list
    for i in range(elect_df.shape[0]):
        if bool(~pd.isnull(elect_df.iloc[i,1])) and \
            all(~pd.isnull(elect_df.iloc[i, 2:5])):
            elect_cleaned_list.append(elect_df.iloc[i])
    
    # convert nested list to df with original columns 
    elect_cleaned_df = pd.DataFrame(elect_cleaned_list, columns = elect_df.columns) 
            
# END PART (I): clean electrode input array
#---------------------------------------------------------------------------

# PART (II): sort cleaned electrode array

    # initialize lists, ele_sorted_list and cur_ele_list, the former is 
    # a nested list with each entry stores a list of electrodes of same type
    # the latter is used for identifying electrodes of same type
    ele_sorted_list = []
    cur_ele_list  = []
    
    # assign first row of cleaned elect df to cur_ele_list
    cur_ele_list.append(elect_cleaned_df.iloc[0, 1:].tolist())

    # each electrode name contains a compo. of string and a compo. of number
    # e.g. 'RPF1' has string compo. of 'RPF' and numeric compo. of '1'
    # use for loop to go through each row in cleaned elect df and use regular 
    # expression to search for matching alphabetic and numeric components in 
    # current and the next rows,
    # if names match and number of the next one is larger than the current
    # one, store the info. within current electrode type, 
    # else create another list for the new electrode type
    for i in range(elect_cleaned_df.shape[0] - 1):
        ele_name = elect_cleaned_df.iloc[i, 1]   # name tag of current electrode
        ele_name1 = elect_cleaned_df.iloc[i+1, 1]   # name tag of next electrode
        
        # use re to find matching pattern, alphabetic and numeric
        alph_compo = re.search('[a-zA-Z]+', ele_name).group()   # alphabetic compo. of current ele.
        num_compo = re.search('[0-9]+', ele_name).group()   # numeric compo. of current ele.
        alph_compo1 = re.search('[a-zA-Z]+', ele_name1).group()   # alphabetic compo. of next ele.
        num_compo1 = re.search('[0-9]+', ele_name1).group()   # numeric compo. of next ele.
        
        # check if alph. compo. agree and if number of curr. ele. is smaller than 
        # the next one, if so, append current row to cur_ele_list, otherwise, 
        # create new list and assign info. of electrode to it
        if alph_compo == alph_compo1 and num_compo < num_compo1:
            cur_ele_list.append(elect_cleaned_df.iloc[i+1, 1:].tolist())
        else:
            # if any compo. doesn't match, it means a new type of electrode is found
            # add that row to a new list
            ele_sorted_list.append(cur_ele_list)   # update ele_sorted_list
            cur_ele_list = []   # update content of curr. list
            cur_ele_list.append(elect_cleaned_df.iloc[i+1, 1:].tolist())   # assign info.
    
    # after finished looping through all rows in cleaned ele. df, append the 
    # 'final' current list to ele_sorted_list
    ele_sorted_list.append(cur_ele_list)
    
    # print msg. after sorting showing how many types of electrodes are found
    print('electrodes sorted, {} of them found'.format(len(ele_sorted_list))) 
    
    # assign sorted ele. list to return statement 
    return ele_sorted_list

# END PART (II): sort cleaned electrode array
# ---------------------------------------------------------------------------

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

############################################################################

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------