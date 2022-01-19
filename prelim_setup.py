
import os 
import re

import numpy as np
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

# MNI2IJK FUNCTION:
# MNI2IJK(MNICOORD_INPUT, IMG_HEADER):
# inputs: 
# MNICOORD = coordinates in mni space
# IMG_HEADER = object containing metadata of image
# outputs:
# IJK = coordinates in image space

# It converts coordinates from mni space (x, y, z) to image space (i, j, k) 
# using image headers (metadata associated with nifti image)
#---------------------------------------------------------------------------

def mni2ijk(mnicoord_input, img_header):
    
    # convert input mni coordinates to numpy array
    mnicoord = np.array(mnicoord_input)
    
    # method 3: get (i,j,k) from (x,y,z) by using a matrix that is constructed 
    # when general affine transformation took place from (i,j,k) to (x,y,z)

    # get the required parameters for the transformation matrix from nifti info
    srow_x = img_header['srow_x']
    srow_y = img_header['srow_y']
    srow_z = img_header['srow_z']
    
    # assemble the trans. matrix with quantities obtained (first 3 entries)
    srowmat = np.array([srow_x[:3], srow_y[:3], srow_z[:3]])
    
    # from nifti documentation, (x,y,z) is related to (i,j,k) by a matrix
    # equation in the form of Ax=B.
    
    # by comparison, matrix A is the trans. matrix, while matrix B is the
    # remaining terms in the equation given.
    
    # so, create matrix B
    matb3 = mnicoord - np.array([srow_x[-1], srow_y[-1], srow_z[-1]])
    
    # solve the linear system using np.linalg.solve and store solution to ijk3
    ijk3 = np.linalg.solve(srowmat, matb3)   # method 3 is used
    
    # convert ijk3 from numpy array back to list 
    ijk = ijk3.tolist()
   
    # format return statement
    return ijk

#--------------------------------------------------------------------------------
# References:

# nifti documentation on coordinate transformation:
# https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/qsform.html

# nifti documentation on rotation matrix:
# https://nifti.nimh.nih.gov/nifti-1/documentation/nifti1fields/nifti1fields_pages/quatern.html/view?searchterm=nifti1_io.c

# info on affine transformation (12 degrees of freedom):
# https://afni.nimh.nih.gov/sscc/staff/rwcox/ISMRM_2006/Syllabus%202006%20-%203340/files/B_05.pdf
   

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

# GET_MIDPT_ELECT FUNCTION:
    
# This function intakes array of sorted electrode contacts and nifti image
# that the midpoints of electrode contacts are to be converted to, and
# generate output array which lists midpoint of each electrode channel in
# both mni and image space.

# GET_MIDPT_ELECT(ELE_SORTED, NIFTI_IMG_INFO):
# inputs:
# ELE_SORTED is the nested list of sorted electrode contacts 
# NIFTI_IMG_INFO is the metadata (info) of nifti image, in which the
# transformation matrix is used to transform coordinates from mni space to
# image space.
# output:
# MIDPT is the output nested list of midpint, with the format given below
# 1st layer, entry for each type of electrodes, 
# e.g. midpt[0] denotes cell for the first type of electrodes identified
# 2nd layer, midpoint of each pair of electrode contacts of the type
# specified, with
# col. 1 = name of channel (formed by joining two electrode contacts)
# col. 2 = midpoint between the two contacts in mni space
# col. 3 = midpoint between the two contacts in image space


def get_midpt_elect(ele_sorted, nifti_img_info):
    
#---------------------------------------------------------------------------
# Part (I): find midpoint of each pair of electrode contacts

    # initialize list for storing midpoint of every channel (electrode pair)
    # in this nested list, every entry contains a list with format given below:
    # 1st entry = name tag of channel
    # 2nd entry = list of midpt corrdinates [x, y, z]
    midpt = []
    
    # scan across every row within every entry in sorted ele. list, 
    # then assemble name tag of channel, 
    # and compute midpoint of x-, y-, and z-coordinates
    for i in range(len(ele_sorted)):   # every entry in sorted ele. list
        cur_list = []   # initialize list for storing info. of curr. ele. type
        # within curr. entry in nested list, scan across every sub-entry, 
        # format name tag for channel and compute midpoint
        for j in range(len(ele_sorted[i]) - 1):
            pair_name = ele_sorted[i][j][0] + '-' + ele_sorted[i][j+1][0]
            mx = 0.5*(ele_sorted[i][j][1] + ele_sorted[i][j+1][1]) 
            my = 0.5*(ele_sorted[i][j][2] + ele_sorted[i][j+1][2]) 
            mz = 0.5*(ele_sorted[i][j][3] + ele_sorted[i][j+1][3])
        
            # append info. of every channel to curr. list (curr. ele. type)
            cur_list.append([pair_name, [mx, my, mz]])      
        
        # after looping through all elements in current entry (all channels in 
        # current type, append curr. list to nested list, midpt)
        midpt.append(cur_list)
        
# END Part (I): find midpoint of each pair of electrode contacts                        
#---------------------------------------------------------------------------
# Part (II): convert coordinates from mni space to image space

    for i in range(len(midpt)):   # every type of electrode
        for j in range(len(midpt[i])):   # each pair of electrode contacts
            # convert from mni to image space, then round number to nearest integer
            ijk = np.rint(mni2ijk(midpt[i][j][-1], nifti_img_info))
            midpt[i][j].append(ijk.tolist())   # append ijk to curr. entry in midpt list

    # assign midpt list to return statement
    return midpt        
                           
# Part (II): convert coordinates from mni space to image space
#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------

############################################################################

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# FILER_MOTION_FUNCTION:

# Given motion parameters, this function removes motion artifacts from fmri 
# time series by solving the system of linear euqations (A*beta = y), where A 
# is the motion parameter matrix defined below, y is the time series for each 
# voxel, and beta is the least-square fitted solutions to the system. After, 
# residuals of signal without motion artifacts are given by y - A*beta
    
# FILTER_MOTION(MO_PARA, SIGNAL):
# inputs:
# MO_PARA is numpy array of motion parameters
# SIGNAL is array of signal required to be motion-filtered, dim. of signal
# array is ni X nj X nk X nt, where ni, nj, nk denote dim. in x, y, and z
# directions, and nt denotes dim. in time domain.
# and outputs
# RESIDUALS, the time series of signal with motion artifacts removed, with
# dim. the same as SIGNAL

def filter_motion(mo_para, signal):
    
    # initialize output var., residual
    residuals  = []
    
    # The 24 regressors are the 6 motion time series, 
    # those same 6 time series shifted by one scan, 
    # and the squares of the previous 12 regressors

    # create 6 time series shifted by one scan, horiz. conca. with
    # motion para. read to form reg_shifted 
    mo_shifted = np.concatenate((np.tile(0, (1, 6)), mo_para[:-1, :]), axis=0)   # motion para. shifted by one scan
    reg_shifted = np.concatenate((mo_para, mo_shifted), axis=1)   # combine motion para. read with shifted para.

    # square the 12 regressors and add them all toegther to form reg_sq
    mo_sq = np.square(reg_shifted)   # square of the previous 12 motion regressors
    reg_sq = np.concatenate((reg_shifted, mo_sq), axis=1)   # add mo_sq to list of regressors
    
    # create constant term to form all regressors required
    mo_constant = np.ones((reg_sq.shape[0], 1))   # constant term
    reg_mo = np.concatenate((reg_sq, mo_constant), axis=1)   # assemble all motion regressors
    
    # initialize residuals array for storing motion-filtered signals
    residuals = np.zeros(signal.shape)
    
    # dimensions of input singal 
    ni = signal.shape[0]   # dim. in x
    nj = signal.shape[1]   # dim. in y
    nk = signal.shape[2]   # dim. in z
    
    # loop through input signal array, perform motion filtering in each voxel
    for i in range(ni):   # for each voxel in x
        for j in range(nj):   # for each voxel in y
            for k in range(nk):   # for each voxel in z
            
            # solve system of eq. (A*beta = y) for beta, where A is the matrix
            # of motion regressors and y is the time series, select 0th entry of 
            # output to get the least sq. soluton, beta
            
                beta = np.linalg.lstsq(reg_mo, signal[i, j, k, :], rcond=None)[0]
                
                # comput the residuals by subtracting fitted regressors
                # i.e. time series - motion regressor * beta (matrix multip.)
                residuals[i, j, k, :] = signal[i, j, k, :] - np.matmul(reg_mo, beta)
                
    return residuals
                


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------