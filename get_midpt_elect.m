% This function intakes array of sorted electrode contacts and nifti image
% that the midpoints of electrode contacts are to be converted to, and
% generate output array which lists midpoint of each electrode channel in
% both mni and image space.

% OP_MIDPT = GET_MIDPT(ELE_SORTED, NIFTI_IMG_INFO)
% ELE_SORTED is the array of sorted electrode contacts 
% NIFTI_IMG_INFO is the metadata (info) of nifti image, in which the
% transformation matrix is used to transform coordinates from mni space to
% image space.
% OP_MIDPT is the output array of midpint, with the format listed below
% 1st layer, cell for each type of electrodes, 
% e.g. op_midpt{1} denotes cell for the first type of electrodes identified
% 2nd layer, midpoint of each pair of electrode contacts of the type
% specified, with
% col. 1 = name of channel (formed by joining two electrode contacts)
% col. 2 = midpoint between the two contacts in mni space
% col. 3 = midpoint between the two contacts in image space

function op_midpt = get_midpt_elect(ele_sorted, nifti_img_info)

%---------------------------------------------------------------------------
% Part (I): find midpoint of each pair of electrode contacts

% initialize midpt cell array
midpt = {};

% use loop to go through each entry in ele_sorted
% get midpoint of each pair of ele. contacts and store value to array
for i = 1:numel(ele_sorted)
    for j = 1:size(ele_sorted{i}, 1) - 1
        midpt{i}{j, 1} = [ele_sorted{i}{j, 1} '-' ele_sorted{i}{j+1, 1}];   % string showing the ele. pair
        mx = 0.5*(ele_sorted{i}{j, 2} + ele_sorted{i}{j+1, 2});   % midpt. in x
        my = 0.5*(ele_sorted{i}{j, 3} + ele_sorted{i}{j+1, 3});   % midpt. in y
        mz = 0.5*(ele_sorted{i}{j, 4} + ele_sorted{i}{j+1, 4});   % midpt. in z
        midpt{i}{j, 2} = [mx, my, mz];   % assign value to array
    end
end

% END Part (I): find midpoint of each pair of electrode contacts
%---------------------------------------------------------------------------
% Part (II): convert coordinates from mni space to image space

for i = 1:numel(ele_sorted)   % each electrode
    for j = 1:size(midpt{i}, 1)   % each pair of electrode contacts
        % convert from mni to image space
        midpt{i}{j, 3} = round(mni2ijk(midpt{i}{j, 2}, nifti_img_info)); 
    end
end

% END Part (II): convert coordinates from mni space to image space
%---------------------------------------------------------------------------
% Part (III): output variables calculated

op_midpt = midpt;

% END Part (III): output variables calculated
%---------------------------------------------------------------------------

end