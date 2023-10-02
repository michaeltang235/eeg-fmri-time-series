clear all
close all

% overlaying map of falff on anatomical image across different time segments,
% with number of IEDs per min. as the title

% falff of channels recording IED type interested is denoted in
% filled markers, while that of channels having other or no IED types is
% denoted in unfilled markers.

% the 4 consecutive time segments were selected by users manually
% user might have to manually adjust z-level at which anatomical image is
% selected so that dots of falff allign well with tracks of electrodes 

% users may have to manually adjust (around line 319) z-level (zlv_anat_ijk) 
% at which anatomical image is selected (which is done by .nii header, but not so perfect) 
% for max. alignment with tracks left by electrodes implanted.
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT
% enter subject number (str, e.g. sub31, sub27, ...)
subnum = '34';

% enter run index interested
run_ind = 1;

% enter event index interested
ev_ind = 1;

% enter index of 1st time segment interested (it forms 4 consecutive
% time segments using the index entered, e.g. seg_ind_start = 1 generates
% the following sequence: 1, 2, 3, 4)
seg_ind_start = 10;

% enter path to directory where all input files are located
% directname = ['/work/levan_lab/mtang/fmri_project/'];
directname = 'C:\Users\siumichael.tang\Downloads\fmri_project';

% enter if output should be written to path (1 = yes, 0 = no)
op_results = 1;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% format path to directory and name(s) of input file(s)
% anat. img. of subject interested
fname_anat = [directname, filesep, 'sub', subnum, filesep, 'anat_img_new_norm'];   % sub-directory where anat. img. is
filename_anat = 'w3*.nii';   % pattern of anat. img. filename

% dynamic falff vs during-scan ied .mat file (falff of channels with ieds
% recorded during scan)
fname_falff = [directname, filesep, 'matrices', filesep, 'dynamic_analysis_falff'];  % directory of file(s) below
filename_name_falff_mat = 'dynamic_analysis_falff.mat';

% dynamic falff map .mat file (falff of all channels)
filename_name_falff_map = 'dynamic_analysis_falff_map.mat';

% obtain paths required in directories entered above
[~, anat_img_path, ~] = get_path(fname_anat, filename_anat);   % path of anat. img.
[~, falff_mat_path, ~] = get_path(fname_falff, filename_name_falff_mat);   % path of falff vs during-scan ied .mat file
[~, falff_map_path, ~] = get_path(fname_falff, filename_name_falff_map);   % path of falff map .mat file

% load files using paths found
anat_img = niftiread(anat_img_path{1});   % anat. img., using first available path found
falff_mat = load(falff_mat_path{1});   % only one .mat file is here, so select the first one
falff_map = load(falff_map_path{1});   % only one .mat file is here, so select the first one

% acess falff across time segments in subejct interested

% obtain list of subjects in falff .mat input file, then get index of the 
% subject interested in the list
sub_list = fieldnames(falff_mat.grand);
sub_ind_req = find(strcmp(sub_list, ['sub', subnum]));

% obtain list of runs in the subject interested
run_list_req = fieldnames(falff_mat.grand.(sub_list{sub_ind_req}));

% obtain list of subjects in falff map input file, then get index of the 
% subject interested in the list
sub_list_map = fieldnames(falff_map.grand);
sub_ind_req_map = find(strcmp(sub_list_map, ['sub', subnum]));

% access the following usig subject number and run index input
% (1) mid point array of all channels
% (2) channel ch_ext_ds_ied_falff array, which has the format given below,
% 1st col. = subject number
% 2nd col. = ied type 
% 3rd col. = channel name
% 4th col. = during scan ied rate per minute
% 5th col. = time series feature (e.g. ReHo, falff)
% (3) falff map of all channels in subject and run interested, which is
% stored inside sig_box_all array, with format givene below
% 1st col. = channel name
% 2nd col. = falff in every time segment

% falff_mat_req = falff_mat.grand.(sub_list{sub_ind_req}).(run_list_req{run_ind});
midpt_all = falff_mat.grand.(sub_list{sub_ind_req}).(run_list_req{run_ind}).midpt_all;
ch_ext_ds_ied_falff = falff_mat.grand.(sub_list{sub_ind_req}).(run_list_req{run_ind}).ch_ext_ds_ied_falff;
sig_box_all = falff_map.grand.(sub_list{sub_ind_req}).(run_list_req{run_ind}).sig_box_all;

% get unique event types from ch_ext_ds_ied_falff (channel-specific falff across time)
ev_col = [];   % event column 2 of ch_ext_ds_ied_falff array
for i = 1:size(ch_ext_ds_ied_falff)
    ev_col = [ev_col; str2num(ch_ext_ds_ied_falff{i, 2})];
end
uni_ev_type = unique(ev_col);   % unique event type

% sort ch_ext_ds_ied_falff array by event type
% initialize array with format given below,
% 1st layer = unique event type
% 2nd layer = same as ch_ext_ds_ied_falff, but for the uni. ev. type
ch_ext_ds_ied_falff_sorted = {};
ncol = size(ch_ext_ds_ied_falff, 2);   % get number of col. est.
row_num = 1;   % initialize row number 
for item = 1:numel(uni_ev_type)   % for every unique event type
    cur_uni_ev_type = uni_ev_type(item);   % current unique event type
    for i = 1:size(ch_ext_ds_ied_falff, 1)   % for every row in channel array
        cur_ev_type = ch_ext_ds_ied_falff{i, 2};   % ev. type from col. 2
        if isequal(str2num(cur_ev_type), cur_uni_ev_type)   % if event types matched
            ch_ext_ds_ied_falff_sorted{item}(row_num, 1:ncol) = ch_ext_ds_ied_falff(i, :);   % obtain all info.
            row_num = row_num + 1;   % set row number for next channel
        end
    end
    row_num = 1;   % set row number for next event type
end

% get midpt coordinates of channels, append midpt to rightmost col. of
% sorted channels array
% initialize array, ch_ext_ds_ied_falff_midpt, with format given below
% same format as ch_ext_ds_ied_falff_sorted, with midpt appeneded to
% rightmost col.
ch_ext_ds_ied_falff_midpt = ch_ext_ds_ied_falff_sorted;
for item = 1:numel(ch_ext_ds_ied_falff_sorted)  % for each uni. ev. type
    ncol = size(ch_ext_ds_ied_falff_sorted{item}, 2);   % get num. of col. est.
    for i = 1:size(ch_ext_ds_ied_falff_sorted{item}, 1)   % for every channel under curr. uni. ev. type
        cur_ch_name = ch_ext_ds_ied_falff_sorted{item}{i, 3};   % curr. ch. name
        [~, row_req, ~] = intersect(midpt_all(:, 1), cur_ch_name);   % get req. index in midpt_all
        % append midpt coord. to rightmost col. of array
        ch_ext_ds_ied_falff_midpt{item}{i, ncol+1} = midpt_all{row_req, end};
    end
end

% get most frequently occured z-level (mode) from all z-coord. obtained for
% every event type. If all z-levels have the same freq., then mode returns
% the first available z-level in z_coord
% initialize most_freq_z array, with format given below,
% 1st layer = unique event type
% 2nd layer has format given below,
% 1st col. = event type
% 2nd col. = most freq. occurred z-level
most_freq_z = cell(1, numel(ch_ext_ds_ied_falff_midpt));
for item = 1:numel(ch_ext_ds_ied_falff_midpt)   % for every event type
    most_freq_z{item}{1} = str2num(ch_ext_ds_ied_falff_midpt{item}{1, 2});   % assign event type from col. 2
    ch_z_coord = [];   % initialize array storing z-coordinate of all channels under curr. ev. type
    for ch = 1:size(ch_ext_ds_ied_falff_midpt{item}, 1)
        ch_z_coord = [ch_z_coord; ch_ext_ds_ied_falff_midpt{item}{ch, end}(3)];   % obtain z-level from 3rd entry in rightmost col.
    end
    most_freq_z{item}{2} = mode(ch_z_coord);   % use mode to get most freq. occurred z-level
end

% select channels with z-level falling within most freq. z +/- 1 for every event type
% initialize ch_ds_ied_falff_sel array with same format as ch_onset_ds_ied_falff_midpt
ch_ds_ied_falff_sel = cell(1, numel(most_freq_z));
for i = 1:numel(most_freq_z)   % for every unique event type
    cur_uni_ev_type = most_freq_z{i}{1};   % obtain curr. uni. ev. type from 1st entry of curr. cell in most_freq_z
    row_num = 1;   % initializ row number
    for j = 1:numel(ch_ext_ds_ied_falff_midpt)   % for every event type (cell) in channel array
        % obtain curr. ev. type from col. 2 of channel array 
        cur_ch_ev_type = str2num(ch_ext_ds_ied_falff_midpt{j}{1, 2});
        if isequal(cur_uni_ev_type, cur_ch_ev_type)   % if event types in both arrays matched
            for k = 1:size(ch_ext_ds_ied_falff_midpt{j}, 1)   % for every channel in curr. ev. type of channel array
                cur_z = ch_ext_ds_ied_falff_midpt{j}{k, end}(3);   % obtain z-level of curr. channel 
                cur_most_freq_z = most_freq_z{i}{end};   % obtain most-freq. z-level of curr. event type
                if cur_z >= cur_most_freq_z-1 && cur_z <= cur_most_freq_z+1   % if curr. z-level falls into most freq. z +/- 1
                    % concatenate info to curr. cell of event type
                    ch_ds_ied_falff_sel{i}(row_num, :) = ch_ext_ds_ied_falff_midpt{j}(k, :);
                    row_num = row_num + 1;   % increment row number by 1 for next channel with same z-level 
                end
            end
        end
    end
end

% obtain a list of ieds in time segments interested from 1st channel in the
% array (ch_ds_ied_falff_sel) obtained above
ds_ied_list = ch_ds_ied_falff_sel{ev_ind}{1, 4}(seg_ind_start:seg_ind_start+3);   % 4th col. of ch_ds_ied_falff_sel

% obtain a list of falff in time segments interested for every channel selected
% under the event type interested
% initialize falff_list_sel with format given below
% 1st col. = channel name
% 2nd col. = falff in time segments interested
falff_list_sel = {};
for row = 1:size(ch_ds_ied_falff_sel{ev_ind}, 1)
    falff_list_sel{row, 1} = ch_ds_ied_falff_sel{ev_ind}{row, 3};   % 3rd col. of ch_ds_ied_falff_sel
    falff_list_sel{row, 2} = cell2mat(ch_ds_ied_falff_sel{ev_ind}{row, 5}(seg_ind_start:seg_ind_start+3));   % 5th col. of ch_ds_Ied_falff_sel
end

%---------------------
% % make red-yellow-white color matrix
% nrow_rw = 11;   % num. of rows from red to yellow (regular)
% ryw_cmap = zeros(nrow_rw*2 - 1, 3);   % specify dimensions of color matrix
% 
% ryw_cmap(1, :) = [1 0 0];   % assign rgb values of red to first row 
% ryw_cmap(nrow_rw, :) = [1 1 0];   % assign rgb values of yellow to middle row
% ryw_cmap(end, :) = [1 1 1];   % assign rgb values of white to last row
% 
% step_size = 1/(nrow_rw-1);   % step size from red to yellow
% for i = 2:nrow_rw-1  % from red to yellow
%     ryw_cmap(i, :) = [1, ryw_cmap(i-1, 2) + step_size, 0];
% end
% 
% step_size = 1/(nrow_rw-1);   % step size from yellow to white
% for i = nrow_rw+1:nrow_rw*2 - 1  % from yellow to white
%     ryw_cmap(i, :) = [1, 1, ryw_cmap(i-1, 3) + step_size];
% end
% 
% %------------------------------------------------

% % %-----------------------------------------------------
% % OPTION 6: 21 rows in color matrix, (dark red, red, yellow, white)
% 
% % make dark red-red-yellow-white color matrix
nrow_rw = 11;   % num. of rows from red to yellow
ryw_cmap = zeros(nrow_rw*2 - 1, 3);   % specify dimensions of color matrix

% ryw_cmap(1, :) = [0.5 0 0 ];   % assign rgb values of dark red to first row
ryw_cmap(1, :) = [0.1 0 0 ];   % assign rgb values of dark red to first row
ryw_cmap(floor((nrow_rw*2-1)*1/3), :) = [1 0 0];   % assign rgb values of red to first row
ryw_cmap(floor((nrow_rw*2-1)*2/3), :) = [1 1 0];   % assign rgb values of yellow to middle row
ryw_cmap(end, :) = [1 1 1];   % assign rgb values of white to last row

step_size = 0.7/(floor((nrow_rw*2-1)*1/3) -1);   % step size from dark red to red
for i = 2:floor((nrow_rw*2-1)*1/3)-1  % from dark red to red
    ryw_cmap(i, :) = [ryw_cmap(i-1, 1) + step_size, 0, 0];
end

step_size = 1/(floor((nrow_rw*2-1)*2/3) - floor((nrow_rw*2-1)*1/3));   % step size from red to yellow
for i = floor((nrow_rw*2-1)*1/3)+1:floor((nrow_rw*2-1)*2/3)-1  % from red to yellow
    ryw_cmap(i, :) = [1, ryw_cmap(i-1, 2) + step_size, 0];
end

step_size = 1/(nrow_rw*2 - 1 - floor((nrow_rw*2-1)*2/3));  % step size from yellow to white
for i = floor((nrow_rw*2-1)*2/3)+1:nrow_rw*2 - 2  % from yellow to white
    ryw_cmap(i, :) = [1, 1, ryw_cmap(i-1, 3) + step_size];
end
% 
% %------------------------------------------------

% % %-----------------------------------------------------
% % OPTION 6: 21 rows in color matrix, (dark red, red, yellow, white)
% % (small portion of red)
% 
% % make dark red-red-yellow-white color matrix
% nrow_rw = 11;   % num. of rows from red to yellow
% ryw_cmap = zeros(nrow_rw*2 - 1, 3);   % specify dimensions of color matrix
% 
% ryw_cmap(1, :) = [0.5 0 0 ];   % assign rgb values of dark red to first row
% ryw_cmap(floor((nrow_rw*2-1)*1/6), :) = [1 0 0];   % assign rgb values of red to first row
% ryw_cmap(floor((nrow_rw*2-1)*2/6), :) = [1 1 0];   % assign rgb values of yellow to middle row
% ryw_cmap(end, :) = [1 1 1];   % assign rgb values of white to last row
% 
% step_size = 0.5/(floor((nrow_rw*2-1)*1/6) -1);   % step size from dark red to red
% for i = 2:floor((nrow_rw*2-1)*1/6)-1  % from dark red to red
%     ryw_cmap(i, :) = [ryw_cmap(i-1, 1) + step_size, 0, 0];
% end
% 
% step_size = 1/(floor((nrow_rw*2-1)*2/6) - floor((nrow_rw*2-1)*1/6));   % step size from red to yellow
% for i = floor((nrow_rw*2-1)*1/6)+1:floor((nrow_rw*2-1)*2/6)-1  % from red to yellow
%     ryw_cmap(i, :) = [1, ryw_cmap(i-1, 2) + step_size, 0];
% end
% 
% step_size = 1/(nrow_rw*2 - 1 - floor((nrow_rw*2-1)*2/6));  % step size from yellow to white
% for i = floor((nrow_rw*2-1)*2/6)+1:nrow_rw*2 - 2  % from yellow to white
%     ryw_cmap(i, :) = [1, 1, ryw_cmap(i-1, 3) + step_size];
% end
% 
% %------------------------------------------------

%-----------------------------------------------------------------------------
% make plot

% open figure window for curr. event type
% f1{ev_ind} = figure;
% f1{ev_ind} = figure('units','normalized','outerposition',[0 0 1 0.5]);
f1{ev_ind} = figure('units','normalized','outerposition',[0 0 1 1]);

% get dimensions of func. images
ni = size(anat_img, 1);   % length in x-direct.
nj = size(anat_img, 2);   % length in y-direct.

% create x-and y-vectors for imagesc
xvec = [1, ni];
yvec = [1, nj];

% create 2d matrices in x and y for contour
x2d = repmat([1:ni]', 1, nj);
y2d = repmat([1:nj], ni, 1);

% overlay falff map onto anatomical image for each time segment selected
for pn_num = 1:4

% this requires overlaying falff map on mri image, which is done by
% plotting the two varaibles on separate axes, then merging the axes
% toegther

% begin plotting at current panel, store subplot object to handle
sp = subplot(1, 4, pn_num);

% assign current subplot axes to ax_fmri, axes for mri image
ax_mri = sp;

% with the list of channels falling within the range of most freq. z-level
% obtained, get equivalent z-level in anatomical image for the event type 
% interested using the first channel in the list
ch_name_req = ch_ds_ied_falff_sel{ev_ind}{1, 3};   % get channel name
ch_midpt_ijk = get_anat_ijk([directname, filesep, 'sub', subnum], subnum, ch_name_req);   % obtain channel coordinates in ijk space
zlv_anat_ijk = ch_midpt_ijk(3);   % obtain z-level of curr. channel coordinates

% zlv_anat_ijk = 22;   % use this for sub33_run_ind_2_ev_ind_3
zlv_anat_ijk = 15;   % use this for sub34_run_ind_1_ev_ind_1

% plot a slice of anatomical image at the most freq. z level using imagesc
% on mri image axes
imagesc(ax_mri, xvec, yvec, anat_img(:,:,zlv_anat_ijk)');

% initialize another axis for falff
ax_falff = axes;

%----------------------------
% plot falff of channels within the same z-range as anat. img. and not 
% recording the IEDs interested with unfilled markers

% plot falff of all other channels within the same z-range

for row = 1:size(sig_box_all, 1)   % for every channel in sig_box_all
    cur_ch_name = sig_box_all{row, 1};   % get current channel name

    % obtain channel image space coordinates in normalized anatomical image
    cur_ch_ijk = get_anat_ijk([directname, filesep, 'sub', subnum], subnum, cur_ch_name);

    % if z-coordinate of curr. channel falls within z-level +/- 1 of anat.
    % image, and channel name is not in the list of falff chanenls interestd, 
    % plot their falff in current time segment
    if cur_ch_ijk(3) >= zlv_anat_ijk-1 && cur_ch_ijk(3) <= zlv_anat_ijk+1 ...
            && ~ismember(cur_ch_name, falff_list_sel(:, 1))
        % plot falff of current channel in curr. time segment, unfilled
        % marker
        scatter(ax_falff, cur_ch_ijk(1), cur_ch_ijk(2), 30, ...
            sig_box_all{row, end}{seg_ind_start+pn_num-1}, 'LineWidth', 0.8); hold on
    end
end
%----------------------------

% plot falff of channels selected (those recording the ied rate shown on title)
% in current time segment with filled markers
for ch = 1:size(falff_list_sel, 1)   % for every channel selected
    ch_name = falff_list_sel{ch, 1};   % curr. channel name

    % get channel coordinates in normalized anat. img. space
    ch_ijk = get_anat_ijk([directname, filesep, 'sub', subnum], subnum, ch_name);
    xpos_falff = ch_ijk(1);   % 1st entry for x-coordinate
    ypos_falff = ch_ijk(2);   % 2nd entry for y-coordinates

    % plot falff of current channel in curr. time segment
    scatter(ax_falff, xpos_falff, ypos_falff, 30, falff_list_sel{ch, end}(pn_num), 'filled'); hold on
end

%----------------------------

% synchroniz limits of both axes
linkaxes([ax_mri, ax_falff])

% adjust positions of falff axes to that of current subplot panel
ax_falff.Position = sp.Position;

% set ticks on y and x-axes
ax_mri.XTick = [20, 137];   % set positions of x-ticks
ax_mri.XTickLabel = {'L', 'R'};   % set x-ticks str
ax_mri.XAxis.FontSize = 15;   % set x-axis font size
ax_mri.YTick = [];   % set y-ticks as empty

% set visibility to off, and empty x and y-ticks for falff axes
ax_falff.Visible = 'off';   % visibility off
ax_falff.XTick = [];   % empty x-ticks
ax_falff.YTick = [];   % empty y-ticks

% set colormap for each of the axes
colormap(ax_mri, 'gray')   % gray colormap for fmri axes
% colormap(ax_falff, 'hot')   % hot colormap for falff axes
colormap(ax_falff, ryw_cmap)   % hot colormap for falff axes3

% get limits of falff of all selected channels in all segments interested
falff_min = min(min(cell2mat(falff_list_sel(:, end))));
falff_max = max(max(cell2mat(falff_list_sel(:, end))));

% set color bar limits on ax_falff axes
caxis(ax_falff, [falff_min, falff_max])

% set box on
box on;

% set interpreter of tick labels on both axes 
ax_mri.TickLabelInterpreter = 'latex';

% set line width of axes (frame of plot)
ax_mri.LineWidth = 1;

% set y-axis direction of fmri axis as normal 
ax_mri.YDir = 'normal';

% set data aspect ratio as 1:1 (both axes have the same data units)
daspect(ax_mri, [1 1 1])
daspect(ax_falff, [1 1 1])

% get current number of ied (current segment)
num_ied_curr_pan = ds_ied_list(pn_num);

% format title 
title_str = ['\# IED = ', num2str(num_ied_curr_pan), ' [/min]'];

% set title on fmri axes
title(ax_mri, title_str, 'interpreter', 'latex', 'fontsize', 13);

% place and format colorbar after the 4th panel (segment) is plotted
if pn_num == 4
    cb_falff = colorbar(ax_falff);   % assign colorbar to falff axis
    cb_falff.Title.String = 'fALFF';   % set title
    cb_falff.Title.Interpreter = 'latex';   % set title interpreter
    cb_falff.Title.FontSize = 12;    % set title font size
%     cb_falff.Ticks = [ceil(falff_min*10)/10:0.1:floor(falff_max*10)/10];   % set ticks
    cb_falff.Ticks = [ceil(falff_min*10)/10:0.1:floor(falff_max*10)/10];   % set ticks
    cb_falff.TickLabelInterpreter = 'latex';   % set tick label interpreter
    cb_falff.FontSize = 10;
end

% assign current mri and falff axes to variables according to panel number
if pn_num == 1   % the 1st panel
    ax1_mri = ax_mri;
    ax1_falff = ax_falff;
end
    
if pn_num == 2   % the 2nd panel
    ax2_mri = ax_mri;
    ax2_falff = ax_falff;
end

if pn_num == 3   % the 3rd panel
    ax3_mri = ax_mri;
    ax3_falff = ax_falff;
end

if pn_num == 4   % the 4th panel
    ax4_mri = ax_mri;
    ax4_falff = ax_falff;
end

end   % end for pn_num = 1:4

% adjust panel positions
hdiff = 0.04/10;   % set horizontal diff. btw. panels

% set pannel positions
ax2_mri.Position(1) = ax1_mri.Position(1) + ax1_mri.Position(3) + hdiff;
ax2_falff.Position(1) = ax1_falff.Position(1) + ax1_falff.Position(3) + hdiff;

ax3_mri.Position(1) = ax2_mri.Position(1) + ax2_mri.Position(3) + hdiff;
ax3_falff.Position(1) = ax2_falff.Position(1) + ax2_falff.Position(3) + hdiff;

ax4_mri.Position(1) = ax3_mri.Position(1) + ax3_mri.Position(3) + hdiff;
ax4_falff.Position(1) = ax3_falff.Position(1) + ax3_falff.Position(3) + hdiff;

% set colorbar positions
cb_falff.Position(1) = ax4_mri.Position(1) + ax4_mri.Position(3) + 1.2*hdiff;
cb_falff.Position(3) = 0.75*cb_falff.Position(3);

% output figure to path
fname_fig = [directname, filesep, 'sub', subnum, filesep, 'plots', filesep, 'falff_map_anat_img', filesep, 'may09_2023'];
filename_fig = ['falff_anat_sub_', subnum, '_', 'run_ind_', num2str(run_ind), ...
   '_ev_ind_', num2str(ev_ind), '_seg_ind_start_', num2str(seg_ind_start), '_new_cb'];

if op_results == 1
saveas(f1{ev_ind}, fullfile(fname_fig, filename_fig), 'epsc');
saveas(f1{ev_ind}, fullfile(fname_fig, filename_fig), 'png');
end