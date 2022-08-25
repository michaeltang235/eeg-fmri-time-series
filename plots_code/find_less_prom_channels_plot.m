close all
clear all

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
% BEGIN USER INPUT

% enter subject number
subnum = '27';

% enter path where input struct. is stored at
directname = ['/work/levan_lab/mtang/fmri_project/', 'sub', subnum];
% directname = ['/work/levan_lab/mtang/fmri_project/matrices', filesep 'static_analysis_reho'];   % direct. of output matrix

% format directory and filename of input files
directname_input = [directname filesep 'eeg_data'];
filename_input = 'find_less_prom_channels.mat'; 

% enter path where figures are saved at
directname_op = [directname filesep 'plots' filesep 'find_less_prom_channels'];

% enter if output needs to be saved at path indicated abve
op_results = 1;

% END USER INPUT
%---------------------------------------------------------------------------
%---------------------------------------------------------------------------

% load structure
st = load(fullfile(directname_input, filename_input));

% get session numbers
fd_run = fieldnames(st.opst);

% make plot for every session
for run_ind = 1:numel(fd_run)
% run_ind = 1;

% get arrays storing mean time series of all prom. ch. of curr. ev. type
% and avg. waveform of every channel w.r.t. each event type
mean_eeg_ts = st.opst.(fd_run{run_ind}).mean_eeg_ts;  
avg_waveform_ch = st.opst.(fd_run{run_ind}).avg_waveform_ch;   

% access info for every event type under current run
for ch_type_ind = 1:numel(mean_eeg_ts)
    
% get event def. of curr. run
event_def_curr = st.opst.(fd_run{run_ind}).event_def_prom_ch{ch_type_ind};

% get number of target chnnaels under curr. event type and curr. session
ch_ind = 1:size(avg_waveform_ch{ch_type_ind}(:, 1), 1);

% format time axis  for plotting
t_axis = 1:1:length(mean_eeg_ts{ch_type_ind});

% format of the following arrays are given below. 
% ts_mat: top row stores mean time series of all prom. ch. of curr. event
% type
% ------- other rows stores avg. time series of every target channels
% max_mat: max. of every time series in each row
% spac_mat: spacing matrix to be added to ts_mat to lift every time series
% from the previous one for plotting
% all_ts_mat: normalized time series - spac_mat; (normalized time series
% means that the amplitude is within [0, 1], for plotting purposes only
ts_mat = [mean_eeg_ts{ch_type_ind}; cell2mat(avg_waveform_ch{ch_type_ind}(ch_ind, end))];
max_mat = repmat(max(ts_mat, [], 2), [1, size(ts_mat, 2)]);  
spac_mat = repmat([0:1:size(ts_mat, 1)-1]', [1, size(ts_mat, 2)]);
all_ts_mat = ts_mat./max_mat - spac_mat;

% open figure window
f_grand = figure('units','normalized','outerposition',[0 0 0.5 1]);

% plot mean time series of all prom. ch. of curr. event type, denoted in black 
% then plot avg. waveform of other channels
plot(t_axis, all_ts_mat(1, :), 'k'); hold on
for i = 2:size(all_ts_mat, 1)
plot(t_axis, all_ts_mat(i, :)); hold on
end

% get current axes
ax = gca;

% format x-axis limits
xmin = min(t_axis);
xmax = max(t_axis);
ax.XLim = [xmin, xmax];

% format y-axis limits
ax.YLim = [-(size(spac_mat, 1)-1), 2];

% yticks = [0:1:size(ts_mat, 1)-1];
% yticklabels = ['mean'; avg_waveform_ch{ch_type_ind}(ch_ind, 1)];

% note: YTICKS AND YTICK LABELS HAVE TO BE IN THE SAME ORDER, 
% want to place mean eeg time series on top, then add each avg. waveform of
% other channels below
% get channel name of every time series plotted
ch_name_subset = avg_waveform_ch{ch_type_ind}(ch_ind, 1);
yticks = [-size(ts_mat, 1)+1:1:0];   % format yticks
yticklabels = {};   % format ytick labels for target channels
for i = size(all_ts_mat, 1)-1:-1:1
    yticklabels = [yticklabels; ch_name_subset{i}];
end
yticklabels = [yticklabels; 'mean wf-' num2str(event_def_curr)];   % add mean eeg ts to yticklabels

% format yticks and ytick labels
ax.YTick = yticks;
ax.YTickLabel = yticklabels;

% set tick label interpreter
ax.TickLabelInterpreter = 'latex';

% set aspect ratio of plot
pbr = 1.5;
pbaspect([xmax-xmin, pbr*(xmax-xmin), 1]); % multiple y-axis by the factor

% set title for curr. plot
title(['sub' subnum, '\_', fd_run{run_ind}]);
ax.Title.Interpreter = 'latex';

% format event def. for output filename 
ev_def_op = st.opst.(fd_run{run_ind}).prom_ch_type{ch_type_ind}{1, 1};

% format filename for current plot
filename_op = sprintf('sub%s_%s_ev_def_%s', subnum, fd_run{run_ind}, num2str(ev_def_op));
            
% save figure to file 
if op_results == 1
    saveas(f_grand, fullfile(directname_op, filename_op), 'epsc');
    saveas(f_grand, fullfile(directname_op, filename_op), 'jpeg');
end


end   % end for ch_type_ind = 1:numel(mean_eeg_ts)

end   % end for run_ind = 1:numel(fd_run)
