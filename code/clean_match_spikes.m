% This function first cleans spike rates array (A), then match its channels 
% to that of another array (B), usually B could be a sorted t-statistics
% array, sorted signals array, etc

% SPIKES = CLEAN_MATCH_SPIKES(SPIKES_AR, REF_AR, REF_AR_ALL)
% SPIKES_AR is spike rate array 
% REF_AR is the reference array to which channels of spike rates are
% matched
% REF_AR_ALL is the single layer of reference array
% SPIKES is a structure that contains the following
% (i) cleaned spike rates array
% (ii) spike rates array with channels matched with ref. array
% (iii) matched spike rates in single layer with ids same as ref. array
% assigned

function spikes = clean_match_spikes(spikes_ar, ref_ar, ref_ar_all)

%---------------------------------------------------------------------------
% PART (I): clean input spike rates array

% sometimes the file contains extra rows that aren't showing any electrode
% pairs, i.e. rows showing the name tags (e.g. LO1, L02, ...) instead of
% (L01 - L02), 
% select only rows that contain info. about electrode pairs (e.g. L01-L02)

% get number of rows in spikerates array
srni = size(spikes_ar, 1);   % length in i dimension, spike rates

% create empty array named spikerates_cleaned
spikesar_cleaned = {};

% assign the first row (header info., i.e. Channel, spike/min, ...)
% of input file to spikerates_cleaned
spikesar_cleaned(1,:) = spikes_ar(1,:);

% set row index as 1 for iteration in loop
rowind = 1;

% search str '-' in each row of input spike rate file, if found, assign all
% columns of that row to spikesar_cleaned array, then increment rownum by
% 1 and continue the loop
for i = 1:srni   % for each row
    if ~isempty(regexp(spikes_ar{i,1}, '-'))   % if str '-' is found in certain name tag (first col.)
        spikesar_cleaned(rowind,:) = spikes_ar(i, :);   % assign all info. to array
        rowind = rowind + 1;   % increment rownum by 1
    end
end

% END PART (I): clean input spike rates array
%---------------------------------------------------------------------------
% PART (II): match names of channels in cleaned spike rates array to 
% that of second input array (ref_array)

% use second input array as reference, 
% define search str as each name tag in ref. array
% find corresponding row number in cleaned spike rates array,  
% if found, assign related info. to array,
% if not, assign value of 0 to that specific name tag
% note: assume ref. array has NO MISSING electrode contacts

% initialize spike_sorted cell array
spikes_matched = {};

for i = 1:numel(ref_ar)   % for each type of electrode
    for j = 1:length(ref_ar{i}(:,1))   % for each pair of electrode contacts
        search_str = ref_ar{i}{j,1};   % search str accor. to name tag in sig_box
        % get index (row number) of name tag in spikesar_cleaned
        if ~isempty(find(strcmp(search_str, spikesar_cleaned(:,1))))
            indreq = find(strcmp(search_str, spikesar_cleaned(:,1)));   % row index required
            % assign values to spike_matched
            spikes_matched{i}(j,:) = spikesar_cleaned(indreq,:); 
        else
            % if specific electrode contacts (name tag) isn't found in
            % spikesar_cleaned array, assign value of 0 to that name
            spikes_matched{i}(j,:) = {ref_ar{i}{j,1} [0] [0]};   
            sprintf('i=%d, j=%d, not found, but added 0 to name tag', i, j);           
        end
    end
end

% END PART (II): match names of channels in cleaned spike rates array to 
% that of second input array (ref_array)
%---------------------------------------------------------------------------
% PART (III): place spike rates in single layer, with ids matched with
% ref_ar_all

% ONLY EXECUTE THIS PART IF 3RD INPUT ARG. IS NOT EMPTY
if ~isempty(ref_ar_all)
% in ref_ar_all (3rd input arg.) (A), data are assigned in single layer, 
% with id assigned to each of them.
% place matched spike rates in single layer (B) and assign the same id of 
% each channel from (A) to the array (B)

% initialize array listing spike rates of all channels in single layer,
% with the following format
% col. 1 = id of channel
% col. 2 = name tag of channel
% col. 3 = spikes/min
% col. 4 = spike+HFO/min
spikes_all = {};
rownum = 1;   % initialize rownum as 1 

for type = 1:numel(spikes_matched)   % for each type of channel in spike rates array
    for i = 1:size(spikes_matched{type}, 1)   % for each channel in current type
        spikes_name_tag = spikes_matched{type}{i, 1};   % name tag (str) of current channel
        % compare strings in the two arrays, 1 for identical strings, 
        % 0 otherwise, then use find to get row index of the only non-zero entry
        row_ind = find(strcmp(spikes_name_tag, ref_ar_all(:, 2)));   % row index required
        spikes_all{rownum, 1} = ref_ar_all{row_ind, 1};   % id of current channel
        spikes_all{rownum, 2} = ref_ar_all{row_ind, 2};   % name of current channel
        spikes_all{rownum, 3} = spikes_matched{type}{i, 2};   % spike/min of current channel
        spikes_all{rownum, 4} = spikes_matched{type}{i, 3};   % spike+HFO/min of current channel
        rownum = rownum + 1;   % increment row number by 1 for next channel
    end
end

end   % end if isempty(ref_ar_all)

% END PART (III): place spike rates in single layer, with ids matched with
% ref_ar_all
%---------------------------------------------------------------------------

% PART (IV): store quantities calcu. in spikes struct.

spikes = struct;   % initialize structure
spikes.spikes_cleaned = spikesar_cleaned;   % store cleaned spike rates array
spikes.spikes_matched = spikes_matched;   % store matched spike rates array
spikes.spikes_all = spikes_all;   % store matched spike rates array in the form of single layer 

% END PART (IV): store quantities calcu. in spikes struct.
%---------------------------------------------------------------------------

