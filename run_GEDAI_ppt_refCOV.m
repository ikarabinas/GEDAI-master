% Set parameters for GEDAI, including the file path to a previously computed 
% individual participant leadfield-derived refCOV, and run script to run GEDAI on file

datadir = '/athena/grosenicklab/scratch/imk2003/acc_tmseeg/eeg_data/RELAX_GEDAI/RELAX_twICA_GEDAI-dlpfc/';
refCOV_path = '/athena/grosenicklab/scratch/imk2003/eeg_sources_data/GEDAI_refCOV';

%% Load data files needed for plotting
filename = 'm155_dlpfc_day3_reststate1_post_20230315_123516.set';

% Load raw .set file
EEG = pop_loadset(fullfile(datadir, filename));

% Load participant-specific refCOV
refCOV = match_and_load_file(refCOV_path, filename, "refCOV");

% Run GEDAI with the custom refCOV and otherwise default params
% Retain the cleaned continous EEG and artifact outputs for later use
[EEGclean, EEGartifacts] = GEDAI(EEG, 'auto', 12, 0.5, refCOV, true, false);

% Save figure
fig = gcf;
saveas(fig, fullfile(datadir, 'SENSAI_plots', sprintf('SENSAI_vis_%s.png', filename)));
fprintf('Saving SENSAI plot fig to %s/SENSAI_plots', datadir);

%%____________________________________________________________
function datafile = match_and_load_file(dirpath, filename, filetype)

if filetype == "refCOV"
    prefix_match_number = 2;
    file_suffix = '*.mat';
else
    prefix_match_number = 5;
    file_suffix = '*.set';
end

% Build prefix from first 5 underscore-delimited parts of EEG filename
[~, eeg_filename_noext, ~] = fileparts(filename);
eeg_parts = strsplit(eeg_filename_noext, '_');
eeg_prefix = strjoin(eeg_parts(1:prefix_match_number), '_');

% Load file
dirList = dir(fullfile(dirpath, file_suffix));
matchIdx = find(startsWith({dirList.name}, eeg_prefix));

if isempty(matchIdx)
    error('No matching file found in %s for %s', dirpath, filetype);
elseif length(matchIdx) > 1
    error('Multiple matching files found in %s for %s', dirpath, filetype);
end

% Define datapath
datapath = fullfile(dirpath, dirList(matchIdx).name);

if filetype == "refCOV"
    disp(['Loading custom ppt refCOV from: ' datapath]);
    data = load(datapath);
    datafile = data.refCOV;
else
    disp(['Loading .set file: ' datapath]);
    datafile = pop_loadset(datapath);
end

end