% Use this script to iteratively extract GEDAI SENSAI metrics (SSSI, NSSI, LDA accuracy) from already-cleaned files,
% run SENSAI_visualization.m to plot 2D SENSAI scatter, and save the dataset's SENSAI metrics to csv.
% 
% Individual participant refCOV were previously computed and are loaded in here, but this script could be
% modified to work with the 'precomputed' or 'interpolated' GEDAI denoising methods as well.

% Load eeglab
eeglab;

datadir = '/athena/grosenicklab/scratch/imk2003/acc_tmseeg/eeg_data/RELAX_GEDAI/RELAX_GEDAI-dlpfc';
savepath = fullfile(datadir, 'RELAXProcessed');
cleaned_path = fullfile(savepath, 'Cleaned_Data');
artifacts_path = fullfile(savepath, 'GEDAI_artifacts');
refCOV_path = '/athena/grosenicklab/scratch/imk2003/eeg_sources_data/GEDAI_refCOV';
[~, preprocessing_label] = fileparts(datadir);  % save datadir folder name for SENSAI metrics csv name

%% Load data files needed for plotting
% Get list of all cleaned files
all_cleaned_files = dir(fullfile(cleaned_path, '*GEDAI*.set'));
fprintf('Total cleaned files found: %d\n', length(all_cleaned_files));

% Create new dir for SENSAI figs if it does not exist
SENSAI_figdir = fullfile(datadir, 'RELAXProcessed', 'SENSAI_plots');
%mkdir(SENSAI_figdir)

% Turn off auto fig display
set(0, 'DefaultFigureVisible', 'off')

% Initialize struct
SENSAI_lda_ssi = struct();

% Iterate over all cleaned files in datadir
for i = 1:length(all_cleaned_files)
    filename = all_cleaned_files(i).name;
    
    % Strip extra junk from filename to get format ppt_day_reststatex_recording
    [~, eeg_filename_noext, ~] = fileparts(filename);
    eeg_parts = strsplit(eeg_filename_noext, '_');
    filename = strjoin(eeg_parts(1:5), '_');
 
    % Load files
    % Load raw .set file
    EEG = match_and_load_file(datadir, filename, "raw");
    
    % Load clean EEG (post RELAX twICA and GEDAI denoising)
    EEGclean = match_and_load_file(cleaned_path, filename, "clean");
    
    % Load artifact EEG (saved from GEDAI processing)
    EEGartifacts  = match_and_load_file(artifacts_path, filename, "artifacts");
    
    % Load participant-specific refCOV
    refCOV = match_and_load_file(refCOV_path, filename, "refCOV");
    
    % Apply minimal preprocessing to raw file (as in GEDAI) and set relevant GEDAI defaults
    % Ensure double input 
    EEG.data=double(EEG.data);
    
    % Apply avg reference as preprocessing for raw file as is done in GEDAI.m
    EEGavRef = GEDAI_nonRankDeficientAveRef(EEG); % non rank-deficient average referencing (Makoto's plugin)
    
    broadband_epoch_size = 1;

    %% Manifold classicication and plotting -- copied over from GEDAI.m section from 04/2026 update
    % --- Manifold Classification (Broadband) BEFORE Cleaning ---
    % Uses 50% overlapping 1-second epochs for denser coverage in the scatter plot
        %fprintf('\nGenerating SENSAI Plot for Broadband data...\n');

    % 50% overlapping epoch parameters
    epoch_samples    = round(EEGavRef.srate * broadband_epoch_size); % 1-second window
    epoch_step       = floor(epoch_samples / 2);                     % 50% overlap
    pnts_original    = size(EEGavRef.data, 2);
    eeg_data_temp    = EEGavRef.data;

    % Pad so the last window is complete
    last_start = floor((pnts_original - epoch_samples) / epoch_step) * epoch_step + 1;
    last_end   = last_start + epoch_samples - 1;
    if last_end > pnts_original
        samples_to_pad = last_end - pnts_original;
        reflection_segment = eeg_data_temp(:, end - samples_to_pad + 1 : end);
        eeg_data_temp = [eeg_data_temp, fliplr(reflection_segment)];
    end

    % Compute starts for all 50%-overlapping windows
    num_epochs = floor((size(eeg_data_temp, 2) - epoch_samples) / epoch_step) + 1;
    COV_emp_array_before = cell(num_epochs, 1);
    for epo = 1:num_epochs
        i_start = (epo - 1) * epoch_step + 1;
        i_end   = i_start + epoch_samples - 1;
        COV_emp_array_before{epo} = cov(eeg_data_temp(:, i_start:i_end)');
    end


    % --- Manifold Classification (Broadband) AFTER Cleaning ---
    % Uses the same 50% overlapping 1-second epoch parameters as the BEFORE block
        % fprintf('\nGenerating SENSAI Plot for Final Reconstructed Data (Before/After)...\n');
    eeg_data_temp      = EEGclean.data;
    artifact_data_temp = EEGartifacts.data;

    % Pad each signal independently so the last window is complete
    pnts_clean = size(eeg_data_temp, 2);
    last_start_clean = floor((pnts_clean - epoch_samples) / epoch_step) * epoch_step + 1;
    last_end_clean   = last_start_clean + epoch_samples - 1;
    if last_end_clean > pnts_clean
        pad_clean = last_end_clean - pnts_clean;
        eeg_data_temp = [eeg_data_temp, fliplr(eeg_data_temp(:, end - pad_clean + 1 : end))];
    end

    pnts_art = size(artifact_data_temp, 2);
    last_start_art = floor((pnts_art - epoch_samples) / epoch_step) * epoch_step + 1;
    last_end_art   = last_start_art + epoch_samples - 1;
    if last_end_art > pnts_art
        pad_art = last_end_art - pnts_art;
        artifact_data_temp = [artifact_data_temp, fliplr(artifact_data_temp(:, end - pad_art + 1 : end))];
    end

    % 50% overlapping windows – reuse num_epochs from the BEFORE block
    num_epochs_after = floor((size(eeg_data_temp, 2) - epoch_samples) / epoch_step) + 1;
    COV_emp_array_after     = cell(num_epochs_after, 1);
    COV_emp_array_artifacts = cell(num_epochs_after, 1);
    for epo = 1:num_epochs_after
        i_start = (epo - 1) * epoch_step + 1;
        i_end   = i_start + epoch_samples - 1;
        COV_emp_array_after{epo}     = cov(eeg_data_temp(:, i_start:i_end)');
        COV_emp_array_artifacts{epo} = cov(artifact_data_temp(:, i_start:i_end)');
    end

    % Align BEFORE array to the same epoch count for a fair comparison
    if num_epochs_after < num_epochs
        COV_emp_array_before = COV_emp_array_before(1:num_epochs_after);
    end

    [mean_ssi_before, mean_lda_accuracy, mean_ssi_after, mean_ssi_artifacts] = SENSAI_visualization(refCOV, COV_emp_array_before, COV_emp_array_after, COV_emp_array_artifacts);
    
    % Save SENSAI cleaning metrics
    SENSAI_lda_ssi(i).filename = filename;
    SENSAI_lda_ssi(i).ssi_before = mean_ssi_before;
    SENSAI_lda_ssi(i).lda_accuracy = mean_lda_accuracy;
    SENSAI_lda_ssi(i).sssi = mean_ssi_after;
    SENSAI_lda_ssi(i).nssi = mean_ssi_artifacts;

    % Save figure
    fig = gcf;
    saveas(fig, fullfile(SENSAI_figdir, sprintf('SENSAI_vis_%s.png', filename)));
    fprintf('Saving SENSAI plot fig to: %s\n', SENSAI_figdir);
    close(fig)
end

% Save SENSAI metrics to file as .csv
csv_filename = sprintf('%s_SENSAI_Metrics.csv', preprocessing_label);
T = struct2table(SENSAI_lda_ssi);
writetable(T, fullfile(savepath, csv_filename));
fprintf('SENSAI metrics saved to file as %s.\n', csv_filename);


% Turn auto fig display back on
set(0, 'DefaultFigureVisible', 'on')

%%-------------------------------------------------------------------
function datafile = match_and_load_file(dirpath, filename, filetype)

if filetype == "refCOV"
    [~, eeg_filename_noext, ~] = fileparts(filename);
    eeg_parts = strsplit(eeg_filename_noext, '_');
    filename = strjoin(eeg_parts(1:2), '_');
    file_suffix = '*.mat';
else
    file_suffix = '*.set';
end

% Load file
dirList = dir(fullfile(dirpath, file_suffix));
matchIdx = find(startsWith({dirList.name}, filename));

if isempty(matchIdx)
    error('No matching file found in %s for %s. Filename: %s', dirpath, filetype, filename);
elseif length(matchIdx) > 1
    error('Multiple matching files found in %s for %s. Filename: %s', dirpath, filetype, filename);
end

% Define datapath
datapath = fullfile(dirpath, dirList(matchIdx).name);

if filetype == "refCOV"
    disp(['Loading custom ppt refCOV from: ' datapath]);
    data = load(datapath);
    datafile = data.refCOV;
else
    datafile = pop_loadset(datapath);
end

end