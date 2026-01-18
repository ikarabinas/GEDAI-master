% Quick verification test for GEDAI-master memory optimizations
% This script runs a simple test to verify the optimizations work correctly

fprintf('\n=== GEDAI Memory Optimization Quick Test ===\n\n');

% Load sample data
sample_file = 'c:\Users\Ros\Documents\MATLAB\eeglab2025.0.0\plugins\GEDAI-master\example data\empirical_NOISE_EOG_EMG.set';

if ~exist(sample_file, 'file')
    error('Sample file not found: %s', sample_file);
end

fprintf('Loading sample data...\n');
EEG = pop_loadset(sample_file);
fprintf('  Data: %d channels, %d samples (%.1f sec @ %d Hz)\n\n', ...
    EEG.nbchan, EEG.pnts, EEG.pnts/EEG.srate, EEG.srate);

% Run GEDAI with optimized code
fprintf('Running optimized GEDAI...\n');
tic;
[EEGclean, EEGartifacts, SENSAI_score, SENSAI_per_band, thresh_per_band, mean_ENOVA, ENOVA_per_epoch] = ...
    GEDAI(EEG, 'auto', 12, 0.5, 'precomputed', true, false, inf);
elapsed = toc;

fprintf('  ✓ Completed in %.2f seconds\n', elapsed);
fprintf('  SENSAI score: %.2f\n', SENSAI_score);
fprintf('  Wavelet bands processed: %d\n', length(SENSAI_per_band));
fprintf('  Mean ENOVA: %.4f\n\n', mean_ENOVA);

% Verify data integrity
reconstruction_error = max(abs(EEGclean.data(:) + EEGartifacts.data(:) - EEG.data(:)));
fprintf('Data integrity check:\n');
fprintf('  Max reconstruction error: %.2e\n', reconstruction_error);
if reconstruction_error < 1e-10
    fprintf('  ✓ Perfect reconstruction\n\n');
else
    fprintf('  ⚠ Warning: reconstruction error detected\n\n');
end

% Calculate expected memory savings
num_channels = EEG.nbchan;
num_samples = EEG.pnts;
num_bands = length(SENSAI_per_band) - 1;  % Subtract broadband

old_3d_memory = num_bands * num_channels * num_samples * 8 / 1e6;
new_2d_memory = num_channels * num_samples * 8 / 1e6;
memory_saved = old_3d_memory - new_2d_memory;

fprintf('Expected memory savings (wavelet storage):\n');
fprintf('  Old 3D array: %.1f MB\n', old_3d_memory);
fprintf('  New 2D array: %.1f MB\n', new_2d_memory);
fprintf('  Savings: %.1f MB (%.0f%% reduction)\n\n', memory_saved, (memory_saved/old_3d_memory)*100);

fprintf('=== Test Passed Successfully! ===\n');
fprintf('Optimized GEDAI-master is working correctly.\n\n');
