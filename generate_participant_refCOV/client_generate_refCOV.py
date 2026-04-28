import os
import subprocess
import pandas as pd
import argparse
from eeg_file_operations import raw_reststate_data_finder  # custom function for file loading. Can be modified to take an input path for general use.
import papermill as pm

'''
Compute participant reference covariance structure from leadfield for individualized GEDAI denoising
'''
# Modify input paths here
subject_dir = '/home/imk2003/Documents/freesurfer_outputs/'
sources_savedir = '/athena/grosenicklab/scratch/imk2003/eeg_sources_data/'

# Parse ppt_id, tms_target, and analysis_day variables provided by submission script
# The passed variables (tms_target, diagnosis, etc.) can be modified to fit your individual study.
# Using individual head circumference measurements (if available) can be helpful for improving forward model fit.
parser = argparse.ArgumentParser()
parser.add_argument('--subject_list', required=True)
parser.add_argument('--ppt_id', required=True)
parser.add_argument('--tms_target', required=True)
parser.add_argument('--diagnosis', required=True)
parser.add_argument('--head_circumference', required=True)
args = parser.parse_args()
subject_list_csv = args.subject_list
ppt_id = args.ppt_id
tms_target = args.tms_target
diagnosis = args.diagnosis
head_circumference = args.head_circumference

# Read subject counts csv and extract MRI id and tms target info
subjects_df = pd.read_csv(subject_list_csv)
subjects_df.set_index('record_id', inplace=True)  # Set record_id as index for filtering
mri_id = subjects_df.loc[ppt_id, 'mri_id']

# Passes paths to ppt's pre or post-tx EEG for analysis day, tms_target, mri_id, and Freesurfer subjects_dir
# through environment variables to compute_inverse_solution.py and runs script
def run_refCOV(ppt_id, eeg_filepath):
    plotting_notebook = 'generate_refCOV.ipynb'
    filename = os.path.basename(eeg_filepath)
    filename = os.path.basename(eeg_filepath)
    if 'crossover' in ppt_id:
        parts = filename.split('_')
        parts.insert(1, 'crossover')
        filename = '_'.join(parts[:6])
    else:
        filename = '_'.join(filename.split('_')[:5])

    # Run forward solution computation on raw EEG for the purpose of generating refCOV for GEDAI
    # and plot covariance figs as check
    try:
        print(f'Running refCOV generation for {filename}')
        
        # Create directory for saving papermill notebook metadata
        os.makedirs('refCOV_plots', exist_ok=True)
        
        # Execute notebook via papermill for interactive plotting with the 'notebook' backend
        pm.execute_notebook(
            plotting_notebook,
            f'refCOV_plots/{ppt_id}_{tms_target}_generate_refCOV.ipynb',  # Save executed version
            parameters=dict(
                ppt_id=ppt_id,
                subject=mri_id,
                head_circumference=float(head_circumference),
                tms_target=tms_target,
                filename=filename,
                
                # Read paths to EEG data, source estimates dir, and Freesurfer subjects dir
                savedir=sources_savedir,
                subjects_dir=subject_dir,
                
                # Read paths to stc, raw preprocessed EEG, preprocessed phantom EEG
                raw_path=eeg_filepath)
        )
        print(f'Successfully generated refCOV for {filename}')
        
    except FileNotFoundError as e:
        print(f'File not found for {filename}: {e}')
        print(f'Check that {plotting_notebook} exists and paths are correct')
        
    except Exception as e:
        print(f'Error executing notebook for {filename}')
        print(f'Error type: {type(e).__name__}')
        print(f'Error message: {e}')
        print(f'Check executed_notebooks/{ppt_id}_{plotting_notebook} for details')


# Obtain a list of all raw resting state EEG files available for this participant
ppt_raw_filepaths = raw_reststate_data_finder(ppt_id, tms_target, diagnosis)

# Leadfield and therefore reference covariance structure are the same for a given participant over their recordings
# So just select one participant recording and create leadfield + refCOV for ppt based on that
ppt_raw_filepath = ppt_raw_filepaths[0]
print(f'Initializing refCOV estimation for {ppt_raw_filepath}.')

# Run forward model estimation 
try:
    run_refCOV(ppt_id, ppt_raw_filepath)
except Exception as e:
    print(e)

