import os
import subprocess
import pandas as pd
import argparse
import logging
import papermill as pm

'''
Compute participant reference covariance structure from leadfield for individualized GEDAI denoising
'''
# Modify input paths here
subjects_dir = '/home/imk2003/Documents/freesurfer_outputs/'  # path to Freesurfer subjects directory
sources_savedir = '/athena/grosenicklab/scratch/imk2003/eeg_sources_data/'  # data save dir
store_eeg_datapath = '/athena/grosenicklab/store/tms_eeg/'

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

# Given a data path, generates a list of all raw resting state EEG files for an identified participant
# And returns the first available file path. The leadfield and therefore reference covariance structure
# should be the same for a given participant over time.
def raw_reststate_data_finder(ppt_id, tms_target, diagnosis, eeg_datapath=store_eeg_datapath):
    reststate_filepaths = []
    tms_target = tms_target.lower()
    diagnosis = diagnosis.lower()
    data_folder = f'{diagnosis}_{tms_target}'
    datapath = os.path.join(eeg_datapath, data_folder)

    # Set lookup ID. Avoids errors with matching when 'crossover' is in ppt_id.
    lookup_id = ppt_id.split('_')[0] if 'crossover' in ppt_id else ppt_id

    try:
        ppt_path = [entry.path for entry in os.scandir(datapath) if lookup_id in entry.name][0]
        day_paths = [entry.path for entry in os.scandir(ppt_path) if entry.name.startswith(lookup_id)]

        for path in day_paths:
            reststate_filepaths.extend(
                entry.path for entry in os.scandir(path) if 'reststate' in entry.name
            )

    except IndexError:  # No matching participant folder in store_eeg_datapath
        logging.warning(f'No EEG data found for participant: {ppt_id}. Check participant ID and file paths.')
    except OSError as e:
        logging.error(f'File system error while searching for participant {ppt_id}: {e}')

    if not reststate_filepaths:  # Successfully indexed into ppt_path but found no matching files
        logging.warning(f'No files found for {ppt_id}.')
        return None

    return reststate_filepaths[0]


# Passes paths to ppt's pre or post-tx EEG for analysis day, tms_target, mri_id, and Freesurfer subjects_dir
# through environment variables to compute_inverse_solution.py and runs script
def run_refCOV(ppt_id, eeg_filepath):
    plotting_notebook = 'generate_refCOV.ipynb'
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
                subjects_dir=subjects_dir,
                
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


# Identify the path to the raw EEG file to be denoised with GEDAI
ppt_raw_filepath = raw_reststate_data_finder(ppt_id, tms_target, diagnosis)
print(f'Initializing refCOV estimation for {ppt_raw_filepath}.')

# Run forward model estimation 
try:
    run_refCOV(ppt_id, ppt_raw_filepath)
except Exception as e:
    print(e)

