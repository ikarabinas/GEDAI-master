#!/bin/bash
# Batch submission of forward modeling computation to estimate individual participant reference covariance (refCOV).
# Takes a csv with participant IDs, file identifiers, and corresponding forward modeling parameters as input.

# Load an environment with mne-python installed
source ~/.bashrc
mamba activate mne_offscreen

# Specify paths to: 
# 1. A csv listing ppt IDs and corresponding info needed for file identification and forward modeling (e.g. group tms_target, head_circumference)
# 2. The path to the client script for running refCOV computation ('client_generate_refCOV.py')
# 3. The path to a save folder
# 4. The directory for log saving
# 5. The path to a job tracking file ('refCOV_job_tracking.tsv')

SUBJECT_LIST="/home/imk2003/Documents/updated_subject_list_rofc.csv"
COMPUTE_REFCOV_SCRIPT="/home/imk2003/Documents/MATLAB/eeglab/plugins/GEDAI-master/generate_custom_refCOV/client_generate_refCOV.py"
SAVEPATH="/athena/grosenicklab/scratch/imk2003/eeg_sources_data/"
LOG_DIR="/athena/grosenicklab/scratch/imk2003/eeg_sources_data/generate_refCOV_logs/"
TRACKING_FILE="refCOV_job_tracking.tsv"

# Set to true to simulate submission without running sbatch
DRY_RUN=false

# Create log dir and write tracking file
mkdir -p "$LOG_DIR"
echo -e "subject_id\tjob_id\ttimestamp" > "$TRACKING_FILE"

# Read subject ID, group tms_target, and head circumference info from csv
# Match the csvcut command numbers to column numbers for specified variables.
csvcut -c 2,5,7,8 "$SUBJECT_LIST" | tail -n +2 | while IFS=',' read -r subject_id tms_target diagnosis head_circumference; do
    subject_id=$(echo "$subject_id" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    tms_target=$(echo "$tms_target" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    diagnosis=$(echo "$diagnosis" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    head_circumference=$(echo "$head_circumference" | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    
    # Skip entries with blank record_id, tms_target, or head circumference
    if [ -z "$subject_id" ] || [ -z "$tms_target" ] || [ -z "$head_circumference" ] || [ -z "$diagnosis" ]; then
        echo "Skipping ${subject_id}. Missing key info for this participant in the provided csv."
        continue
    fi

    # Reformat tms target to lowercase
    tms_target="${tms_target,,}"  

    # Skip subjects with an existing refCOV file in savepath
    if compgen -G "${SAVEPATH}/GEDAI_refCOV/${subject_id}*refCOV.mat" > /dev/null; then
        echo "Skipping ${subject_id}_${tms_target}. refCOV file already exists."
        continue
    fi


    # Build SLURM command. Allocate 16-20G for creating refCOV
    SLURM_CMD="sbatch --mem=16G --cpus-per-task=2 \
        --job-name=stc_${subject_id} \
        --partition=sackler-cpu,scu-cpu \
        --time=3:00:00 \
        --output=${LOG_DIR}/${subject_id}_${tms_target}_%j.out \
        --error=${LOG_DIR}/${subject_id}_${tms_target}_%j.err \
        --wrap=\"python -u ${COMPUTE_REFCOV_SCRIPT} --subject_list ${SUBJECT_LIST} --ppt_id ${subject_id} --tms_target ${tms_target} --diagnosis ${diagnosis} --head_circumference ${head_circumference}\""

    # Submit job or test with dry run
    if [ "$DRY_RUN" = true ]; then
        echo "[DRY RUN] Would submit: $SLURM_CMD"
    else
        job_output=$(eval "$SLURM_CMD")

        # Job submission error handling
        if [[ "$job_output" != Submitted* ]]; then
            echo "ERROR submitting ${subject_id}_${tms_target}. sbatch said: $job_output"
            continue
        fi

        # Job tracking
        job_id=$(echo "$job_output" | awk '{print $4}')
        timestamp=$(date +"%Y-%m-%d %H:%M:%S")
        echo -e "${subject_id}_${tms_target}\t${job_id}\t${timestamp}" >> "$TRACKING_FILE"
        echo "Submitted ${subject_id}_${tms_target} as Job ${job_id}"
    fi

done
