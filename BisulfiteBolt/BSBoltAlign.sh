# adjust the maximum number of jobs to run in parallel based on computing resources
max_jobs=3

echo "Found $(ls *.fq.gz 2>/dev/null | wc -l) files to process"
echo "Using $max_jobs parallel jobs"
echo "Starting processing at $(date)"
echo "========================================"

# Function to wait for jobs when we hit the limit
wait_for_jobs() {
    while [ $(jobs -r | wc -l) -ge $max_jobs ]; do
        sleep 1
    done
}

# Process each file. Assuming paired end reads. If using single end reads, remove lines 23-45 and line 54 (-F2) from bsbolt Align
for sample_id in *_1_val_1.fq.gz; do
    wait_for_jobs

    {
        base_id=$(basename "$sample_id" "_1_val_1.fq.gz")

        # Properly construct R2 filename
        # If your files are like: FD_22_WT_004_R1_val_1.fq.gz
        # Then R2 should be:     FD_22_WT_004_R2_val_2.fq.gz
        r2_file="${base_id}_2_val_2.fq.gz"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting: $base_id"
        echo "R1 file: $sample_id"
        echo "R2 file: $r2_file"


        # Check if both files exist before proceeding
        if [[ ! -f "$sample_id" ]]; then
            echo "ERROR: R1 file missing: $sample_id"
            exit 1
        fi

        if [[ ! -f "$r2_file" ]]; then
            echo "ERROR: R2 file missing: $r2_file"
            echo "Checked these locations:"
            echo "  ${base_id}_2_val_2.fq.gz"
            echo "  ${base_id}_val_2.fq.gz"
            exit 1
        fi

        echo "Output: ${base_id}_output"

        start_time=$(date +%s)

        # Use full paths and existing files. Adjust -t based on the number of threads to use for alignment
        bsbolt Align -t 8 \
            -F1 "/PATH/$sample_id" \
            -F2 "/PATH/$r2_file" \
            -DB /PATH/BSBolt \
            -O /PATH/BSBolt/"${base_id}_output"

        end_time=$(date +%s)
        duration=$((end_time - start_time))

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed: $base_id (took ${duration}s)"
        echo "----------------------------------------"
    } &
done

echo "All jobs started, waiting for completion..."
wait

echo "========================================"
echo "All processing completed at $(date)"
