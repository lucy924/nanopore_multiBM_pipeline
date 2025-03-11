#!/bin/bash -e

# Define input and output variables
output_bam=${snakemake_output[0]}
input_dir=${snakemake_input[0]}
temp_dir=${snakemake_params[0]}
threads=${snakemake_params[1]}
log_file=${snakemake_log[0]}

# Temporary file to hold intermediate merge results
mkdir -p ${temp_dir}
temp_bam="${temp_dir}/temp_merged.bam"
temp_int_bam="${temp_dir}/temp_merged_new.bam"

# Get the list of BAM files
bam_files=($input_dir/*.bam)

# Check if there are any BAM files
if [ ${#bam_files[@]} -eq 0 ]; then
    echo "No BAM files found in $input_dir" >&2
    exit 1
fi

# Start merging process
echo "Merging BAM files sequentially..." > "$log_file"

# Initialize with the first file
cp "${bam_files[0]}" "$temp_bam"
# Note this was initially "mv" and I think it probably corrupted the first file and maybe the second

# Merge remaining BAM files sequentially
for bam in "${bam_files[@]:1}"; do
    samtools merge --threads "$threads" -f -o "$temp_int_bam" "$temp_bam" "$bam" 2>> "$log_file"
    mv "$temp_int_bam" "$temp_bam"
done

# Move final merged file to the output path
mv "$temp_bam" "$output_bam"
rm -r ${temp_dir}

echo "Merging completed successfully." >> "$log_file"
