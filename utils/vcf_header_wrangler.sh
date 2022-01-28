#!/bin/bash

while getopts i:o:r: flag; do
  case "${flag}" in
  i) input_folder=${OPTARG} ;;
  o) output_folder=${OPTARG} ;;
  r) reference=${OPTARG} ;;
  *) echo "usage: $0 [-i] [-o] [-r]" >&2
    exit 1;;
  esac
done

function ProgressBar() {
  # Process data
  ((_progress = (${1} * 100 / ${2} * 100) / 100))
  ((_done = (${_progress} * 4) / 10))
  ((_left = 40 - $_done))

  # Build progressbar string lengths
  _done=$(printf "%${_done}s")
  _left=$(printf "%${_left}s")

  # Display and update the progress bar
  printf "\rProgress : [${_done// /#}${_left// /-}] ${_progress}%%"
}

total_operations=$(find $input_folder/*.vcf -maxdepth 0 | wc -l)
printf "Input Folder:\t $input_folder\n"
printf "Output Folder:\t $output_folder\n"
printf "Total files found:\t $total_operations\n"

file_counter=1
start=$(date +%s)
for FILE in $input_folder/*; do
  base_file_name="$(basename -- $FILE)"
  converted_file_name="$output_folder/$base_file_name"
  awk '!/^##/' $FILE >$converted_file_name
  cat $reference $converted_file_name >$converted_file_name
  ProgressBar ${file_counter} ${total_operations}
  ((file_counter++))
done

end=$(date +%s)
runtime=$((end - start))
printf "\nTime taken:\t $runtime seconds."
printf "\nDone.\n"
