#!/bin/bash

while getopts i:o:e:s: flag; do
  case "${flag}" in
  i) input_file=${OPTARG} ;;
  o) output_file=${OPTARG} ;;
  e) environment=${OPTARG} ;;
  s) strands=${OPTARG} ;;
  *) echo "usage: $0 [-i] [-o] [-e] [-s]" >&2
    exit 1;;
  esac
done

if [ -z "$environment" ]; then
  echo "No environment set..."
else
  # Check meme suit is installed
  ENVS=$(conda env list | awk -v awkvar="$environment" '{print awkvar}')
  if [[ $ENVS == *"$environment"* ]]; then
    echo "Activating environment..."
    source activate $environment
  else
    echo "Error: Please provide a valid virtual environment. For a list of valid virtual environment, please see 'conda env list' "
    exit
  fi
fi

if ! type "oligo-analysis" >/dev/null; then
  echo "RSAT is not installed. Please install use RSAT tools. Exiting program."
  exit
fi

# Check whether to use both strands
if [[ "${strands}" == "double" ]]; then
  strand_param="-2str"
else
  strand_param="-1str"
fi

# Generate background file
oligo-analysis -quick -v 1 "${strand_param}" -i "${input_file}" -format fasta -seqtype dna -noov -l 2 -type dna -return freq,occ -o "${output_file}"
