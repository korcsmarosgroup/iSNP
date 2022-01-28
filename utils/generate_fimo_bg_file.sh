#!/bin/bash

while getopts i:o:e:m: flag
do
  case "${flag}" in
    i) input_file=${OPTARG};;
    o) output_file=${OPTARG};;
    e) environment=${OPTARG};;
    m) order=${OPTARG};;
    *) echo "usage: $0 [-i] [-o] [-e] [-m]" >&2
         exit 1 ;;
  esac
done

if [ -z "$environment" ];
then
  echo "No environment set..."
else
  # Check meme suit is installed
  ENVS=$(conda env list | awk -v awkvar="$environment" '{print awkvar}' )
  if [[ $ENVS = *"$environment"* ]]; then
     echo "Activating environment..."
     source activate $environment
  else
     echo "Error: Please provide a valid virtual environment. For a list of valid virtual environment, please see 'conda env list' "
     exit
  fi;
fi

if ! type "fasta-get-markov" > /dev/null; then
  echo "MEME Suite is not installed. Please install use FIMO tools. Exiting program."
  exit
fi

# Generate background file
fasta-get-markov -m "${order}" "${input_file}" "${output_file}"
