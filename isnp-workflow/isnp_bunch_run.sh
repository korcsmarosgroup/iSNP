#!/bin/bash

POSITIONAL=()
while [[ $# -gt 0 ]]; do
	key="$1"
	case $key in
	-t|--thread)
		ThreadNum="$2"
		shift
		shift
		;;
	--no-docker-build)
		NoDockerBuild=yes
		shift
		;;
	--vcf-files)
		VcfFiles="$2"
		shift
		shift
		;;
	-c|--clear-output)
		ClearOutput=yes
		shift
		;;
	-s|--separate-tasks-in-containers)
		Separate=--separate
		shift
		;;
	-h|--help)
		echo "Usage:"
		echo " "
		echo " isnp_bunch_run.sh [OPTIONS] [<relative path to directory contains vcf files, without input-data dir>]"
		echo " "
		echo "options:"
		echo "-t|--thread <number of threads>      : specify the number of 'isnp.py' scripts will be executed parallel "
		echo " "
		echo "--no-docker-build                    : will not build the docker images (expects you already built them)  "
		echo "                                       If you don't provide this option, the docker images will be built"
		echo "                                       once before the starting of parallel execution."
		echo ""
		echo "--vcf-files                          : you can manually give a set of paths to VCF files separated by ','"
		echo "                                       which will be executed by the script (the default behavior is to "
		echo "                                       provide a folder where all the VCF files will be taken)"
		echo ""
		echo "-s|--separate-tasks-in-containers    : using this option, all the analytical tasks will be executed in "
		echo "                                       separate docker containers (default behavior: execute single "
		echo "                                       container for each patient and run all the steps for a patient in "
		echo "                                       it's container)"
		echo " "
		exit 0
		;;
    *)
    	POSITIONAL+=("$1")
    	shift
    	;;
    esac
done
set -- "${POSITIONAL[@]}"

# post process parameters
NoDockerBuild=${NoDockerBuild:-no}
ClearOutput=${ClearOutput:-no}
ThreadNum=${ThreadNum:-1}
IFS=',' read -r -a VcfFileList <<< "$VcfFiles"
VcfNum=${#VcfFileList[@]}
IterDir=$1

if [ $VcfNum == 0 ] && [ -z "$IterDir" ]
then
	echo error
	exit 1
fi

function run_isnp {
	# sleep random number of secs
	#sleepnum=$((20 + RANDOM % 100))
	#echo "Sleep for $sleepnum"
	#sleep $sleepnum

	#vcffile=norwoch-cohort-cleaned/$(basename $1)
	vcffile=$2/$(basename $1)
	outdir=$(basename $vcffile)

	# remove outfiles
	#if [ $ClearOutput == 'yes' ]
	#then
	#	rm -rf output-data/$outdir.log
	#	rm -rf output-data/$outdir.err
	#	rm -rf output-data/$outdir
	#fi

	optional_separate_param=$3
	
	# run isnp
	start=`date +%s`
	python3 isnp.py \
		-i `pwd`/input-data \
		-o `pwd`/output-data/$outdir \
		--patient_vcf $vcffile \
		--debug \
		--no-docker-build \
		$optional_separate_param \
		> output-data/$outdir.log \
		2> output-data/$outdir.err
	wait

	# calculate and print run time
	end=`date +%s`
	runtime=$((end-start))
	echo "Done:  $runtime"
}
export -f run_isnp

# build docker containers
if [ $NoDockerBuild != 'yes' ]
then
	python3 isnp.py --only-build-docker --debug
	wait
fi

if [ $VcfNum -gt 0 ]
then
	printf '%s\n' "${VcfFileList[@]}" | xargs -n 1 -P $ThreadNum -I % bash -c 'run_isnp "$@"' _ % $IterDir $Separate
else
	# run isnp parallel
	#ls input-data/norwoch-cohort-cleaned/*.vcf | parallel -j 8 run_isnp
	#ls input-data/norwoch-cohort-cleaned/*.vcf | xargs -n 1 -P 6 -I % bash -c 'run_isnp "$@"' _ %
	ls input-data/$IterDir/*.vcf | xargs -n 1 -P $ThreadNum -I % bash -c 'run_isnp "$@"' _ % $IterDir $Separate
fi
