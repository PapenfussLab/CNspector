#!/bin/bash

usage="`basename $0` [-h] [-c] [-r] [-o] [-m] [-b] [-n] [-t] [-a] [-g] [-p] -- This script calls the cnvkit batch command for germline/tumour hyb capture analysis.
	where:
                -h show the help
		-c cnvkit path				[REQUIRED. Path to the cnvkit tool folder]
                -r reference				[REQUIRED. reference genome fasta file]
		-o output path				[DEFAULT: CWD]
		-m mode					[DEFAULT: Somatic. Somatic/Germline]
		-b tumour list				[REQUIRED. A list containing tumour samples bam full paths]
		-n normals list				[REQUIRED. A list containing normal samples bam full paths]
		-t target bed				[REQUIRED]
		-a anti target bed			[REQUIRED]
		-g cnvkit access file			[REQUIRED]
		-p cpus					[REQUIRED]"	

if [ "$1" == "-h" ]; then
        echo "$usage"
        exit 0
fi
while getopts ":h:c:r:o:m:b:n:t:a:g:p:" args; do
        case $args in
                h ) echo "$usage"
                    exit ;;
		c ) cnvkit_folder=$OPTARG
		    ;;
		r ) reference_file=$OPTARG
		    ;;
		o ) outdir=$OPTARG
                    ;;
		m ) mode=$OPTARG
                    ;;
		b ) tumour_list=$OPTARG
		    ;;
		n ) normals_list=$OPTARG
		    ;;
		t ) target_bed=$OPTARG
		    ;;
		a ) antitarget_bed=$OPTARG
		    ;;
		g ) access_file=$OPTARG
		    ;;
		p ) cpus=$OPTARG
		    ;;
		\? ) echo "$usage"
                    exit 1 ;;
        esac
done

outdir=${outdir:-$(pwd)}
mode=${mode:-Somatic}

parameter_check ()
{
	echo "ERROR: $1 not found"
	exit 1;
}

if [[ ! -d "$cnvkit_folder" ]] ; then
	parameter_check "cnvkit tool" 
elif [[ ! -f "$reference_file" ]] ; then
	parameter_check "reference fasta file"
elif [[ ! -f "$tumour_list" ]] ; then
	parameter_check "tumour list"
elif [[ ! -f "$normals_list" ]] ; then
	parameter_check "normals list"
elif [[ ! -f "$target_bed" ]] ; then
	parameter_check "target bed"
elif [[ ! -f "$antitarget_bed" ]] ; then
	parameter_check "antitarget bed"
elif [[ ! -f "$access_file" ]] ; then
        parameter_check "sequence-accessible file"
elif [[ ! -n "$cpus" ]] ; then
	parameter_check "#no of cpus parameter"
else	
	cnvkit_output="$outdir/cnvkit_output"
	cnvkit_reference="$outdir/cnvkit_output/reference"
	mkdir -p $cnvkit_reference
	
	echo "MESSAGE: output dir	- $cnvkit_output"
	echo "MESSAGE: reference genome - $reference_file"
	echo "MESSAGE: Mode 		- $mode"
	echo "MESSAGE: target bed	- $target_bed"
	echo "MESSAGE: anti-target bed	- $antitarget_bed"
	echo "MESSAGE: cnvkit access bed- $access_file"
	echo "MESSAGE: cnvkit ref outdir- $cnvkit_reference"

	declare control_bams
	while IFS=$"\n" read -r normals; do
		control_bams="$control_bams${normals[@]} "
	done < "$normals_list"
	control_bams="${control_bams% }"

	declare tumour_bams
	while IFS=$"\n" read -r tumours; do
		tumour_bams="$tumour_bams${tumours[@]} "
	done < "$tumour_list"
	tumour_bams="${tumour_bams% }"

	if [ "$mode" == "Somatic" ] ; then
		$cnvkit_folder/cnvkit.py batch -m hybrid --drop-low-coverage -p $cpus -n $control_bams -f $reference_file -t $target_bed -a antitarget_bed -g $access_file --output-reference $cnvkit_reference/reference.cnn -d $cnvkit_output --scatter --diagram $tumour_bams	

	elif [ "$mode" == "Germline" ] ; then
		$cnvkit_folder/cnvkit.py batch -m hybrid --drop-low-coverage -p $cpus -n $control_bams -f $reference_file -t $target_bed -a antitarget_bed -g $access_file --output-reference $cnvkit_reference/reference.cnn -d $cnvkit_output --scatter --diagram $tumour_bams
	fi
fi
