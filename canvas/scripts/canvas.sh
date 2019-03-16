#!/bin/bash

usage="`basename $0` [-h] [-a] [-c] [-r] [-m] [-o] [-s] [-v] [-t] [-b] [-n] -- This script will run canvas tool for somatic/germline samples. Tools set up to be done prior to running the script.
	where:
                -h show the help
		-a mono path				[REQUIRED. Path to the folder containing mono executable]
		-c canvas path				[REQUIRED. Path to the canvas tool folder]
                -r reference folder path		[REQUIRED. Need to contains genome.fa, kmer.fa, filter13.bed and corresponding .fai files]
		-m manifest file			[REQUIRED]
		-o output path				[DEFAULT: CWD/output/samplename]
		-s sample name				[REQUIRED. Needed for canvas work flow]
		-v sample vcf file			[REQUIRED]
		-t mode					[DEFAULT: Somatic. Somatic/Germline]
		-b sample bam				[REQUIRED]
		-n normals list				[REQUIRED. A list containing normal samples bam full paths]"
	
if [ "$1" == "-h" ]; then
        echo "$usage"
        exit 0
fi
while getopts ":h:a:c:r:m:o:s:v:t:b:n:" args; do
        case $args in
                h ) echo "$usage"
                    exit ;;
		a ) mono_folder=$OPTARG
		    ;;
		c ) canvas_folder=$OPTARG
                    ;;
		r ) reference_folder=$OPTARG
		    ;;
		m ) manifest=$OPTARG
		    ;;
		o ) outdir=$OPTARG
                    ;;
		s ) sampleoutname=$OPTARG
		    ;;
		v ) samplevcf=$OPTARG
		    ;;
		t ) mode=$OPTARG
                    ;;
		b ) test_bam=$OPTARG
		    ;;
		n ) normals_list=$OPTARG
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

if [[ ! -d "$mono_folder" || ! -d "$canvas_folder" ]] ; then
	parameter_check "mono and/or canvas" 
elif [[ ! -d "$reference_folder" ]] ; then
	parameter_check "reference folder"
elif [[ ! -f "$manifest" ]] ; then
	parameter_check "manifest file"
elif [[ ! -n "$sampleoutname" ]] ; then
	parameter_check "sample name"
elif [[ ! -f "$samplevcf" ]] ; then
	parameter_check "sample vcf"
elif [[ ! -f "$test_bam" ]] ; then
	parameter_check "sample bam file"
elif [[ ! -f "$normals_list" ]] ; then
	parameter_check "normal samples list"
else	
	canvas_output="$outdir/output/$sampleoutname"
	filter_bed="$reference_folder/filter13.bed"
	kmer_file="$reference_folder/kmer.fa"
	mkdir -p $canvas_output
	
	echo "MESSAGE: output dir	- $canvas_output"
	echo "MESSAGE: reference genome - $reference_folder/genome.fa"
	echo "MESSAGE: Mode 		- $mode"
	echo "MESSAGE: filter bed file  - $filter_bed"
	echo "MESSAGE: kmer file	- $kmer_file"	

	declare control_bams
	while IFS=$"\n" read -r normals; do
		control="--control-bam=${normals[@]} "
		control_bams="$control_bams$control"
	done < "$normals_list"
	control_bams="${control_bams% }"

	if [ "$mode" == "Somatic" ] ; then
		$mono_folder/mono $canvas_folder/Canvas.exe Somatic-Enrichment -b $test_bam $control_bams --manifest=${manifest} --b-allele-vcf=${samplevcf} -n ${sampleoutname} -o $canvas_output -r $kmer_file -g $reference_folder -f $filter_bed
	
	elif [ "$mode" == "Germline" ] ; then
		$mono_folder/mono $canvas_folder/Canvas.exe Somatic-Enrichment --custom-parameters=CanvasSomaticCaller,--definedpurity=1 -b $test_bam $control_bams --manifest=${manifest} --b-allele-vcf=${samplevcf} -n ${sampleoutname} -o $canvas_output -r $kmer_file -g $reference_folder -f $filter_bed
	fi
fi
