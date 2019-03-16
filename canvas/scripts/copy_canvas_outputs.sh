#!/bin/sh
#tool_output - path to the canvas output folder 
tool_output=$1
dest=canvas_output
mkdir canvas_output
for f in `ls $tool_output/*/TempCNV_*/*raw* $tool_output/*/TempCNV_*/*raw*  $tool_output/*/TempCNV_*/*clean* $tool_output/*/TempCNV_*/*normal* `
do
	echo Copying $f
	cp -v $f $dest/`basename $f`.gz
done
for f in `ls -d $tool_output/* | grep -v vcf`;
do
	dest_file=`echo $f | sed 's/^.*\///g'`
	echo $dest_file
	cp -v  $f/CNV.vcf.gz $dest/$dest_file.vcf.gz
done
