#!/bin/sh

mkdir combination
head -n 1 copywriter/CNspector/index.tsv > combination/index.tsv
cat copywriter/CNspector/index.tsv | awk '{if($2=="reference") print "annotations/"$0}' >> combination/index.tsv
cat copywriter/CNspector/index.tsv | awk '{if($2=="counts" && $4 != "targeted") print "copywriter/CNspector/"$0}' >> combination/index.tsv
cat canvas/CNspector/index.tsv  |  awk '{if($2=="counts" && $4 == "targeted") print "canvas/CNspector/"$0}'   >> combination/index.tsv
cat cnvkit/CNspector/index.tsv  | awk '{if($2 == "baf" || $2=="deletions") print "cnvkit/CNspector/"$0}' >> combination/index.tsv
