							CNSpector
                                                           
CNSpector is a web-based browser to visualize Copy Number Variation (CNV) calls. The browser is an R program for shiny-server that reads an input file containing a table whose rows point to the files to display for a given sample or group of samples (multi-sample mode). Currently it can support output from cnvkit, copywriter and canvas tools. The outputs from these tools are converted in a format accepted by CNSpector.

Installation

Git clone the repository as below or download as required. The folder structure needs to remain the same able to run the scripts or load the references correctly. 
git clone https://github.com/PapenfussLab/CNspector.git 

shiny-server scripts

The following R files are provided for CNSpector browser. It was tested with R 3.3.1, 3.4.0 and 3.5.1 versions and can be launched through Rstudio.

		Helpers.R 
		Server.R
		ui.R

Reference files

The following reference files are provided to support tools output conversion and for the browser in terms of annotations and bin sizes for computations – 
 
 1. annotations/browser
  
  Required for the visualization in the browser.
  
	  chromInfo.txt
	  cytoBand.txt
	  panel_transcripts.bed 
	  targeted_transcripts.bed
	  wg_transcripts.bed

  2. annotations/preprocessing
  
  Required for converting different tools output.
  
	  hg19_targeted.bed
	  hg19_wg_coarse.bed
	  hg19_wg_fine.bed
	  hg19_wg_medium.bed

  <hg_19_targeted.bed/targeted_transcripts.bed/panel_transcripts.bed> - Bed files to be generated as per the panel design.

R package dependencies

  Browser
  
  The browser can be launched through Rstudio and requires the following packages to run the application correctly – 
  
	  library(futile.logger)
	  library(plyr)
	  library(data.table)
	  library(R.utils)
	  library(RColorBrewer)
	  library(ggplot2)
	  library(Cairo)

  Conversion scripts 

	  library(genomation)
	  library(GenomicRanges)
	  library(stats)
	  library(VariantAnnotation)


Convert different CNV tools output – 

Refer to the basic instructions provided in the documentation/CNVTools_README to help run other cnv tools (not limited to this but can be run directly by following the actual tool documentation). Conversion script steps are list here again which converts the different tool outputs to browser accepted input file.

1.	Cnvkit

		a.	The cnvkit_output folder needs to be placed inside the cnvkit folder. 
		b.	Run the following commands to generate the “index.tsv” file for browser
			cd cnvkit
			Rscript scripts/import.R 

2.	Canvas

		a.	The “canvas_output” folder needs to be location inside canvas folder.
		b.	Run the following commands to generate the CNSpector/index.tsv file for browser
			cd canvas
			Rscript scripts/import.R

3.	CopyWriter

		a.	copywriter_output folder needs to be places inside the copywriter folder. 
		b.	Run the following commands to generate CNSpector/index.tsv file for browser.
			cd copywriter
			Rscript scripts/import.R

4.	Multiple tool combination 
	If successful running all the 3 tools, then best resolution from each tool is put together to view in the browser as follows – 

		scripts/make_combination_index.sh
		Output: combination/index.tsv


Rstudio run Instructions 

Open the server.R through Rstudio and run the application to launch the shiny application as follows – 

		http://shiny.server.url/cnb/?file=<path_to_the_file.tsv>

Supported options in the URL link - 

		reference_path=REFERENCE_PATH
		base_path=BASE_PATH
		run=RUN_NAME
		sample=SAMPLE_NAME
		gene_id=GENE_NAME
		locus=CHR:START-END
		debug=TRUE|FALSE
		verbose=TRUE|FALSE
		locked=TRUE|FALSE

Example to load the “index.tsv” file into the browser (specific tool / combination index )

		http://shiny.server.url/s/8c75ccfc78a3183671b0d/p/7743/?file=prefix_path/ CNspector/canvas/CNSpector/index.tsv&reference_path=prefix_path/CNSpector/annotations/browser&base_path=/home/canvas/CNSpector

Debug mode – 
	
		http://shiny.server.url/s/8c75ccfc78a3183671b0d/p/7743/?file=prefix_path/ CNspector/canvas/CNSpector/index.tsv&reference_path=prefix_path/CNSpector/annotations/browser&base_path=/home/canvas/CNSpector&debug=T

