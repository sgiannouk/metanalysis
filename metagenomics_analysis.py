 ###Stavros Giannoukakos###
##University of Granada##

#Version of the program
__version__ = "0.2.0"

import argparse
import operator
from Bio import SeqIO
from Bio.Seq import Seq
import subprocess, gzip
from datetime import datetime
import shutil, fnmatch, glob, sys, os



batch_to_analyse = ["batch1", "batch2"]
inputDir = '/shared/projects/martyna_rrna_illumina_igtp/blood_microbiome'
metadata_file = "/shared/projects/martyna_rrna_illumina_igtp/blood_microbiome/clinical_data.txt"

# batch_to_analyse = ["batch2",  "batch3", "batch4", "batch5"]
# inputDir = '/shared/projects/martyna_rrna_illumina_igtp/balf_microbiome'
# metadata_file = "/shared/projects/martyna_rrna_illumina_igtp/balf_microbiome/clinical_data.txt"

rscripts = f"{os.path.dirname(os.path.realpath(__file__))}/Ranalysis"
# Calling Qiime2 tools from its conda environment
qiime2_env = "/home/stavros/playground/progs/anaconda3/envs/qiime2-2021.2/bin"
# Configuration file needed for FastQ Screen
fastQscreen_config = "/home/stavros/references/fastQscreen_references/fastq_screen.conf"
# SILVA database
silva_taxonomy = "/home/stavros/references/metagenomics_dbs/SILVA/silva-138-99-tax.qza"
silva_reference = "/home/stavros/references/metagenomics_dbs/SILVA/silva-138-99-seqs.qza"
silva_reference_targeted = "/home/stavros/references/metagenomics_dbs/SILVA/silva-138-99-seqs_targeted.qza"
# silva_99_classifier = "/home/stavros/references/metagenomics_dbs/SILVA/silva_138_99_v3v4_scikitv0.23.1.qza"  # scikit-learn version 0.23.1.
silva_99_classifier = "/home/stavros/references/metagenomics_dbs/SILVA/silva_138_99_v3v4_scikitv0.24.2.qza"  # scikit-learn version 0.24.2



usage = "metagenomics_analysis [options] -i <input_directory/input_files>"
epilog = " -- January 2019 | Stavros Giannoukakos -- "
description = "DESCRIPTION"

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter, usage=usage, description=description, epilog=epilog)
# Create required section in help
requiredArgs = parser.add_argument_group('required arguments')
# Input folder option
requiredArgs.add_argument('-i', '--input_dir', dest='input_dir', metavar='', required=False,
						   help="Path of the input directory that contains the raw data.\nBoth forward and reverse reads are expected to be found\nin this directory.")
# Number of threads/CPUs to be used
parser.add_argument('-th', '--threads', dest='threads', default=str(50), metavar='', 
                	help="Number of threads to be used in the analysis")
# Number of threads/CPUs to be used
parser.add_argument('-fp', '--forwardPrimer', dest='forwardPrimer', default="CCTACGGGNGGCWGCAG", metavar='', required=False, 
                	help="Sequence of the forward primer")
# Number of threads/CPUs to be used
parser.add_argument('-rp', '--reversePrimer', dest='reversePrimer', default="GACTACHVGGGTATCTAATCC", metavar='', required=False,
                	help="Sequence of the reverse primer")
# Metadata file
parser.add_argument('-m', '--metadata', dest='metadata', required=False, metavar='', 
                	help="Metadata file containing several info concerning the\ninput data")
# Path of where the output folder should be located
parser.add_argument('-o', '--output_dir', dest='output_dir', default=os.path.dirname(os.path.realpath(__file__)), 
					help="Path of the directory that will host the analysis.\n(default: <current_directory>)")
# Display the version of the pipeline 
parser.add_argument('-v', '--version', action='version', version=f'%(prog)s {__version__}')

# Get the options and return them
args = parser.parse_args()
start_time = datetime.now()  # Tracking time of analysis

""" All the necessary directory that will host the analysis are being created. These 
include 'preprocessed_files' that will host the filtered and quality controlled 
data and reports, where all reports from all software will be stored. """
args.metadata = metadata_file
args.input_dir = inputDir

# Main directories hosting the analysis
analysis_dir = os.path.join(args.output_dir, "blood_metanalysis")
# analysis_dir = os.path.join(args.output_dir, "balf_metanalysis")
reports_dir = os.path.join(analysis_dir, "reports")  # Reports directory
prepr_dir = os.path.join(analysis_dir, "preprocessed_data")  # Save processed .fastq files
inputqiime_dir = os.path.join(analysis_dir, "qiime2_input_data")  # Input data directory
denoise_dir = os.path.join(analysis_dir, "denoising_analysis")  # Denoising output
taxonomy_dir = os.path.join(analysis_dir, "taxonomic_analysis")  # Taxonomic analysis
krona_dir = os.path.join(analysis_dir, "taxonomic_analysis/krona_visualisation")  # Krona taxonomy visual
diversity_dir = os.path.join(analysis_dir, "diversity_analysis")  # Diversity analysis
random_forest_dir = os.path.join(analysis_dir, "random_forest")
downstream_dir = os.path.join(analysis_dir, "downstream_analysis")
# Reporting directories
qc_reports = os.path.join(analysis_dir, "reports/QC_reports")
pipeline_reports = os.path.join(analysis_dir, "reports/pipeline_reports")
if not os.path.exists(pipeline_reports): os.makedirs(pipeline_reports)



def assess_inputdata(batch):
	""" In this function the PE input data will be assesses for valid format and 
	for correct pairing. """
	input_files = []  # Output list that will contain the paired-input files
	for path, subdirs, files in os.walk(f'{args.input_dir}/{batch}'):
		for name in files:
			# Verifying that all reads have their pairs
			if name.endswith((".fastq.gz", ".fq.gz")) and not name.startswith("rep") and not any(x in name.upper() for x in ["_R1_", "_R2_"]):
				sys.exit(f'\nWARNING - No paired read found for read {name} (in {batch})')  

			elif name.endswith((".fastq.gz", ".fq.gz")) and not name.startswith("rep") and "_R1_" in name:  # Obtaining the paired-input files
				inR1 = os.path.join(os.path.abspath(path), name)
				inR2 = inR1.replace("_R1_","_R2_")
				assert (os.path.isfile(inR2)), f'Could not detect the pair of {inR1} ({inR2})'
				input_files.append((inR1, inR2))
	
	input_files = sorted(input_files, key=lambda x: x[0])
	flatlist = [item for sublist in input_files for item in sublist if "_R1_" in item]
	return input_files, flatlist

def quality_control(batch, R1list):
	""" Running FastQ  Screen software to identify possible  contaminations in our samples. 
	Additionally, use AfterQC to make a preliminary quality check of the processed PE reads. 
	Then MultiQC  will summarise the QC reports from  all samples into a summary report """
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} QUALITY CONTROL OF THE PREPROCESSED SAMPLES ({batch})')
	


	prepr_qc_reports = f'{qc_reports}/{batch}'
	mfiltered_data = ' '.join(R1list)  # Obtaining the preprocessed reads
	if not os.path.exists(prepr_qc_reports): os.makedirs(prepr_qc_reports)
	

	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  1|4 FastQScreen - Checking random reads for possible contamination: in progress ..')
	fastQscreen = ' '.join([
	"fastq_screen",  # Call fastQ screen to check contamination in the processed data
	"--threads", args.threads,  # Number of threads to use
	"--outdir",  prepr_qc_reports,  # Directory in which the output files will be saved
	"--quiet",  # Suppress all progress reports on stderr and only report errors
	"--conf", fastQscreen_config,  # Location of the required configuration file
	mfiltered_data, mfiltered_data.replace("_R1_", "_R2_"),  # Input PE files
	"2>>", os.path.join(pipeline_reports, "preprocessing1_fastQscreen-report.log")])  # Output fastQ screen report
	subprocess.run(fastQscreen, shell=True)

	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  2|4 FastQC - Quality Control reports for the forward reads are being generated: in progress ..')
	fastQC_frw = ' '.join([
	"fastqc",  # Call fastQC to quality control all processed data
	"--threads", args.threads,  # Number of threads to use
	"--quiet",  # Print only log warnings
	"--outdir", prepr_qc_reports,  # Create all output files in this specified output directory
	mfiltered_data, mfiltered_data.replace("_R1_", "_R2_"),  # String containing all samples that are about to be checked
	"2>>", os.path.join(pipeline_reports, "preprocessing2_fastQC_frw-report.log")])  # Output fastQC report
	subprocess.run(fastQC_frw, shell=True)

	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  3|4 FastP - Quality Control reports for the input reads are being generated: in progress ..')
	for files in R1list:
		fastP = ' '.join([
		"fastp",  # Call fastQC to quality control all processed data
		"--thread", args.threads,  # Number of threads to use
		"--in1", files,  # Input read1 file
		"--in2", files.replace("_R1_", "_R2_"),  # Input read2 file
		"--disable_adapter_trimming",  # Adapter trimming is disabled
		"--disable_quality_filtering",  # Quality filtering is disabled
		"--disable_length_filtering",  # Length filtering is disabled
		"--overrepresentation_analysis",  # Enable overrepresented sequence analysis
		"--html", f"{prepr_qc_reports}/{os.path.basename(files)[:-9]}_fastp.html",  # Create HTML file in this specified output directory
		"--json", f"{prepr_qc_reports}/{os.path.basename(files)[:-9]}_fastp.json",  # Create json output file in this specified output directory
		"2>>", os.path.join(pipeline_reports, "preprocessing3_fastP-report.log")])  # Output fastP report
		subprocess.run(fastP, shell=True) 

	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  4|4 MultiQC - Summarised Quality Control report for all the input reads is being generated: in progress ..')
	multiQC = " ".join([
	"multiqc",  # Call MultiQC
	"--quiet",  # Print only log warnings
	"--outdir", prepr_qc_reports,  # Create report in the FastQC reports directory
	"--filename", f"summarised_report_{batch}",  # Name of the output report 
	prepr_qc_reports,  # Directory where all FastQC and FastP reports reside
	"2>>", os.path.join(pipeline_reports, "preprocessing4_multiQC-report.log")])  # Output multiQC report
	subprocess.run(multiQC, shell=True)


	# os.system(f'rm {prepr_qc_reports}/*fastqc.zip')  # Removing zipped fastqc files
	os.system(f'mv {prepr_qc_reports}/summarised_report*.html {reports_dir}')  # Moving summary reports in the reports folder
	# os.system(f'rm -r {qc_reports}/*/summarised_report_*data')  # Removing MultiQC temporary folder
	return

def primer_removal(batch, pairedReads):
	""" An initial very mild base quality trimming will be performed. In this step, we are trying to 
	discard very troublesome bases (whose quality is below Q20). That way we remove obvious trash and 
	trying to improve the chances of a proper merge. """
	print(f'\t{datetime.now().strftime("%d.%m.%Y %H:%M")} PREPROCESSING THE INPUT SAMPLES ({batch})')


	batch_prepr_dir = f'{prepr_dir}/{batch}'
	if not os.path.exists(batch_prepr_dir): os.makedirs(batch_prepr_dir)  # Generation of the directories

	for i, read in enumerate(pairedReads, 1):
		# Obtaining the forward and the reverse reads
		forwardRead = read[0]  # Forward
		reverseRead = read[1]  # Reverse

		# Forward and the reverse reads that will be exported as individual files
		forwardRead_output = os.path.join(batch_prepr_dir, os.path.basename(forwardRead).replace([x for x in [".fastq.gz", ".fq.gz"]\
							 if forwardRead.endswith(x)][0], ".fastq.gz"))
		reverseRead_output = os.path.join(batch_prepr_dir, os.path.basename(reverseRead).replace([x for x in [".fastq.gz", ".fq.gz"]\
							 if reverseRead.endswith(x)][0], ".fastq.gz"))
		
		# Obtaining the sample name
		sample_name = os.path.basename(forwardRead).split("_")[0]

		print(f'{i}/{len(pairedReads)}. BBDUK - Quality filtering and primer trimming of {sample_name}')
		# Calculating the minimum and maximum acceptable length after primer and quality trimming
		# The thresholds have been chosen based on https://doi.org/10.1186/s12859-019-3187-5
		# FYI the average length of V3-V4 is 443bp
		minlength, maxlength = calculateMinMax(forwardRead)
		
		bbduk = ' '.join([
		# "/home/stavros/playground/progs/anaconda3/bin/bbduk.sh",  # Call BBDuck (BBTools) to preprocess the raw data
		"bbduk.sh",  # Call BBDuck (BBTools) to preprocess the raw data
		f"threads={args.threads}",  # Set number of threads to use
		f"in={forwardRead}",  # Input of the forward file
		f"in2={reverseRead}",  # Input of the reverse file
		f"out={forwardRead_output}",  # Export edited forward read to file
		f"out2={reverseRead_output}",  # Export edited reverse read to file
		f"literal={args.forwardPrimer},{args.reversePrimer}",  # Providing the forward and reverse primers
		f"minlength={minlength}", #220 Pairs (or reads) will be discarded if both are shorter than this after trimming
		f"maxlength={maxlength}", # ~ 560 Pairs (or reads) will be discarded only if both are longer than this after trimming
		# "copyundefined=t",  # Process non-AGCT IUPAC reference bases by making all possible unambiguous copies
		"ordered=t",  # Keeps the reads in the same order as we gave them to the software
		"mink=2",  # Specifies the smallest word size it will check against either edge of a read
		"rcomp=t",  # Look for reverse-complements of kmers in addition to forward kmers
		"ktrim=l",  # Trim to the left of reads to remove bases matching reference kmers
		"trimq=15",  # Regions with average quality BELOW this will be trimmed
		"qtrim=r",  # Trim read ends (left end only) to remove bases with quality below trimq
		"maq=15",  # Reads with average quality (after trimming) below this will be discarded (15)
		"mm=f",  # Looking for exact kmers and not mask the middle bases
		"k=10",  # Setting the kmer size we want to search for
		"2>>", os.path.join(pipeline_reports, "removeprimers_bbduk-report.log")])  # Output trimming report
		print("\n\n", file=open(os.path.join(pipeline_reports, "removeprimers_bbduk-report.log"),"a"))
		subprocess.run(bbduk, shell=True)
	return 

def denoising(batch):
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} QIIME2 ANALYSIS ({batch})')


	# Creating input_qiime2data directory to save the input data qiime2 artifacts 
	if not os.path.exists(inputqiime_dir): os.makedirs(inputqiime_dir)  


	""" Importing and denoising the preprocessed PE reads """
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Importing the preprocessed reads to the Qiime2 Artifact: in progress ..')
	importingSamplesToQiime2 =	' '.join([
	f"{qiime2_env}/qiime tools import",  # Run QIIME IMPORT to import data and create a new QIIME 2 Artifact
	"--type", "\'SampleData[PairedEndSequencesWithQuality]\'",  # The semantic type of the artifact that will be created upon importing
	"--input-format", "CasavaOneEightSingleLanePerSampleDirFmt",
	"--input-path", f"{prepr_dir}/{batch}",  # Path to the directory that should be imported
	"--output-path", f"{inputqiime_dir}/input_data_{batch}.qza",  # Path where output artifact should be written
	"2>>", os.path.join(pipeline_reports, "denoising1_importingSamplesToQiime2-report.log")])  # Output denoising report
	subprocess.run(importingSamplesToQiime2, shell=True)

	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Generating interactive positional quality plots: in progress ..')
	importSamplesQC =	' '.join([
	f"{qiime2_env}/qiime demux summarize",  # Calling qiime demux summarize to quality of each sample
	"--quiet",  # Silence output if execution is successful
	"--i-data", f"{inputqiime_dir}/input_data_{batch}.qza",  # Path where the input artifact is written
	"--o-visualization", f"{inputqiime_dir}/{batch}_data_QC.qzv",  # Output reports
	"2>>", os.path.join(pipeline_reports, "denoising2_importSamplesQC-report.log")])  # Output importSamplesQC report
	subprocess.run(importSamplesQC, shell=True)
	export(f"{inputqiime_dir}/{batch}_data_QC.qzv")

	""" Denoising is an attempt to correct reads with sequencing errors and then 
	remove chimeric sequences originating from different DNA templates. """
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Denoising, dereplicating and filtering chimera sequences from the paired-end data: in progress ..')
	dada2 = ' '.join([
	f"{qiime2_env}/qiime dada2 denoise-paired",  # Call qiime dada2 to denoise the preprocessed data
	"--p-n-threads", args.threads,  # Number of threads to use
	"--quiet",  # Silence output if execution is successful
	"--p-trunc-len-f", "0",  # No truncation will be performed cause we have already trimmed low quality ends
	"--p-trunc-len-r", "0",  # No truncation will be performed cause we have already trimmed low quality ends
	"--output-dir", f"{denoise_dir}/{batch}",  # Output results to a directory
	"--i-demultiplexed-seqs", f"{inputqiime_dir}/input_data_{batch}.qza",  # The paired-end demultiplexed sequences to be denoised
	"2>>", os.path.join(pipeline_reports, "denoising3_dada2-report.log")])  # Output denoising report
	subprocess.run(dada2, shell=True)

	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Feature table summary: in progress ..')
	featureTableSummary = ' '.join([
	f"{qiime2_env}/qiime feature-table summarize",  # Calling qiime2 feature-table summarize function
	"--quiet",  # Silence output if execution is successful
	"--i-table", f"{denoise_dir}/{batch}/table.qza",  # The feature table to be summarized
	"--m-sample-metadata-file", args.metadata,  # Metadata file
	"--o-visualization", f"{denoise_dir}/{batch}/feature_table.qzv",  # Output results to directory
	"2>>", os.path.join(pipeline_reports, "denoising4_featureTableSummary-report.log")])  # Output featureTableSummary report
	subprocess.run(featureTableSummary, shell=True)
	export(f"{denoise_dir}/{batch}/feature_table.qzv")

	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Visualising the denoising stats: in progress ..')
	denoisingStats = ' '.join([
	f"{qiime2_env}/qiime metadata tabulate",  # Calling qiime2 metadata tabulate
	"--m-input-file", f"{denoise_dir}/{batch}/denoising_stats.qza",  # Input stats file
	"--o-visualization", f"{denoise_dir}/{batch}/denoising_stats.qzv",  # Output results to directory
	"2>>", os.path.join(pipeline_reports, "denoising5_denoisingStats-report.log")])  # Output report
	subprocess.run(denoisingStats, shell=True)
	export(f"{denoise_dir}/{batch}/denoising_stats.qzv")
	return 

def decontam(batch):
	""" Using the statistical method Decontam to identify contaminating DNA features from the 
	Negative Controls and remove them in order to capture a more accurate picture of sampled 
	communities in the metagenomics data."""
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} DECONTAMINATING ({batch})')



	# Running decontam
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Decontam - Running the decontamination task on R: in progress ..')
	decontam = " ".join([
	"Rscript",  # Call Rscript
	f"{rscripts}/decontam.R",  # Calling the diffExpr_ExplAnalysis.R script
	f"{denoise_dir}/{batch}/table.qza",  # Input feature matrix
	args.metadata,  # Input metadata
	f"{denoise_dir}/{batch}",  # Output directory
	"2>>", os.path.join(pipeline_reports, "decontam1_decontam-report.txt")])  # Directory where all reports reside
	subprocess.run(decontam, shell=True)
	
	# Converting the .tsv to .biom
	biomconvert = ' '.join([
	"biom convert",  # Calling biom convert
	"--input-fp", f"{denoise_dir}/{batch}/decontam_filtered_table.tsv",  # The input BIOM table
	"--output-fp", f"{denoise_dir}/{batch}/decontam_filtered_table.biom",
	"--to-hdf5",  # Output to hdf5 format
	"2>>", os.path.join(pipeline_reports, "decontam2_biomconvert-report.log")])  # Output report
	subprocess.run(biomconvert, shell=True)
	
	# Importing the decontam feature matrix
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Decontam - Importing the decontam results to Qiime2: in progress ..')
	importFeatureMat = ' '.join([
	f"{qiime2_env}/qiime tools import",  # Calling qiime tools import
	"--input-path", f"{denoise_dir}/{batch}/decontam_filtered_table.biom",  # Path to representative_sequences.qza file that should be exported
	"--output-path", f"{denoise_dir}/{batch}/decontam_filtered_table.qza",  # Output reports
	"--type \'FeatureTable[Frequency]\'",  # Output 
	"2>>", os.path.join(pipeline_reports, "decontam3_iimportFeatureMat-report.log")])  # Output report
	subprocess.run(importFeatureMat, shell=True)
	
	# Removing all Negative Control samples from the analysis 
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Decontam - Removing the Negative Control samples from the analysis: in progress ..')
	remove_negctrls = ' '.join([
	f"{qiime2_env}/qiime feature-table filter-samples",  # Calling qiime feature-table filter-samples
	"--i-table", f"{denoise_dir}/{batch}/decontam_filtered_table.qza",  # Table containing feature ids used for id-based filtering
	"--m-metadata-file", args.metadata,  # Metadata file
	"--p-where", "\"[SampleType]='Negative control'\"",    # Removing samples indicated as 'Negative controls' in the metadata
	"--p-exclude-ids",  # If true, the samples selected by 'metadata' or 'where' parameters will be excluded from the filtered table
	"--o-filtered-table", f"{denoise_dir}/{batch}/decontam_filtered_table_noNegCtrls.qza",  # The resulting filtered table
	"2>>", os.path.join(pipeline_reports, "decontam4_exportReprSeqs-report.log")])  # Output report
	subprocess.run(remove_negctrls, shell=True)
	
	# Filtering out any features with less than X counts
	min_frequency = min_freq(batch)
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Decontam - Filtering out any features with less than {min_frequency} counts: in progress ..')
	filter_feature_table = ' '.join([
	f"{qiime2_env}/qiime feature-table filter-features",
	"--p-min-frequency", min_frequency,  # Minimum frequency that a feature must have to be retained
	"--i-table", f"{denoise_dir}/{batch}/decontam_filtered_table_noNegCtrls.qza",
	"--o-filtered-table", f"{denoise_dir}/{batch}/decontam_filtered_table_noNegCtrls_noLowlyExprFeatures.qza",
	"2>>", os.path.join(pipeline_reports, "decontam5_filter_feature_table-report.log")])
	subprocess.run(filter_feature_table, shell=True)

	# Exporting the decontam representative sequences 
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Decontam - Applying the decontamination results to the representative sequences: in progress ..')
	exportReprSeqs = ' '.join([
	f"{qiime2_env}/qiime feature-table filter-seqs",  # Calling qiime feature-table filter-seqs
	"--i-data", f"{denoise_dir}/{batch}/representative_sequences.qza",  #  The sequences from which features should be filtered
	"--i-table", f"{denoise_dir}/{batch}/decontam_filtered_table_noNegCtrls_noLowlyExprFeatures.qza",  # Table containing feature ids used for id-based filtering
	"--o-filtered-data", f"{denoise_dir}/{batch}/representative_sequences_filt_new.qza",  # The resulting filtered sequences
	"2>>", os.path.join(pipeline_reports, "decontam6_exportReprSeqs-report.log")])  # Output report
	subprocess.run(exportReprSeqs, shell=True)

	os.system('rm {0}/{1}/*table.biom {0}/{1}/*filtered_table.tsv {0}/{1}/*noNegCtrls.qza'.format(denoise_dir, batch))  # Removing unnecessary files
	return

def merge_denoised_data(batch_to_analyse):
	""" Merging different batches of representative 
	sequences and feature tables """
	if len(batch_to_analyse) > 1:
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Combining all {len(batch_to_analyse)} sets of data feature tables: in progress ..')
		feature_tables = ' '.join(glob.glob(f'{denoise_dir}/*/decontam_filtered_table_noNegCtrls_noLowlyExprFeatures.qza'))
		# Combine two sets of data feature tables
		combFeatureTable =	' '.join([
		f"{qiime2_env}/qiime feature-table merge",  # Calling qiime feature-table merge
		"--i-tables", feature_tables,  # Feature tables to be merged
		"--o-merged-table", f"{denoise_dir}/table.qza",  # Output reports
		"2>>", os.path.join(pipeline_reports, "merge1_combFeatureTable-report.log")])  # Output importSamplesQC report
		subprocess.run(combFeatureTable, shell=True)

		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Combining the representative sequences of all {len(batch_to_analyse)} sets of data: in progress ..')
		sequence_tables = ' '.join(glob.glob(f'{denoise_dir}/*/representative_sequences_filt_new.qza'))
		# Combine the representative sequences of the two sets of data
		combSeqsTable =	' '.join([
		f"{qiime2_env}/qiime feature-table merge-seqs",  # Calling qiime feature-table merge-seqs
		"--i-data", sequence_tables,  # Representative sequences to be merged
		"--o-merged-data", f"{denoise_dir}/representative_sequences.qza",  # Output reports
		"2>>", os.path.join(pipeline_reports, "merge2_combSeqsTable-report.log")])  # Output importSamplesQC report
		subprocess.run(combSeqsTable, shell=True)		
	else:
		subprocess.run(f'mv {denoise_dir}/*/decontam_filtered_table_noNegCtrls_noLowlyExprFeatures.qza {denoise_dir}/table.qza', shell=True)  # Moving table.qza to the denoising directory
		subprocess.run(f'mv {denoise_dir}/*/representative_sequences_filt_new.qza {denoise_dir}/representative_sequences.qza', shell=True)  # Moving representative_sequences.qza to the denoising directory

	
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Final feature table summary: in progress ..')
	finfeatureTableSummary = ' '.join([
	f"{qiime2_env}/qiime feature-table summarize",  # Calling qiime2 feature-table summarize function
	"--quiet",  # Silence output if execution is successful
	"--i-table", f"{denoise_dir}/table.qza",  # The feature table to be summarized
	"--m-sample-metadata-file", args.metadata,  # Metadata file
	"--o-visualization", f"{denoise_dir}/feature_table.qzv",  # Output results to directory
	"2>>", os.path.join(pipeline_reports, "merge3_finfeatureTableSummary-report.log")])  # Output featureTableSummary report
	subprocess.run(finfeatureTableSummary, shell=True)
	export(f"{denoise_dir}/feature_table.qzv")
	return

def taxonomic_assignmnet():
	""" We will train the Naive Bayes classifier using SILVA (132) reference sequences 
	and classify the representative sequences from the input dataset """
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} TAXONOMY ASSIGNMENT')

	# Importing SILVA reference taxonomy sequences
	if not os.path.exists(taxonomy_dir): os.makedirs(taxonomy_dir)  # Creating the directory which will host the analysis
	if not os.path.exists(silva_99_classifier):
		print("WARNING: Pre-trained classifier DOES NOT exist. Training with SILVA db: in progress..")
		
		""" It has been shown that taxonomic classification accuracy of 16S rRNA gene sequences 
		improves when a Naive Bayes classifier is trained on only the region of the target 
		sequences that was sequenced. Here we will extract the reference sequences. """
		# Extract sequencing-like reads from a reference database
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  1/2 | Extract sequencing-like reads from the reference database: in progress ..')
		extractRefReads = ' '.join([
		f"{qiime2_env}/qiime feature-classifier extract-reads",
		"--quiet",  # Silence output if execution is successful
		"--p-n-jobs", args.threads,  # Number of separate processes to run
		"--p-trunc-len", "466",  # Minimum amplicon length
		"--p-f-primer", args.forwardPrimer,  # Forward primer sequence
		"--p-r-primer", args.reversePrimer,  # Reverse primer sequence
		"--i-sequences", silva_reference,  # Input reference seq artifact
		"--o-reads", silva_reference_targeted,  # Output ref sequencing-like reads
		"2>>", os.path.join(pipeline_reports, "taxonomyNB1_extractRefReads-report.log")])  # Output extractRefReads report
		subprocess.run(extractRefReads, shell=True)

		""" We can now train a Naive Bayes classifier as follows, using 
		the reference reads and taxonomy that we just created """
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  2/2 | Training the Naive Bayes classifier using the reference reads: in progress ..')
		trainClassifier = ' '.join([
		f"{qiime2_env}/qiime feature-classifier fit-classifier-naive-bayes",
		"--quiet",  # Silence output if execution is successful
		"--i-reference-reads", silva_reference_targeted,
		"--i-reference-taxonomy", silva_taxonomy,
		"--o-classifier", silva_99_classifier, 
		"2>>", os.path.join(pipeline_reports, "taxonomyNB2_trainClassifier-report.log")])  # Output trainClassifier report
		subprocess.run(trainClassifier, shell=True)
		
	
	""" Assign the taxonomy """
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Assigning the taxonomy to each ASV: in progress ..')
	assignTaxonomy = ' '.join([
	f"{qiime2_env}/qiime feature-classifier classify-sklearn",
	"--quiet",  # Silence output if execution is successful
	"--p-n-jobs", args.threads,  # Number of threads to use
	"--i-classifier", silva_99_classifier,  # Using Silva classifier
	"--i-reads", os.path.join(denoise_dir, "representative_sequences.qza"),  # The output sequences
	"--o-classification", os.path.join(taxonomy_dir, "taxonomic_classification.qza"),
	"2>>", os.path.join(pipeline_reports, "taxonomy1_assignTaxonomy-report.log")])  # Output assignTaxonomy report
	subprocess.run(assignTaxonomy, shell=True)

	outputClassifications = ' '.join([
	f"{qiime2_env}/qiime metadata tabulate",
	"--quiet",  # Silence output if execution is successful
	"--m-input-file", os.path.join(taxonomy_dir, "taxonomic_classification.qza"),
	"--o-visualization", os.path.join(taxonomy_dir, "taxonomic_classification.qzv"),
	"2>>", os.path.join(pipeline_reports, "taxonomy2_outputClassifications-report.log")])  # Output outputClassifications report
	subprocess.run(outputClassifications, shell=True)
	export(os.path.join(taxonomy_dir, "taxonomic_classification.qzv"))

	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Generating the taxonomic barplot: in progress ..')
	barplotOfTaxonomy = ' '.join([
	f"{qiime2_env}/qiime taxa barplot", 
	"--quiet",  # Silence output if execution is successful
	"--i-table", os.path.join(denoise_dir, "table.qza"),
	"--i-taxonomy", os.path.join(taxonomy_dir, "taxonomic_classification.qza"),
	"--m-metadata-file", args.metadata,  # Metadata file
	"--o-visualization", os.path.join(taxonomy_dir, "taxonomy_barplot.qzv"),
	"2>>", os.path.join(pipeline_reports, "taxonomy3_barplotOfTaxonomy-report.log")])  # Output barplotOfTaxonomy report
	subprocess.run(barplotOfTaxonomy, shell=True)
	export(os.path.join(taxonomy_dir, "taxonomy_barplot.qzv"))
	return

class phylogenetics():
	def __init__(self):
		""" This pipeline will start by creating a sequence alignment using MAFFT,
	  	after which any alignment columns that are phylogenetically uninformative
	  	or  ambiguously aligned  will be removed  (masked). The resulting masked
	  	alignment will be used to infer a phylogenetic tree and then subsequently
	  	rooted at its  midpoint. """
		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} PHYLOGENETIC DIVERSITY ANALYSIS')
		if not os.path.exists(diversity_dir): os.makedirs(diversity_dir)  # Creating the directories which will host the analysis
		self.phylogenetic_tree()
		self.alpha_rarefaction()
		self.core_diversity_analysis()
		return

	def max_depth_threshold(self):
		""" Calculating the maximum depth and step for the rarefraction experiment"""
		depth_file = os.path.join(denoise_dir, "feature_table/sample-frequency-detail.csv")
		if not os.path.exists(depth_file):
			sys.exit(f"The file {depth_file} does NOT exist...")

		depths = {}
		with open(depth_file) as fin:
			for line in fin:
				depths[line.split(",")[0].strip()] = int(float(line.split(",")[1].strip()))

		# Ordering dictionary by decreasing value
		depths = sorted(depths.items(), key = operator.itemgetter(1), reverse = True)
		if depths[-1][1] == 0:
			min_depth = str(depths[-2][1])
		else:
			min_depth = str(depths[-1][1])
		return min_depth

	def phylogenetic_tree(self):
		""" The tree provides an inherent structure to the data, allowing us to consider an evolutionary relationship between organisms """
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Generating a rooted tree for phylogenetic diversity analysis: in progress ..')
		phylogeneticdiversity_dir =	' '.join([
		f"{qiime2_env}/qiime phylogeny align-to-tree-mafft-fasttree", # Calling qiime2 align-to-tree-mafft-fasttree function
		"--quiet",  # Silence output if execution is successful
		"--p-n-threads", args.threads,  # Number of threads to use
		"--i-sequences", os.path.join(denoise_dir, "representative_sequences.qza"),  # he sequences to be used for creating a phylogenetic tree
		"--o-alignment", os.path.join(diversity_dir, "aligned_representative_sequences.qza"),  # The aligned sequences
		"--o-masked-alignment", os.path.join(diversity_dir, "masked_aligned_representative_sequences.qza"),  # The masked alignment
		"--o-tree", os.path.join(diversity_dir, "unrooted_tree.qza"),  # The unrooted phylogenetic tree
		"--o-rooted-tree", os.path.join(diversity_dir, "rooted_tree.qza"),  # The rooted phylogenetic tree
		"2>>", os.path.join(pipeline_reports, "phylogenetics1_phylogeneticdiversity_dir-report.log")])  # Output phylogeneticdiversity_dir report
		subprocess.run(phylogeneticdiversity_dir, shell=True)
		return

	def alpha_rarefaction(self):
		""" A key quality control step is to plot rarefaction curves for all 
		the samples to determine if performed sufficient sequencing """
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Performing a rarefraction curves analysis: in progress ..')
		# min_depth, steps = maxDepthNpStepsThresholds()
		rarefactionCurvesAnalysis =	' '.join([
		f"{qiime2_env}/qiime diversity alpha-rarefaction",  # Calling qiime2 diversity alpha-rarefaction function
		"--quiet",  # Silence output if execution is successful
		"--p-max-depth", self.max_depth_threshold(),  # The maximum rarefaction depth
		"--p-steps 50",
		"--m-metadata-file", args.metadata,  # Metadata file
		"--i-table", os.path.join(denoise_dir, "table.qza"),  # Input feature table
		"--i-phylogeny", os.path.join(diversity_dir, "rooted_tree.qza"),  #  Input phylogeny for phylogenetic metrics
		"--o-visualization", os.path.join(diversity_dir, "rarefaction_curves.qzv"),  # Output visualisation
		"2>>", os.path.join(pipeline_reports, "phylogenetics2_rarefactionCurvesAnalysis-report.log")])  # Output rarefactionCurvesAnalysis report
		subprocess.run(rarefactionCurvesAnalysis, shell=True)
		export(os.path.join(diversity_dir, "rarefaction_curves.qzv"))
		return

	def core_diversity_analysis(self):
		""" Common alpha and beta-diversity metrics and ordination plots (such as PCoA plots for weighted UniFrac distances) 
		This command will also rarefy all samples to the sample sequencing depth before calculating these metrics 
		(X is a placeholder for the lowest reasonable sample depth; samples with depth below this cut-off will be excluded) """
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Calculating multiple diversity metrics: in progress ..')
		diversityMetrics =	' '.join([
		f"{qiime2_env}/qiime diversity core-metrics-phylogenetic",  # Calling qiime2 diversity core-metrics-phylogenetic function
		"--quiet",  # Silence output if execution is successful
		"--i-table", os.path.join(denoise_dir, "table.qza"),  # The input feature table
		"--i-phylogeny", os.path.join(diversity_dir, "rooted_tree.qza"),  # The rooted phylogenetic tree
		"--p-sampling-depth", self.max_depth_threshold(),  # The total frequency that each sample should be rarefied to prior to computing diversity metrics
		"--m-metadata-file", args.metadata,  # Metadata file
		"--p-n-jobs-or-threads", args.threads, # The number of CPUs to be used for the computation
		"--output-dir", os.path.join(diversity_dir, "core_metrics_results"),  # Output directory that will host the core metrics
		"2>>", os.path.join(pipeline_reports, "phylogenetics3_diversityMetrics-report.log")])  # Output diversityMetrics report
		subprocess.run(diversityMetrics, shell=True)
		export(os.path.join(diversity_dir, "core_metrics_results"))
		return

def downstream_analysis():
	
	if "balf" in args.input_dir:
		downstream_rscript = os.path.join(rscripts, "downstream_analysis_balf.R")
		accepted_protocol = "standard"
	elif "blood" in args.input_dir:
		downstream_rscript = os.path.join(rscripts, "downstream_analysis_blood.R")
		accepted_protocol = "beads"


	if not os.path.exists(downstream_dir): os.makedirs(downstream_dir)

	# Running downstream analysis
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Downstream - Running additional downstream analysis on R: in progress ..')
	danalysis = " ".join([
	"Rscript",  # Call Rscript
	downstream_rscript,  # Calling the necessary R script
	f"{denoise_dir}/table.qza",  # The input feature table
	f"{taxonomy_dir}/taxonomic_classification.qza",  # Input taxonomic classification
	f"{diversity_dir}/rooted_tree.qza",  #  Input phylogenetic rooted tree
	args.metadata,  # Metadata file
	downstream_dir,  # Output dir
	"nonCancer,cancer",  # Groups for Differential Abundance Testing
	accepted_protocol,  # Accepted protocol
	args.threads,  # Number of cores to use
	# "2>>", os.path.join(pipeline_reports, "downstream_analysis-report.txt")
	])  # Directory where all reports reside
	subprocess.run(danalysis, shell=True)
	return

def ml_approach():

	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} MACHINE LEARNING APPROACHES')


	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Random Forst classifiers for predicting sample characteristics: in progress ..')
	random_forest = ' '.join([
	f"{qiime2_env}/qiime sample-classifier classify-samples",
	"--i-table", os.path.join(denoise_dir, "table.qza"),
	"--m-metadata-file", args.metadata,
	"--m-metadata-column", "Indication_I",
	"--p-random-state", "500",
	# "--p-n-jobs", str(args.threads),
	"--output-dir", random_forest_dir,
	# "2>>", os.path.join(pipeline_reports, "qiime2_random_forest-report.log")
	])  # Output denoising report
	# subprocess.run(random_forest, shell=True)

	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Generating heatmap: in progress ..')
	heatmap = ' '.join([
	f"{qiime2_env}/qiime sample-classifier heatmap",
	"--i-table", os.path.join(denoise_dir, "table.qza"),
	"--i-importance", os.path.join(random_forest_dir, "feature_importance.qza"),
	"--m-sample-metadata-file", args.metadata,
	"--m-sample-metadata-column", "Indication_I",
	"--p-group-samples",
	"--p-cluster", "samples",
	"--p-feature-count", "50",
	"--p-color-scheme", "viridis",
	"--o-heatmap", os.path.join(random_forest_dir, "heatmap.qzv"),
	"--o-filtered-table", os.path.join(random_forest_dir, "filtered-table.qza"),
	# "2>>", os.path.join(pipeline_reports, "qiime2_heatmap-report.log")
	])  # Output denoising report
	subprocess.run(heatmap, shell=True)
	export(os.path.join(random_forest_dir, "heatmap.qzv"))
	return

def calculateMinMax(forwardRead):
	primers_averageLength = int((len(args.forwardPrimer) + len(args.reversePrimer))/2)

	read_length = 0
	with gzip.open(forwardRead, "rt") as handle:
	    for read in SeqIO.parse(handle, "fastq"):
	        read_length = len(read.seq)
	        break

	min_readLength = read_length - (primers_averageLength * 4)
	max_readLength = (read_length * 2) - (primers_averageLength * 2)
	return (min_readLength, max_readLength)

def min_freq(batch):
	""" Based on the summary we will calculate a cut-off for how frequent a variant needs to be for it to be retained. 
	We will remove all ASVs that have a frequency of less than 0.001% of the mean sample depth. This cut-off excludes """
	exportFile = f"{denoise_dir}/{batch}/feature_table.qzv"
	if not os.path.exists(exportFile[:-4]):
		subprocess.run(f"{qiime2_env}/qiime tools export --input-path {exportFile} --output-path {exportFile[:-4]}", shell=True)

	mean_freq = 0
	for path, subdir, folder in os.walk(exportFile[:-4]):
		for name in folder:
			file = os.path.join(path, name)
			if file.endswith("sample-frequency-detail.csv"):
				samples = 0
				with open(file) as fin:
					for line in fin:
						samples += 1
						mean_freq += float(line.split(",")[1])
				mean_freq = mean_freq/(samples-1)

	# Calculation of the cut-off point
	cut_off = str(int(mean_freq * 0.001))
	return cut_off

def export(exportFile):
	if os.path.isfile(exportFile):
		subprocess.run(f"{qiime2_env}/qiime tools export --input-path {exportFile} --output-path {exportFile[:-4]}", shell=True)
	else:
		for path, subdir, folder in os.walk(exportFile):
			for files in folder:
				file = os.path.join(path, files)
				subprocess.run(f"{qiime2_env}/qiime tools export --input-path {file} --output-path {file[:-4]}", shell=True)

	if exportFile.endswith(".qzv"):
		os.system(f'rm {exportFile}')	
	return 

def summarisation():
	""" Finalising the metagenomics analysis by removing all the 
	unnecessary files and reporting log files containing several 
	comments. """
	for path, subdir, folder in os.walk(analysis_dir):
		for name in folder:
			file = os.path.join(path, name)
			if os.stat(file).st_size == 0:  
				print("Removing:\t", file)  
				os.remove(file)
			elif file.endswith("input_data.qza"):
				print("Removing:\t", file)
				os.remove(file)
			elif file.endswith('bbduk-report.log') or file.endswith('fastP-report.log')\
				or file.endswith('fastQscreen-report.log'):
				os.system(f'mv {file} {qc_reports}')
	
	num = 0
	# Mentioning the remaining log files which have warnings!
	for files in os.listdir(reports_dir):
		if files.endswith('-report.log'):
			num+=1
	print(f'ATTENTION -- {num} reports contain warnings!!!')
	# shutil.rmtree(prepr_dir)  # removing the dir containing the preprocessed fastq files
	return 


def main():

	print(f'\t{datetime.now().strftime("%d.%m.%Y %H:%M")} METAGENOMICS ANALYSIS')
	print(f'\n- In total samples from {len(batch_to_analyse)} batches will be processed in this analysis..')
	
	for batch in batch_to_analyse:
		### if batch == 'batch14':
		pairedReads, R1list = assess_inputdata(batch)  # Assessing all input fastq files and obtaining the pairs to be analysed
		print(f'- In total {len(R1list)} samples will be processed from {batch}..')
		
		
		# quality_control(batch, R1list)  # Checking the quality of the reads

		# # Preprocessing of the input data
		# primer_removal(batch, pairedReads)  # Performing quality trimming and removal of all primers on both reads
		
		# denoising(batch)  # Qiime2 analysis - denoising

		# decontam(batch)  # Performing decontamination
	
	# merge_denoised_data(batch_to_analyse)

	# taxonomic_assignmnet()

	# phylogenetics()

	# downstream_analysis()

	# ml_approach()

	# summarisation()

if __name__ == "__main__": main()