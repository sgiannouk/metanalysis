###Stavros Giannoukakos###
##University of Granada##

# To activate this environment, use
#
#     $ conda activate qiime2
#
# To deactivate an active environment, use
#
#     $ conda deactivate

#Version of the program
__version__ = "0.2.0"

import argparse
import operator
from Bio import SeqIO
import subprocess, gzip
from Bio.Seq import Seq
from datetime import datetime
import shutil, fnmatch, glob, sys, os

########### TEST ############
batch_to_analyse = ["batch14",  "batch16"]
inputDir = '/shared/projects/martyna_rrna_illumina_igtp/test_data'
metadata_file = "/shared/projects/martyna_rrna_illumina_igtp/test_data/clinical_data.txt"
#############################

# batch_to_analyse = ["batch2",  "batch3"]
# inputDir = '/shared/projects/martyna_rrna_illumina_igtp'
# metadata_file = "/shared/projects/martyna_rrna_illumina_igtp/clinical_data.txt"
# Configuration file needed for FastQ Screen
fastQscreen_config = "/home/stavros/references/fastQscreen_references/fastq_screen.conf"
# SILVA database
silva_reference = "/home/stavros/playground/progs/16S_subsidiary_files/SILVA_132/rep_set/rep_set_16S_only/99/silva_132_99_16S.fna"
silva_taxinomy = "/home/stavros/playground/progs/16S_subsidiary_files/SILVA_132/taxonomy/16S_only/99/taxonomy_7_levels.txt"
silva_99_classifier = "/home/stavros/playground/progs/16S_subsidiary_files/silva_132_99_v3v4_scikitv0.21.2.qza"  # scikit-learn version 0.21.2.
# silva_99_classifier = "/home/stavros/playground/progs/16S_subsidiary_files/silva_132_99_v3v4.qza"  # scikit-learn version 0.18.0.


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
parser.add_argument('-th', '--threads', dest='threads', default=str(70), metavar='', 
                	help="Number of threads to be used in the analysis")
# Number of threads/CPUs to be used
parser.add_argument('-fp', '--forwardPrimer', dest='forwardPrimer', default="CCTACGGGNGGCWGCAG", metavar='', required=False, 
                	help="Sequence of the forward primer")
# Number of threads/CPUs to be used
parser.add_argument('-rp', '--reversePrimer', dest='reversePrimer', default="GACTACHVGGGTATCTAATCC", metavar='', required=False,
                	help="Sequence of the reverse primer")
# Metadata file
parser.add_argument('-m', '--metadata', dest='metadata', required=False, metavar='', 
                	help="Metadata file containing several info conserning the\ninput data")
# Path of where the output folder should be located
parser.add_argument('-o', '--output_dir', dest='output_dir', default=os.path.dirname(os.path.realpath(__file__)), 
					help="Path of the directory that will host the analysis.\n(default: <current_directory>)")
# Display the version of the pipeline 
parser.add_argument('-v', '--version', action='version', version='%(prog)s {0}'.format(__version__))

# Get the options and return them
args = parser.parse_args()
start_time = datetime.now()  # Tracking time of analysis

""" All the necessary directory that will host the analysis are being created. These 
include 'preprocessed_files' that will host the filtered and quality controlled 
data and reports, where all reports from all software will be stored. """
args.metadata = metadata_file
args.input_dir = inputDir

# Main directories hosting the analysis
analysis_dir = os.path.join(args.output_dir, "test")
reports_dir = os.path.join(analysis_dir, "reports")  # Reports directory
prepr_dir = os.path.join(analysis_dir, "preprocessed_data")  # Save processed .fastq files
inputqiime_dir = os.path.join(analysis_dir, "qiime2_input_data")  # Input data directory
denoise_dir = os.path.join(analysis_dir, "denoising_analysis")  # Denoising output
taxonomy_dir = os.path.join(analysis_dir, "taxonomic_analysis")  # Taxonomic analysis
diversity_dir = os.path.join(analysis_dir, "diversity_analysis")  # Diversity analysis

adonis_dir = os.path.join(diversity_dir, "adonis_analysis")
statistical_analysis_dir = os.path.join(diversity_dir, "statistical_analysis")

abundanceAnalysis = os.path.join(analysis_dir, "differential_abundance_analysis")
biplots = os.path.join(analysis_dir, "PCoA_biplot_analysis")
random_forest_dir = os.path.join(analysis_dir, "random_forest")

# Reporting directories
qc_reports = os.path.join(analysis_dir, "reports/QC_reports")
pipeline_reports = os.path.join(analysis_dir, "reports/pipeline_reports")




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
				assert (os.path.isfile(inR2)), 'Could not detect the pair of {0} ({1})'.format(inR1, inR2)
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
	if not os.path.exists(pipeline_reports): os.makedirs(pipeline_reports)
	if not os.path.exists(prepr_qc_reports): os.makedirs(prepr_qc_reports)
	


	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  1|4 FastQScreen - Checking random reads for possible contamination: in progress ..')
	fastQscreen = ' '.join([
	"fastq_screen",  # Call fastQ screen to check contamination in the processed data
	"--threads", args.threads,  # Number of threads to use
	"--outdir",  prepr_qc_reports,  # Directory in which the output files will be saved
	"--quiet",  # Suppress all progress reports on stderr and only report errors
	"--conf", fastQscreen_config,  # Location of the required configuration file
	mfiltered_data, mfiltered_data.replace("_R1_", "_R2_"),  # Input PE files
	"2>>", os.path.join(pipeline_reports, "1_preprocessing_fastQscreen-report.log")])  # Output fastQ screen report
	subprocess.run(fastQscreen, shell=True)

	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  2|4 FastQC - Quality Control reports for the forward reads are being generated: in progress ..')
	fastQC_frw = ' '.join([
	"fastqc",  # Call fastQC to quality control all processed data
	"--threads", args.threads,  # Number of threads to use
	"--quiet",  # Print only log warnings
	"--outdir", prepr_qc_reports,  # Create all output files in this specified output directory
	mfiltered_data, mfiltered_data.replace("_R1_", "_R2_"),  # String containing all samples that are about to be checked
	"2>>", os.path.join(pipeline_reports, "2_preprocessing_fastQC_frw-report.log")])  # Output fastQC report
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
		"2>>", os.path.join(pipeline_reports, "3_preprocessing_fastP-report.log")])  # Output fastP report
		subprocess.run(fastP, shell=True) 

	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  4|4 MultiQC - Summarised Quality Control report for all the input reads is being generated: in progress ..')
	multiQC = " ".join([
	"multiqc",  # Call MultiQC
	"--quiet",  # Print only log warnings
	"--outdir", prepr_qc_reports,  # Create report in the FastQC reports directory
	"--filename", f"summarised_report_{batch}",  # Name of the output report 
	prepr_qc_reports,  # Directory where all FastQC and FastP reports reside
	"2>>", os.path.join(pipeline_reports, "4_preprocessing_multiQC-report.log")])  # Output multiQC report
	subprocess.run(multiQC, shell=True)


	os.system(f'rm {prepr_qc_reports}/*/*fastqc.zip')  # Removing zipped fastqc files
	os.system(f'mv {prepr_qc_reports}/summarised_report*.html {reports_dir}')  # Moving summary reports in the reports folder
	# os.system(f'rm -r {qc_reports}/*/summarised_report_*data'.format())  # Removing MultiQC temporary folder
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
		"2>>", os.path.join(pipeline_reports, "removepr_bbduk-report.log")])  # Output trimming report
		print("\n\n", file=open(os.path.join(pipeline_reports, "removepr_bbduk-report.log"),"a"))
		subprocess.run(bbduk, shell=True)
	return 

def denoising(batch):
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} QIIME2 ANALYSIS ({batch})')


	# Creating input_qiime2data direcotry to save the input data qiime2 artifacts 
	if not os.path.exists(inputqiime_dir): os.makedirs(inputqiime_dir)  

	""" Importing and denoising the preprocessed PE reads """
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} Importing the preprocessed reads to the Qiime2 Artifact: in progress ..')
	importingSamplesToQiime2 =	' '.join([
	"qiime tools import",  # Run QIIME IMPORT to import data and create a new QIIME 2 Artifact
	"--type", "\'SampleData[PairedEndSequencesWithQuality]\'",  # The semantic type of the artifact that will be created upon importing
	"--input-format", "CasavaOneEightSingleLanePerSampleDirFmt",
	"--input-path", f"{prepr_dir}/{batch}",  # Path to the directory that should be imported
	"--output-path", f"{inputqiime_dir}/input_data_{batch}.qza",  # Path where output artifact should be written
	"2>>", os.path.join(pipeline_reports, "1_denoising_importingSamplesToQiime2-report.log")])  # Output denoising report
	subprocess.run(importingSamplesToQiime2, shell=True)

	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} Generating interactive positional quality plots: in progress ..')
	importSamplesQC =	' '.join([
	"qiime demux summarize",  # Calling qiime demux summarize to quality of each sample
	"--quiet",  # Silence output if execution is successful
	"--i-data", f"{inputqiime_dir}/input_data_{batch}.qza",  # Path where the input artifact is written
	"--o-visualization", f"{inputqiime_dir}/{batch}_data_QC.qzv",  # Output reports
	"2>>", os.path.join(pipeline_reports, "2_denoising_importSamplesQC-report.log")])  # Output importSamplesQC report
	subprocess.run(importSamplesQC, shell=True)
	export(f"{inputqiime_dir}/{batch}_data_QC.qzv")

	""" Denoising is an attempt to correct reads with sequencing errors and then 
	remove chimeric sequences originating from different DNA templates. """
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} Denoising, dereplicating and filtering chimera sequences from the paired-end data: in progress ..')
	dada2 = ' '.join([
	"qiime dada2 denoise-paired",  # Call qiime dada2 to denoise the preprocessed data
	"--p-n-threads", args.threads,  # Number of threads to use
	"--quiet",  # Silence output if execution is successful
	"--p-trunc-len-f", "0",  # No truncation will be performed cause we have already trimmed low quality ends
	"--p-trunc-len-r", "0",  # No truncation will be performed cause we have already trimmed low quality ends
	"--output-dir", f"{denoise_dir}/{batch}",  # Output results to a directory
	"--i-demultiplexed-seqs", f"{inputqiime_dir}/input_data_{batch}.qza",  # The paired-end demultiplexed sequences to be denoised
	"2>>", os.path.join(pipeline_reports, "3_denoising_dada2-report.log")])  # Output denoising report
	subprocess.run(dada2, shell=True)

	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} Feature table summary: in progress ..')
	featureTableSummary = ' '.join([
	"qiime feature-table summarize",  # Calling qiime2 feature-table summarize function
	"--m-sample-metadata-file", args.metadata,  # Metadata file
	"--quiet",  # Silence output if execution is successful
	"--i-table", f"{denoise_dir}/{batch}/table.qza",  # The feature table to be summarized
	"--o-visualization", f"{denoise_dir}/{batch}/feature_table.qzv",  # Output results to directory
	"2>>", os.path.join(pipeline_reports, "4_denoising_featureTableSummary-report.log")])  # Output featureTableSummary report
	subprocess.run(featureTableSummary, shell=True)
	export(f"{denoise_dir}/{batch}/feature_table.qzv")

	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} Visualising the denoising stats: in progress ..')
	denoisingStats = ' '.join([
	"qiime metadata tabulate",  # Calling qiime2 metadata tabulate
	"--m-input-file", f"{denoise_dir}/{batch}/denoising_stats.qza",  # Input stats file
	"--o-visualization", f"{denoise_dir}/{batch}/denoising_stats.qzv",  # Output results to directory
	"2>>", os.path.join(pipeline_reports, "5_denoising_denoisingStats-report.log")])  # Output denoisingVisualStats report
	subprocess.run(denoisingStats, shell=True)
	export(f"{denoise_dir}/{batch}/denoising_stats.qzv")
	return 

def decontam(batch):
	""" Using the statistical method Decontam to identify contaminating DNA features from the 
	Negative Controls and remove them in order to capture a more accurate picture of sampled 
	communities in the metagenomics data."""

	# We will export the represena
	return

def merge_denoised_data(batch_to_analyse):
	""" Merging different batches of representative 
	sequences and feature tables """
	if len(batch_to_analyse) > 1:
		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} Combining all {len(batch_to_analyse)} sets of data feature tables: in progress ..')
		feature_tables = ' '.join(glob.glob(f'{denoise_dir}/*/table.qza'))
		# Combine two sets of data feature tables
		combFeatureTable =	' '.join([
		"qiime feature-table merge",  # Calling qiime feature-table merge
		"--i-tables", feature_tables,  # Feature tables to be merged
		"--o-merged-table", f"{denoise_dir}/table.qza",  # Output reports
		"2>>", os.path.join(pipeline_reports, "1_merge_combFeatureTable-report.log")])  # Output importSamplesQC report
		# subprocess.run(combFeatureTable, shell=True)
		
		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} Combining the representative sequences of all {len(batch_to_analyse)} sets of data: in progress ..')
		sequence_tables = ' '.join(glob.glob(f'{denoise_dir}/*/representative_sequences.qza'))
		# Combine the representative sequences of the two sets of data
		combSeqsTable =	' '.join([
		"qiime feature-table merge-seqs",  # Calling qiime feature-table merge-seqs
		"--i-data", sequence_tables,  # Representative sequences to be merged
		"--o-merged-data", f"{denoise_dir}/representative_sequences.qza",  # Output reports
		"2>>", os.path.join(pipeline_reports, "2_qiime2_combSeqsTable-report.log")])  # Output importSamplesQC report
		subprocess.run(combSeqsTable, shell=True)
	else:
		os.system(f'mv {denoise_dir}/*/table.qza {denoise_dir}')  # Moving table.qza to the denoising directory
		os.system(f'mv {denoise_dir}/*/representative_sequences.qza {denoise_dir}')  # Moving representative_sequences.qza to the denoising directory
	return

def taxonomic_assignmnet():
	""" We will train the Naive Bayes classifier using SILVA (132) reference sequences 
	and classify the representative sequences from the input dataset """
	print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} TAXONOMY ASSIGNMENT')

	# Importing SILVA reference taxonomy sequences
	if not os.path.exists(taxonomy_dir): os.makedirs(taxonomy_dir)  # Creating the directory which will host the analysis
	if not os.path.exists(silva_99_classifier):
		print("Pre-trained classirier DOES NOT exist. Training SILVA db: in progress..")
		
		print("1/3 | Importing the reference sequences and the corresponding taxonomic classifications of SILVA123 database: in progress ..")
		importSilvaReference = ' '.join([
		"qiime tools import",  # Import function
		"--type", "\'FeatureData[Sequence]\'",  # Type of imported data
	  	"--input-path", silva_reference,  # Input SILVA 132 database
	  	"--output-path", os.path.join(taxonomy_dir, "silva132_99_ASVs.qza"),  # Output file
	  	"2>>", os.path.join(pipeline_reports, "1_taxonomyNB_importSilvaReference-report.log")])  # Output importSilvaReference report
		subprocess.run(importSilvaReference, shell=True)

		# Importing SILVA reference taxonomy annotation
		importSilvaRefTaxonomy = ' '.join([
		"qiime tools import",  # Import function
		"--type", "\'FeatureData[Taxonomy]\'",  # Type of imported data
		"--input-format", "HeaderlessTSVTaxonomyFormat",  # Type of input file
	  	"--input-path", silva_taxinomy,  # Input annotation file
	  	"--output-path", os.path.join(taxonomy_dir, "silva132_99_ASV_taxonomy.qza"),  # Output artifact 
		"2>>", os.path.join(pipeline_reports, "1_taxonomyNB_importSilvaRefTaxonomy-report.log")])  # Output importSilvaRefTaxonomy report
		subprocess.run(importSilvaRefTaxonomy, shell=True)

		""" It has been shown that taxonomic classification accuracy of 16S rRNA gene sequences 
		improves when a Naive Bayes classifier is trained on only the region of the target 
		sequences that was sequenced. Here we will extract the reference sequences. """
		# Extract sequencing-like reads from a reference database
		print("2/3 | Extract sequencing-like reads from the reference database: in progress ..")
		extractRefReads = ' '.join([
		"qiime feature-classifier extract-reads",
		"--quiet",  # Silence output if execution is successful
		"--p-trunc-len", "466",  # Minimum amplicon length
		"--p-f-primer", args.forwardPrimer,  # Forward primer sequence
		"--p-r-primer", args.reversePrimer,  # Reverse primer sequence
		"--i-sequences", os.path.join(taxonomy_dir, "silva132_99_ASVs.qza"),  # Input reference seq artifact
		"--o-reads", os.path.join(taxonomy_dir, "silva132_reference_sequences.qza"),  # Output ref sequencing-like reads
		"2>>", os.path.join(pipeline_reports, "2_taxonomyNB_extractRefReads-report.log")])  # Output extractRefReads report
		subprocess.run(extractRefReads, shell=True)

		""" We can now train a Naive Bayes classifier as follows, using 
		the reference reads and taxonomy that we just created """
		print("3/3 | Training the Naive Bayes classifier using the reference reads: in progress ..")
		trainClassifier = ' '.join([
		"qiime feature-classifier fit-classifier-naive-bayes",
		"--quiet",  # Silence output if execution is successful
		"--i-reference-reads", os.path.join(taxonomy_dir, "silva132_reference_sequences.qza"),
		"--i-reference-taxonomy", os.path.join(taxonomy_dir, "silva132_99_ASV_taxonomy.qza"),
		"--o-classifier", silva_99_classifier, 
		"2>>", os.path.join(pipeline_reports, "3_taxonomyNB_trainClassifier-report.log")])  # Output trainClassifier report
		subprocess.run(trainClassifier, shell=True)
		export(os.path.join(taxonomy_dir, "classifier.qza"))
		# Deleting secondry unnecessary files
		os.system("rm {0} {1}".format(os.path.join(taxonomy_dir, "silva132_99_ASVs.qza"),\
									  os.path.join(taxonomy_dir, "silva132_99_ASV_taxonomy.qza"),\
									  os.path.join(taxonomy_dir, "silva132_reference_sequences.qza")))
	
	""" Assign the taxonomy """
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Assigning the taxonomy to each ASV: in progress ..')
	assignTaxonomy = ' '.join([
	"qiime feature-classifier classify-sklearn",
	"--quiet",  # Silence output if execution is successful
	"--p-n-jobs", args.threads,  # Number of threads to use
	"--i-classifier", silva_99_classifier,  # Using Silva classifier
	"--i-reads", os.path.join(denoise_dir, "representative_sequences.qza"),  # The output sequences
	"--o-classification", os.path.join(taxonomy_dir, "taxonomic_classification.qza"),
	"2>>", os.path.join(pipeline_reports, "1_taxonomy_assignTaxonomy-report.log")])  # Output assignTaxonomy report
	subprocess.run(assignTaxonomy, shell=True)

	outputClassifications = ' '.join([
	"qiime metadata tabulate",
	"--quiet",  # Silence output if execution is successful
	"--m-input-file", os.path.join(taxonomy_dir, "taxonomic_classification.qza"),
	"--o-visualization", os.path.join(taxonomy_dir, "taxonomic_classification.qzv"),
	"2>>", os.path.join(pipeline_reports, "2_taxonomy_outputClassifications-report.log")])  # Output outputClassifications report
	subprocess.run(outputClassifications, shell=True)
	export(os.path.join(taxonomy_dir, "taxonomic_classification.qzv"))

	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Generating the taxonimic barplot: in progress ..')
	barplotOfTaxonomy = ' '.join([
	"qiime taxa barplot", 
	"--quiet",  # Silence output if execution is successful
	"--i-table", os.path.join(denoise_dir, "table.qza"),
	"--i-taxonomy", os.path.join(taxonomy_dir, "taxonomic_classification.qza"),
	"--m-metadata-file", args.metadata,  # Metadata file
	"--o-visualization", os.path.join(taxonomy_dir, "taxonomy_barplot.qzv"),
	"2>>", os.path.join(pipeline_reports, "3_taxonomy_barplotOfTaxonomy-report.log")])  # Output barplotOfTaxonomy report
	subprocess.run(barplotOfTaxonomy, shell=True)
	export(os.path.join(taxonomy_dir, "taxonomy_barplot.qzv"))

	# Filtering out any features with less than X counts
	min_frequency = min_freq()
	print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")}  Filtering out any features with less than {min_frequency} counts: in progress ..')
	filter_feature_table = ' '.join([
	"qiime feature-table filter-features",
	"--p-min-frequency", min_frequency,  # Minimum frequency that a feature must have to be retained
	"--i-table", os.path.join(denoise_dir, "table.qza"),
	"--o-filtered-table", os.path.join(denoise_dir, "table_filtered.qza"),
	"2>>", os.path.join(pipeline_reports, "4_taxonomy_filter_feature_table-report.log")])
	subprocess.run(filter_feature_table, shell=True)
	return

class phylogenetics():
	def __init__(self, threads, metadata):
		""" This pipeline will start by creating a sequence alignment using MAFFT,
	  	after which any alignment columns that are phylogenetically uninformative
	  	or  ambiguously aligned  will be removed  (masked). The resulting masked
	  	alignment will be used to infer a phylogenetic tree and then subsequently
	  	rooted at its  midpoint. """
		print(f'\n\t{datetime.now().strftime("%d.%m.%Y %H:%M")} PHYLOGENETIC DIVERSITY ANALYSIS')
		if not os.path.exists(diversity_dir): os.makedirs(diversity_dir)  # Creating the directories which will host the analysis
		if not os.path.exists(statistical_analysis_dir): os.makedirs(statistical_analysis_dir)
		if not os.path.exists(adonis_dir): os.makedirs(adonis_dir)
		self.phylogenetic_tree()
		self.rarefication()
		self.alpha_rarefaction()
		self.core_diversity_analysis()
		self.alpha_diversity_analysis()
		self.beta_diversity_analysis()
		return

	def max_depth_and_steps_thresholds(self):
		""" Calculating the maximum depth and step for the rarefraction experiment"""
		depth_file = os.path.join(denoise_dir, "feature_table/sample-frequency-detail.csv")
		if not os.path.exists(depth_file):
			sys.exit("The file {0} does NOT exist...".format(depth_file))

		depths = {}
		with open(depth_file) as fin:
			for line in fin:
				depths[line.split(",")[0].strip()] = int(float(line.split(",")[1].strip()))
		
		# Ordering dictionary by decreasing value
		depths = sorted(depths.items(), key = operator.itemgetter(1), reverse = True)
		if "neg" in depths[-1][0].lower():
			if "neg" in depths[-2][0].lower():
				min_depth = depths[-3][1]
			else:
				min_depth = depths[-2][1]
		else:
			min_depth = depths[-1][1]
		
		min_depth = str(depths[-1][1])
		# steps = str(int(float(min_depth)/50))
		return min_depth

	def phylogenetic_tree(self):
		""" The tree provides an inherent structure to the data, allowing us to consider an evolutionary relationship between organisms """
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Generating a rooted tree for phylogenetic diversity analysis: in progress ..')
		phylogeneticdiversity_dir =	' '.join([
		"qiime phylogeny align-to-tree-mafft-fasttree", # Calling qiime2 align-to-tree-mafft-fasttree function
		"--quiet",  # Silence output if execution is successful
		"--p-n-threads", args.threads,  # Number of threads to use
		"--i-sequences", os.path.join(denoise_dir, "representative_sequences.qza"),  # he sequences to be used for creating a phylogenetic tree
		"--o-alignment", os.path.join(diversity_dir, "aligned_representative_sequences.qza"),  # The aligned sequences
		"--o-masked-alignment", os.path.join(diversity_dir, "masked_aligned_representative_sequences.qza"),  # The masked alignment
		"--o-tree", os.path.join(diversity_dir, "unrooted_tree.qza"),  # The unrooted phylogenetic tree
		"--o-rooted-tree", os.path.join(diversity_dir, "rooted_tree.qza"),  # The rooted phylogenetic tree
		"2>>", os.path.join(pipeline_reports, "1_phylogenetics_phylogeneticdiversity_dir-report.log")])  # Output phylogeneticdiversity_dir report
		subprocess.run(phylogeneticdiversity_dir, shell=True)
		return

	def rarefication(self):
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Performing random subsampling ({min_depth[0]}): in progress ..')
		randomSubsampling = ' '.join([
		"qiime feature-table subsample",  # Calling qiime2 feature-table subsample
		"--i-table", os.path.join(denoise_dir, "table.qza"),  # The feature table to be subsampled
		"--p-axis", "\'feature\'",  # A random set of features will be selected to be retained
		"--p-subsampling-depth", self.max_depth_and_steps_thresholds(),  # The total number of features to be randomly sampled
		"--o-sampled-table", os.path.join(denoise_dir, "subsampled_table.qza"),  # Output results 
		"2>>", os.path.join(pipeline_reports, "2_phylogenetics_randomSubsampling-report.log")])  # Output randomSubsampling report
		subprocess.run(randomSubsampling, shell=True)
		return

	def alpha_rarefaction(self):
		""" A key quality control step is to plot rarefaction curves for all 
		the samples to determine if performed sufficient sequencing """
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Performing a rarefraction curves analysis: in progress ..')
		# min_depth, steps = maxDepthNpStepsThresholds()
		rarefactionCurvesAnalysis =	' '.join([
		"qiime diversity alpha-rarefaction",  # Calling qiime2 diversity alpha-rarefaction function
		"--quiet",  # Silence output if execution is successful
		"--p-max-depth", self.max_depth_and_steps_thresholds(),  # The maximum rarefaction depth
		"--p-steps 50",
		"--m-metadata-file", args.metadata,  # Metadata file
		"--i-table", os.path.join(denoise_dir, "table.qza"),  # Input feature table
		"--i-phylogeny", os.path.join(diversity_dir, "rooted_tree.qza"),  #  Input phylogeny for phylogenetic metrics
		"--o-visualization", os.path.join(diversity_dir, "rarefaction_curves.qzv"),  # Output visualisation
		"2>>", os.path.join(pipeline_reports, "3_phylogenetics_rarefactionCurvesAnalysis-report.log")])  # Output rarefactionCurvesAnalysis report
		subprocess.run(rarefactionCurvesAnalysis, shell=True)
		export(os.path.join(diversity_dir, "rarefaction_curves.qzv"))
		return

	def core_diversity_analysis(self):
		""" Common alpha and beta-diversity metrics and ordination plots (such as PCoA plots for weighted UniFrac distances) 
		This command will also rarefy all samples to the sample sequencing depth before calculating these metrics 
		(X is a placeholder for the lowest reasonable sample depth; samples with depth below this cut-off will be excluded) """
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Calculate multiple diversity metrics: in progress ..')
		diversityMetrics =	' '.join([
		"qiime diversity core-metrics-phylogenetic",  # Calling qiime 2diversity core-metrics-phylogenetic function
		"--m-metadata-file", args.metadata,  # Metadata file
		"--quiet",  # Silence output if execution is successful
		"--p-n-jobs", args.threads, # The number of CPUs to be used for the computation
		"--i-phylogeny", os.path.join(diversity_dir, "rooted_tree.qza"),  # The rooted phylogenetic tree
		"--i-table", os.path.join(denoise_dir, "table.qza"),  # The input feature table
		"--p-sampling-depth", self.max_depth_and_steps_thresholds(),  # The total frequency that each sample should be rarefied to prior to computing diversity metrics
		"--output-dir", os.path.join(diversity_dir, "core_metrics_results"),  # Output directory that will host the core metrics
		"2>>", os.path.join(pipeline_reports, "4_phylogenetics_diversityMetrics-report.log")])  # Output diversityMetrics report
		subprocess.run(diversityMetrics, shell=True)
		export(os.path.join(diversity_dir, "core_metrics_results"))
		return

	def alpha_diversity_analysis(self):
		""" Visually and statistically compare groups of alpha diversity values """
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Creating Shannon diversity significance boxplot: in progress ..')
		shannonSignificance =	' '.join([
		"qiime diversity alpha-group-significance",  # Calling qiime diversity alpha-group-significance function
		"--m-metadata-file", args.metadata,  # Metadata file
		"--quiet",  # Silence output if execution is successful
		"--i-alpha-diversity", os.path.join(diversity_dir, "core_metrics_results/shannon_vector.qza"),  # Vector of alpha diversity values by sample
		"--o-visualization", os.path.join(statistical_analysis_dir, "shannon_diversity-significance.qzv"),  # Output stats
		"2>>", os.path.join(pipeline_reports, "5_phylogenetics_shannonSignificance-report.log")])  # Output stats report
		subprocess.run(shannonSignificance, shell=True)
		export(os.path.join(statistical_analysis_dir, "shannon_diversity-significance.qzv")) 

		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Creating Faith Phylogenetic Diversity (a measure of community richness) significance boxplot: in progress ..')
		faithSignificance =	' '.join([
		"qiime diversity alpha-group-significance",  # Calling qiime diversity alpha-group-significance function
		"--m-metadata-file", args.metadata,  # Metadata file
		"--quiet",  # Silence output if execution is successful
		"--i-alpha-diversity", os.path.join(diversity_dir, "core_metrics_results/faith_pd_vector.qza"),  # Vector of alpha diversity values by sample
		"--o-visualization", os.path.join(statistical_analysis_dir, "faith_pd-significance.qzv"),  # Output stats
		"2>>", os.path.join(pipeline_reports, "6_phylogenetics_faithSignificance-report.log")])  # Output report
		subprocess.run(faithSignificance, shell=True)
		export(os.path.join(statistical_analysis_dir, "faith_pd-significance.qzv")) 

		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Creating Evenness significance boxplot: in progress ..')
		evennessSignificance =	' '.join([
		"qiime diversity alpha-group-significance",  # Calling qiime diversity alpha-group-significance function
		"--m-metadata-file", args.metadata,  # Metadata file
		"--quiet",  # Silence output if execution is successful
		"--i-alpha-diversity", os.path.join(diversity_dir, "core_metrics_results/evenness_vector.qza"),  # Vector of alpha diversity values by sample
		"--o-visualization", os.path.join(statistical_analysis_dir, "evenness-significance.qzv"),  # Output stats
		"2>>", os.path.join(pipeline_reports, "7_phylogenetics_evennessSignificance-report.log")])  # Output report
		subprocess.run(evennessSignificance, shell=True)
		export(os.path.join(statistical_analysis_dir, "evenness-significance.qzv")) 
		return 

	def beta_diversity_analysis(self):
		""" Visually and statistically compare groups of beta diversity values """
		# Determine whether groups of samples are significantly different from one another using a permutation-based statistical test
		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Analysing sample composition in the context of Indication using PERMANOVA: in progress .')
		permanova_indication =	' '.join([
		"qiime diversity beta-group-significance",  # Calling qiime diversity beta-group-significance function
		"--m-metadata-file", args.metadata,  # Metadata file
		"--p-pairwise",  # Perform pairwise tests between all pairs of groups in addition to the test across all groups
		"--m-metadata-column", "Indication_I",  # Categorical sample metadata column
		"--quiet",  # Silence output if execution is successful
		"--i-distance-matrix", os.path.join(diversity_dir, "core_metrics_results/unweighted_unifrac_distance_matrix.qza"),  # Vector of beta diversity values by sample
		"--o-visualization", os.path.join(statistical_analysis_dir, "unweighted_unifrac_indication_I-significance.qzv"),  # Output stats
		# "2>>", os.path.join(pipeline_reports, "8_phylogenetics_permanova_indication-report.log")
		])  # Output report
		subprocess.run(permanova_indication, shell=True)
		export(os.path.join(statistical_analysis_dir, "unweighted_unifrac_indication_I-significance.qzv")) 

		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Analysing sample composition in the context of Indication_Code using PERMANOVA: in progress ..')
		permanova_indication_code =	' '.join([
		"qiime diversity beta-group-significance",  # Calling qiime diversity beta-group-significance function
		"--m-metadata-file", args.metadata,  # Metadata file
		"--p-pairwise",  # Perform pairwise tests between all pairs of groups in addition to the test across all groups
		"--m-metadata-column", "Indication_I",  # Categorical sample metadata column
		"--quiet",  # Silence output if execution is successful
		"--i-distance-matrix", os.path.join(diversity_dir, "core_metrics_results/weighted_unifrac_distance_matrix.qza"),  # Vector of beta diversity values by sample
		"--o-visualization", os.path.join(statistical_analysis_dir, "weighted_unifrac_indication_I-significance.qzv"),  # Output stats
		# "2>>", os.path.join(pipeline_reports, "9_phylogenetics_permanova_indication_code-report.log")
		])  # Output report
		subprocess.run(permanova_indication_code, shell=True)
		export(os.path.join(statistical_analysis_dir, "weighted_unifrac_indication_I-significance.qzv")) 

		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Adonis action to look at a multivariate mode: in progress ..')
		adonis =	' '.join([
		"qiime diversity adonis",
		"--i-distance-matrix", os.path.join(diversity_dir, "core_metrics_results/unweighted_unifrac_distance_matrix.qza"),
		"--m-metadata-file", args.metadata,
		"--o-visualization", os.path.join(adonis_dir, "unweighted_adonis.qzv"),
		"--p-formula", "Indication_I+Gender",
		# "2>>", os.path.join(pipeline_reports, "10_phylogenetics_adonis-report.log")
		])  # Output report
		subprocess.run(adonis, shell=True)
		export(os.path.join(adonis_dir, "unweighted_adonis.qzv")) 

		print(f'{datetime.now().strftime("%d.%m.%Y %H:%M")} Adonis action to look at a multivariate mode: in progress ..')
		adonis2 =	' '.join([
		"qiime diversity adonis",
		"--i-distance-matrix", os.path.join(diversity_dir, "core_metrics_results/weighted_unifrac_distance_matrix.qza"),
		"--m-metadata-file", args.metadata,
		"--o-visualization", os.path.join(adonis_dir, "weighted_adonis.qzv"),
		"--p-formula", "Indication_I+Gender",
		# "2>>", os.path.join(pipeline_reports, "11_phylogenetics_adonis2-report.log")
		])  # Output report
		subprocess.run(adonis2, shell=True)
		export(os.path.join(adonis_dir, "weighted_adonis.qzv")) 
		return

def differential_abundance():
	""" Following the tutorial for ANCOM, we want to collapse the table to genus level. 
	Before this, we will filter out low-abundance features. First we will removing those 
	that only appear 10 times or fewer across all samples, as well as those that only 
	appear in one sample """
	print("\n\t{0} DIFFERENTIAL ABUNDANCE".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	if not os.path.exists(abundanceAnalysis): os.makedirs(abundanceAnalysis)  # Creating the "qiime2_analysis" directory
	
	""" Applying basic filters to the features """
	print("{0} Applying basic filters to discard low frequent features: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	filteringFeatures = ' '.join([
	"qiime feature-table filter-features",  # Calling qiime2 feature-table filter-features function
	"--quiet",  # Silence output if execution is successful
	"--i-table", os.path.join(denoise_dir, "table.qza"),  # Input feature table to be summarized
	"--p-min-frequency", min_freq(),  # Least frequency that a feature must have to be retained
	"--p-min-samples 2",  # The minimum number of samples that a feature must be observed in to be retained
	"--o-filtered-table", os.path.join(abundanceAnalysis, "feature_table_filtered.qza"),  # Output file
	"2>>", os.path.join(pipeline_reports, "qiime2_filteringFeatures-report.log")])  # Output denoising report
	subprocess.run(filteringFeatures, shell=True)


	""" Collapse table at genus level """
	print("{0} Collapsing the table at genus level: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	collapsingTable = ' '.join([
	"qiime taxa collapse",  # Calling qiime taxa collapse function
	"--quiet",  # Silence output if execution is successful
	"--i-table", os.path.join(abundanceAnalysis, "feature_table_filtered.qza"),  # Input feature table to be summarized
	"--i-taxonomy", os.path.join(taxonomy_dir, "taxonomic_classification.qza"),
	"--p-level 6",  # The minimum number of samples that a feature must be observed in to be retained
	"--o-collapsed-table", os.path.join(abundanceAnalysis, "feature_table_for_ANCOM.qza"),  # Output file
	"2>>", os.path.join(pipeline_reports, "qiime2_collapsingTable-report.log")])  # Output denoising report
	subprocess.run(collapsingTable, shell=True)
	
	addPseudocounts = ' '.join([
	"qiime composition add-pseudocount",  # Calling qiime composition add-pseudocount function
	"--quiet",  # Silence output if execution is successful
	"--i-table", os.path.join(abundanceAnalysis, "feature_table_for_ANCOM.qza"),  # Input table
	"--o-composition-table", os.path.join(abundanceAnalysis, "ANCOM_ready_table.qza"),  # Output file
	"2>>", os.path.join(pipeline_reports, "qiime2_addPseudocounts-report.log")])  # Output denoising report
	subprocess.run(addPseudocounts, shell=True)

	# Run ANCOM
	print("{0} Running ANCOM: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	ancom = ' '.join([
	"qiime composition ancom",  # Calling qiime composition add-pseudocount function
	"--verbose",  # Silence output if execution is successful
	"--i-table", os.path.join(abundanceAnalysis, "ANCOM_ready_table.qza"),  # Input file
	"--m-metadata-file", args.metadata,
	"--m-metadata-column", "Indication_I",
	"--o-visualization", os.path.join(abundanceAnalysis, "ancom_group_results.qzv"),  # Input table
	"2>>", os.path.join(pipeline_reports, "qiime2_ancom-report.log")])  # Output denoising report
	subprocess.run(ancom, shell=True)
	export(os.path.join(abundanceAnalysis, "ancom_group_results.qzv"))
	return

def ml_approach():

	print("\n\t{0} MACHINE LEARNING APPROACHES".format(datetime.now().strftime("%d.%m.%Y %H:%M")))


	print("{0} Random Forst classifiers for predicting sample characteristics: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	random_forest = ' '.join([
	"qiime sample-classifier classify-samples",
	"--i-table", os.path.join(denoise_dir, "table.qza"),
	"--m-metadata-file", args.metadata,
	"--m-metadata-column", "Indication_I",
	# "--p-random-state", "500",
	# "--p-n-jobs", str(args.threads),
	"--output-dir", random_forest_dir,
	# "2>>", os.path.join(pipeline_reports, "qiime2_random_forest-report.log")
	])  # Output denoising report
	# subprocess.run(random_forest, shell=True)

	print("{0} Generating heatmap: in progress ..".format(datetime.now().strftime("%d.%m.%Y %H:%M")))
	heatmap = ' '.join([
	"qiime sample-classifier heatmap",
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

def min_freq():
	""" Based on the summary we will calculate a cut-off for how frequent a variant needs to be for it to be retained. 
	We will remove all ASVs that have a frequency of less than 0.0001% of the mean sample depth. This cut-off excludes """
	exportFile = os.path.join(denoise_dir, "feature_table.qzv")
	if not os.path.exists(exportFile[:-4]):
		subprocess.run("qiime tools export --input-path {0} --output-path {1}".format(exportFile, exportFile[:-4]), shell=True)

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
				mean_freq = mean_freq/samples

	# Calculation of the cut-off point
	cut_off = str(int(mean_freq * 0.0001))
	return cut_off

def export(exportFile):
	if os.path.isfile(exportFile):
		subprocess.run("qiime tools export --input-path {0} --output-path {1}".format(exportFile, exportFile[:-4]), shell=True)
	else:
		for path, subdir, folder in os.walk(exportFile):
			for files in folder:
				file = os.path.join(path, files)
				subprocess.run("qiime tools export --input-path {0} --output-path {1}".format(file, file[:-4]), shell=True)

	if exportFile.endswith(".qzv"):
		os.system('rm {0}'.format(exportFile))	
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
				os.system('mv {0} {1}'.format(file, qc_reports))
	
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
		# if batch == 'batch2':
		pairedReads, R1list = assess_inputdata(batch)  # Assessing all input fastq files and obtaining the pairs to be analysed
		print(f'- In total {len(R1list)} samples will be processed from {batch}..')
		
		# quality_control(batch, R1list)  # Checking the quality of the reads

		# ## Preprocessing of the input data
		# primer_removal(batch, pairedReads)  # Performing quality trimming and removal of all primers on both reads
		
		# ## Qiime2 analysis - denoising
		# denoising(batch)
	
	# merge_denoised_data(batch_to_analyse)

	# ## Downstream analysis
	# taxonomic_assignmnet()

	# phylogenetics(args.threads, args.metadata)

	# differential_abundance()

	# ml_approach()

	# summarisation()

if __name__ == "__main__": main()