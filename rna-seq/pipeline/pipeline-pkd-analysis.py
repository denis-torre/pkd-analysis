#################################################################
#################################################################
############### PKD Analysis Pipeline ################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Ma'ayan Laboratory,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
import sys, glob, os, docker
import pandas as pd
from rpy2.robjects import r, pandas2ri
pandas2ri.activate()

##### 2. Custom modules #####
# Pipeline running
sys.path.append('/Users/denis/Documents/Projects/scripts')
sys.path.append('pipeline/scripts')
import Support3 as S
import PkdAnalysis as P

#############################################
########## 2. General Setup
#############################################
##### 1. Variables #####
fastq_files = glob.glob('/Users/denis/Data/pkd-data/Raw/*.fastq.gz')
transcriptome_path = '/Users/denis/Data/alignment/s1-kallisto.dir/Mus_musculus.GRCm38.cdna.all.idx'
ensembl_annotation_file = '/Users/denis/Data/alignment/data.dir/mouse/GRCm38.p5-Gene_name-Ensembl_transcript_ID.txt'
annotation_file = '/Users/denis/Data/pkd-analysis/rna-seq/Raw/sample.txt'

##### 2. R Connection #####
r.source('pipeline/scripts/pkd-analysis.R')

#######################################################
#######################################################
########## S1. Kallisto
#######################################################
#######################################################

#############################################
########## 1. Alignment
#############################################

@follows(mkdir('s1-kallisto.dir'))

@transform(fastq_files,
		   regex(r'.*/(.*).fastq.gz'),
		   add_inputs(transcriptome_path),
		   r's1-kallisto.dir/\1')

def alignReads(infiles, outfile):

	# Get infiles
	fastq, transcriptome = infiles

	# Get container paths
	fastq_filename = os.path.basename(fastq)
	transcriptome_filename = os.path.basename(transcriptome)

	# Get Docker Client
	client = docker.from_env()

	# Get Volume
	container = client.containers.run(image = 'alignment',
									  volumes = {os.path.dirname(fastq): {'bind': '/input', 'mode': 'rw'}, os.path.join(os.getcwd(), outfile): {'bind': '/output', 'mode': 'rw'}, os.path.dirname(transcriptome): {'bind': '/transcriptome', 'mode': 'rw'}},
									  command = 'kallisto quant --index=/transcriptome/{transcriptome_filename} --output-dir=/output --single --fragment-length=200 --sd=20 /input/{fastq_filename}'.format(**locals()))

	print(container.decode('utf-8'))

#######################################################
#######################################################
########## S2. Expression Quantification
#######################################################
#######################################################

#############################################
########## 1. Merge Kallisto
#############################################

@follows(mkdir('s2-expression.dir'))

def mergeJobs():
	for value in ['rawcounts', 'tpm']:
		outdir = os.path.join('s2-expression.dir', value)
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		infiles = glob.glob('s1-kallisto.dir/*/abundance.tsv')
		outfile = '{outdir}/pkd-{value}_transcript.txt'.format(**locals())
		yield [infiles, outfile]

@files(mergeJobs)

def mergeAbundances(infiles, outfile):

	# Get value
	value = os.path.basename(outfile).split('_')[0].split('-')[1]

	# Read Dataframes
	dataframes = []
	for infile in infiles:
		kallisto_dataframe = pd.read_table(infile).rename(columns={'est_counts': 'rawcounts'})
		kallisto_dataframe['sample'] = infile.split('/')[-2]
		dataframes.append(kallisto_dataframe)

	# Cast
	transcript_expression_dataframe = pd.concat(dataframes).pivot_table(index='target_id', columns='sample', values=value)
	transcript_expression_dataframe.index.name = 'ensembl_transcript_id'

	# Round if counts
	if value == 'rawcounts':
		transcript_expression_dataframe = transcript_expression_dataframe.astype(int)

	# Write
	transcript_expression_dataframe.to_csv(outfile, sep='\t')

#############################################
########## 2. Get Gene Expression
#############################################

@follows(mkdir('s2-expression.dir'))

@transform(mergeAbundances,
		   suffix('_transcript.txt'),
		   add_inputs(ensembl_annotation_file),
		   '.txt')

def aggregateGeneExpression(infiles, outfile):

	# Get infiles
	transcript_rawcount_file, ensembl_annotation_file = infiles

	# Read data
	transcript_expression_dataframe = pd.read_table(transcript_rawcount_file, index_col='ensembl_transcript_id')
	ensembl_annotation_dataframe = pd.read_table(ensembl_annotation_file, index_col='Transcript stable ID')

	# Remove transcript versions
	transcript_expression_dataframe.index = [x.split('.')[0] for x in transcript_expression_dataframe.index]

	# Merge
	merged_dataframe = transcript_expression_dataframe.merge(ensembl_annotation_dataframe, left_index=True, right_index=True, how='left')

	# Get value
	value = outfile.split('-')[-1].split('.')[0]

	# Group by and aggregate
	if value == 'rawcounts':
		expression_dataframe = merged_dataframe.groupby('Gene name').sum()
	elif value == 'tpm':
		expression_dataframe = merged_dataframe.groupby('Gene name').mean()
	expression_dataframe.index.name = 'gene_symbol'

	# Write
	expression_dataframe.to_csv(outfile, sep='\t')

#######################################################
#######################################################
########## S3. Annotations
#######################################################
#######################################################

#############################################
########## 1. Sample Metadata
#############################################

@follows(mkdir('s3-metadata.dir'))

@files(annotation_file,
	   's3-metadata.dir/pkd-sample_metadata.txt')

def getSampleMetadata(infile, outfile):

	# Read dataframe
	sample_metadata_dataframe = pd.read_table(infile)
	print(sample_metadata_dataframe)

	# Set index
	sample_metadata_dataframe.index = [x.split('.')[0] for x in sample_metadata_dataframe['Read1 ']]
	sample_metadata_dataframe.index.name = 'sample_id'

	# Drop columns
	sample_metadata_dataframe.drop(['Sample', 'Genome', 'Experiment', 'Method', 'Read1 ', 'Read2', 'Splice_len', 'Win', 'Unique_align', 'alignment'], axis=1, inplace=True)

	# Write
	sample_metadata_dataframe.to_csv(outfile, sep='\t')

##################################################
##################################################
########## Run pipeline
##################################################
##################################################
pipeline_run([sys.argv[-1]], multiprocess=1, verbose=1)
print('Done!')