
## CONFIGURATION + ARGUMENTS CONSTANTS ##

CONFIG_FILE = "configuration.ini"

ANALYSIS_ORF = "orf"
ANALYSIS_TYPES = [ANALYSIS_ORF]

FASTA_INPUT_FILE = "fasta_input_file"
FASTA_OUTPUT_FILE = "fasta_output_file"
MIN_ORF_LENGTH = "min_orf_len"
MAX_ORF_LENGTH = "max_orf_len"

ARGUMENTS = {
    FASTA_INPUT_FILE: None,
    FASTA_OUTPUT_FILE: None,
    MIN_ORF_LENGTH: None,
    MAX_ORF_LENGTH: None,
}

## DIRECTORY CONSTANTS ##

FASTA_INPUT_DIR = "fasta_inputs/"
FASTA_OUTPUT_DIR = "fasta_outputs/"
CONFIG_DIR = "configuration_files/"

## FILE HANDLING CONSTANTS ##

FASTA_EXTENSION = ".fasta"
FASTA_TYPES = ['fas', 'fasta', 'faa']


