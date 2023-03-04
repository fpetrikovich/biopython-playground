import configparser

from error_helper import exit_with_error
from file_helper import file_exists
from constants import ANALYSIS_ORF, ARGUMENTS, CONFIG_DIR, FASTA_OUTPUT_FILE, FASTA_INPUT_FILE, MIN_ORF_LENGTH, MAX_ORF_LENGTH

def create_config(config_file, analysis):
    if not file_exists(CONFIG_DIR + config_file):
        exit_with_error("File %s does not exists in the %s directory." % (config_file, CONFIG_DIR))

    config = configparser.ConfigParser()
    config.read(CONFIG_DIR + config_file)
    analysis_config = config[analysis]

    return analysis_config

def handle_arguments(analysis, args):
    arguments = ARGUMENTS
    
    config_file = str(args.config_file)
    # Get the arguments from the config file
    config = create_config(config_file, analysis)
    for key in config: 
        value = config.get(key)
        # Convert the value too boolean if necessary
        if value == "yes" or value == "no": value = config.getboolean(key) 
        arguments[key] = value
    
    return arguments

def print_curr_arguments(analysis, arguments):
    print("[INFO] Using the following arguments:")

    if analysis == ANALYSIS_ORF:
        print('\t- %s : %s' % (FASTA_INPUT_FILE, arguments[FASTA_INPUT_FILE]))
        print('\t- %s : %s' % (FASTA_OUTPUT_FILE, arguments[FASTA_OUTPUT_FILE]))
        print('\t- %s : %s' % (MIN_ORF_LENGTH, arguments[MIN_ORF_LENGTH]))
        print('\t- %s : %s' % (MAX_ORF_LENGTH, arguments[MAX_ORF_LENGTH]))
    
    print() # new line
        
    