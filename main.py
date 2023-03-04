import argparse
import time

from arguments_helper import print_curr_arguments, handle_arguments
from file_helper import file_exists, valid_fasta_file, generate_output_path
from constants import ANALYSIS_ORF, CONFIG_FILE, FASTA_INPUT_FILE, FASTA_OUTPUT_FILE, ANALYSIS_TYPES, MIN_ORF_LENGTH, MAX_ORF_LENGTH
from error_helper import exit_with_error
from orf_analysis import analyze_orfs


def main():

    # Parse arguments
    parser = argparse.ArgumentParser(description="Bioinformatics Sequencing")

    # Add arguments
    parser.add_argument('-a', '--analysis', help='name of the analysis to run.',
                        type=str, default = ANALYSIS_ORF)
    parser.add_argument('-cf', '--config_file', help='name of the config file you want to use, configuration.ini by default.',
                        type=str, default = CONFIG_FILE)

    try:
        args = parser.parse_args()

        analysis = args.analysis

        if not analysis in ANALYSIS_TYPES: 
            exit_with_error("Possible analysis to run: {}".format(ANALYSIS_TYPES))

        args_to_use = handle_arguments(analysis, args)
        
        input_file = args_to_use[FASTA_INPUT_FILE]
        output_file = args_to_use[FASTA_OUTPUT_FILE]
        min_orf = args_to_use[MIN_ORF_LENGTH]
        max_orf = args_to_use[MAX_ORF_LENGTH]

    except Exception as e:
        exit_with_error("Invalid arguments. Check if the input values are correct.", e)

    # Param parsing and setup
    try:
        ###### ORF ANALYSIS ARGUMENT HANDLING ######

        if analysis == ANALYSIS_ORF:
            # Check if all arguments are present
            if input_file == None or output_file == None or min_orf == None or max_orf == None:
                exit_with_error("Missing arguments.")
            # Check if the input file is valid
            if not file_exists(input_file) or not valid_fasta_file(input_file):
                exit_with_error("Input file does not exists or is not a valid FASTA file.")
            # Check if filename is a valid fasta file
            if not valid_fasta_file(output_file):
                exit_with_error("Output file is not a valid FASTA file.")
            # Create the output file path
            output_file = generate_output_path(output_file)

    except Exception as e:
        exit_with_error("Error while parsing arguments.", e)


    # Run the exercise with the parsed params
    print("[INFO] Running analysis", analysis, "...")
    print_curr_arguments(analysis, args_to_use)
    start_time = time.time()

    try:
        if analysis == ANALYSIS_ORF:
            analyze_orfs(input_file, output_file, int(min_orf), int(max_orf))

    except Exception as e:
        exit_with_error("Error when running the analysis.", e)
    
    execution_time = int(time.time() - start_time)
    print("[DONE] Execution took %i mins %i seconds" % (int(execution_time / 60), execution_time % 60))


if __name__ == '__main__':
    main()
