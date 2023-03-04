import os.path
from constants import FASTA_TYPES, FASTA_OUTPUT_DIR
from error_helper import exit_with_error

def check_file_is_not_dir(filename):
    length1 = len(filename.split('/'))
    length2 = len(filename.split('\\'))
    return length1 == 1 and length2 == 1

def save_file(file_name, content): 
	f = open(file_name, "w")
	f.write(content)
	f.close()

def file_exists(file):
    return os.path.exists(file)

def valid_fasta_file(file):
  return file.split('.')[-1] in FASTA_TYPES

def generate_output_path(filename):
    if not check_file_is_not_dir(filename):
        exit_with_error("Filename can not be a directory.")
    return FASTA_OUTPUT_DIR + filename    

def delete_file(path):
    try:
        os.remove(path)
    except Exception as e:
        exit_with_error("Error while deleting file.", e)
