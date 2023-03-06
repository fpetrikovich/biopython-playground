import math
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

CODON_SIZE = 3
START_CODON = "ATG"
END_CODONS = ["TAA", "TAG", "TGA"]

def analyze_orfs(fasta_input_file, output_file, min_orf_len, max_orf_len):

    orf_records = []
    
    # Open the input and output to create the handles
    with open(fasta_input_file, 'r') as input_file, open(output_file, 'w') as orf_handle:
        for seq_record in SeqIO.parse(input_file, "fasta"):
            print('Working with Fasta sequence record %s' % seq_record.id)

            orf_counter = 0

            sequence = seq_record.seq
            reverse_sequence = sequence.reverse_complement()

            for strand, seq in [(1, sequence), (-1, reverse_sequence)]:
                for frame in range(0, 3):
                    # Get the possible orfs of the frame
                    frame_orfs = find_all_orfs_in_current_frame(seq, frame, strand, min_orf_len, max_orf_len)

                    for orf in frame_orfs:
                        orf_counter += 1
                        # translate sequence to protein for a given ORF
                        protein_sequence = orf["seq"].translate(to_stop=True)
              
                        orf_records.append(SeqRecord(
                            seq=protein_sequence,
                            id= '%s_ORF%i' % (seq_record.id, orf_counter),
                            description='|Frame %i|Strand %s|Pos[%i - %i]|Nt Len %i|AA Len %i|' % (frame + 1, ("+" if strand == 1 else "-"), orf["start"] + 1, orf["end"] + 1, orf["len"], len(protein_sequence)),
                        ))

        # write to file
        SeqIO.write(orf_records, orf_handle, 'fasta')
    

    
def find_all_orfs_in_current_frame(seq, frame, strand, min_orf_len, max_orf_len):
    """
    Finds all possible ORFs in the current reading frame by checking start and stop codons. 
    Arguments: 
        seq: the complete nucleotide sequence 
        frame: the frame we are analyzing, either 0, 1 or 2
        strand: identification of which DNA strand we are analyzing, positive (1) or negative (-1)
        min_orf_len: minimum length an ORF must have to be considered.
        max_orf_len: maximum length an ORF must have to be considered.
    Returns: list of subsequences of nucleotides between start and stop codons (ORFs).
    """

    # Calculate where the reading frame starts and ends. Must be a multiple of 3.        
    start = frame 
    end = math.floor((len(seq) - frame)/CODON_SIZE) * CODON_SIZE + frame
    nucleotides = seq[start:end]

    base_idx = 0 if strand == 1 else len(seq)-1
    
    # Where the ORFs sequences will be stored
    orfs_list = []

    # Boolean indicating if a start codon was found => inside an ORF
    inside_orf = False
    
    # Loop through the sequence in steps of 3 (codon size)
    for i in range(0, len(nucleotides), CODON_SIZE):
        codon = str(nucleotides[i:i+CODON_SIZE])
        
        # Start of an ORF
        if not inside_orf and codon == START_CODON:
            inside_orf = True
            orf_start_idx = i

        # End of an ORF
        elif inside_orf and codon in END_CODONS:
            inside_orf = False
            orf_end_idx = i + CODON_SIZE
            orf_len = orf_end_idx - orf_start_idx

            # Check if the ORF length is inside the limits
            if min_orf_len <= orf_len and orf_len <= max_orf_len:
                # Save the ORF to the list
                orfs_list.append({
                    "seq": nucleotides[orf_start_idx:orf_end_idx],
                    "start": base_idx + strand * (orf_start_idx + frame),
                    "end": base_idx + strand * ((orf_end_idx - 1) + frame),
                    "len": orf_len
                })

    return orfs_list 

