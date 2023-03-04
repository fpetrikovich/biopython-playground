import math
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

CODON_SIZE = 3
START_CODON = "ATG"
END_CODONS = ["TAA", "TAG", "TGA"]

def analyze_orfs(fasta_input_file, output_file, min_orf_len, max_orf_len):
    
    with open(output_file, 'w') as orf_handle:
        for seq_record in SeqIO.parse(fasta_input_file, "fasta"):
            print('Working with Fasta sequence record %s' % seq_record.id)

            orf_counter = 0

            sequence = seq_record.seq
            reverse_sequence = sequence.reverse_complement()

            for strand, seq in [("+", sequence), ("-", reverse_sequence)]:
                for frame in range(0, 3):
                     # Calculate where the reading frame starts and ends       
                    start = frame
                    end = math.floor((len(seq) - frame)/CODON_SIZE) * CODON_SIZE + frame
                    # Get the possible orfs of the frame
                    frame_orfs = find_all_orfs_in_current_frame(seq[start:end], min_orf_len, max_orf_len)

                    for orf in frame_orfs:
                        orf_counter += 1
                        # translate sequence to protein for a given ORF
                        protein_sequence = orf["seq"].translate(to_stop=True)
              
                        orf_record = SeqRecord(
                            seq=protein_sequence,
                            id= '%s_ORF_%i' % (seq_record.id, orf_counter),
                            description='|Len[%i - %i]|Frame %i|Strand %s|' % (orf["start"], orf["end"], frame, strand),
                        )

                        # write to file
                        SeqIO.write(orf_record, orf_handle, 'fasta')
    

    
def find_all_orfs_in_current_frame(nucleotides, min_orf_len, max_orf_len):
    """
    Finds all posible ORFs in the current reading frame by checking start and stop codons. 
    Start codon: ATG
    Stop codons: TAA, TAG y TGA
    Arguments: 
        nucleotides: string with ARNm nucleotides already transcribed
    Returns: array of substring of nucleotides between start and stop codons.
    """

    start_indexes = []
    orfs_list = []
    inside_orf = False
    
    for i in range(0, len(nucleotides), CODON_SIZE):
        codon = str(nucleotides[i:i+CODON_SIZE])
        
        if codon == START_CODON:
            inside_orf = True
            start_indexes.append(i)

        elif inside_orf and codon in END_CODONS:
            inside_orf = False
            for start_idx in start_indexes:
                orf_len = i - start_idx
                if min_orf_len <= orf_len and orf_len <= max_orf_len:
                    orfs_list.append({
                        "seq": nucleotides[start_idx:i],
                        "start": start_idx,
                        "end": i
                    })
            # Clear the list of start indexes
            start_indexes = []

    return orfs_list 

