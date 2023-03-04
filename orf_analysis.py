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

            for strand, seq in [(1, sequence), (-1, reverse_sequence)]:
                for frame in range(0, 3):
                    # Get the possible orfs of the frame
                    frame_orfs = find_all_orfs_in_current_frame(seq, frame, strand, min_orf_len, max_orf_len)

                    for orf in frame_orfs:
                        orf_counter += 1
                        # translate sequence to protein for a given ORF
                        protein_sequence = orf["seq"].translate(to_stop=True)
              
                        orf_record = SeqRecord(
                            seq=protein_sequence,
                            id= '%s_ORF%i' % (seq_record.id, orf_counter),
                            description='|Len[%i - %i]|Frame %i|Strand %s|' % (orf["start"] + 1, orf["end"] + 1, frame + 1, ("+" if strand == 1 else "-")),
                        )

                        # write to file
                        SeqIO.write(orf_record, orf_handle, 'fasta')
    

    
def find_all_orfs_in_current_frame(seq, frame, strand, min_orf_len, max_orf_len):
    """
    Finds all posible ORFs in the current reading frame by checking start and stop codons. 
    Start codon: ATG
    Stop codons: TAA, TAG y TGA
    Arguments: 
        nucleotides: string with ARNm nucleotides already transcribed
    Returns: array of substring of nucleotides between start and stop codons.
    """

    # Calculate where the reading frame starts and ends       
    start = frame
    end = math.floor((len(seq) - frame)/CODON_SIZE) * CODON_SIZE + frame
    nucleotides = seq[start:end]

    base_idx = 0 if strand == 1 else len(seq)-1
    
    orfs_list = []
    inside_orf = False
    
    for i in range(0, len(nucleotides[start:end]), CODON_SIZE):
        codon = str(nucleotides[i:i+CODON_SIZE])
        
        if not inside_orf and codon == START_CODON:
            inside_orf = True
            start_idx = i

        elif inside_orf and codon in END_CODONS:
            inside_orf = False
            end_idx = i + CODON_SIZE

            orf_len = end_idx - start_idx
            if min_orf_len <= orf_len and orf_len <= max_orf_len:
                orfs_list.append({
                    "seq": nucleotides[start_idx:end_idx],
                    "start": base_idx + strand * (start_idx + frame),
                    "end": base_idx + strand * ((end_idx - 1) + frame)
                })

    return orfs_list 

