"""
230127_writingWithConnor.py
Marcus Viscardi,    January 27, 2023

"""

from Bio import SeqIO

input_file = 'testRef.fasta'
output_file = 'testRef.new.fasta'

if __name__ == '__main__':
    # Read fasta file to fasta_sequences:
    fasta_sequences = SeqIO.parse(open(input_file), 'fasta')
    
    with open(output_file, 'w') as out_fasta:
        for fasta_obj in fasta_sequences:
            name = fasta_obj.id
            seq = fasta_obj.seq
            
            new_sequence = ''
            for nucleotide in seq:
                if nucleotide.lower() == 'c':
                    new_sequence += 't'
                else:
                    new_sequence += nucleotide
            out_fasta.write('>' + name)
            out_fasta.write('\n' + new_sequence)

