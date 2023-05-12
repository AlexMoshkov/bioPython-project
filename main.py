from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction


def get_nucleotides_frequency(seq: Seq) -> dict[str, int]:
    total_count = len(seq)
    return {
        'A': seq.count('A') / total_count,
        'G': seq.count('G') / total_count,
        'C': seq.count('C') / total_count,
        'T': seq.count('T') / total_count,
    }


def get_rna_seq(dna_seq: Seq) -> Seq:
    return dna_seq.transcribe()


def get_protein_seq(rna_seq: Seq) -> Seq:
    return rna_seq.translate()


def get_codon_frequency(dna_seq: Seq) -> dict[str, int]:
    codon_freq = dict()
    for i in range(0, len(dna_seq), 3):
        codon = dna_seq[i:i+3]
        if codon not in codon_freq:
            codon_freq[codon] = 0
        codon_freq[codon] += 1

    total_count = len(dna_seq) // 3
    for codon, count in codon_freq.items():
        codon_freq[codon] = count / total_count
    return codon_freq


def main():
    seq = next(SeqIO.parse('covid_sequence.fasta', "fasta")).seq

    print("Nucleotides frequency: ", get_nucleotides_frequency(seq))
    print("GC: ", gc_fraction(seq))

    rna_seq = get_rna_seq(seq)
    print("RNA seq: ", rna_seq)
    protein_seq = get_protein_seq(rna_seq)
    print("Protein seq: ", protein_seq)

    print("Codon freq: \n", get_codon_frequency(seq))


if __name__ == '__main__':
    main()
