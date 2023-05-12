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


def main():
    seq = next(SeqIO.parse('covid_sequence.fasta', "fasta")).seq

    print("Nucleotides frequency: ", get_nucleotides_frequency(seq))
    print("GC: ", gc_fraction(seq))

    rna_seq = get_rna_seq(seq)
    print("RNA seq: ", rna_seq)
    print("Protein seq: ", get_protein_seq(rna_seq))


if __name__ == '__main__':
    main()
