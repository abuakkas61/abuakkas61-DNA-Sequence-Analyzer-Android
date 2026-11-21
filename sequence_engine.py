from collections import Counter
import math
import re

APP_TITLE = "DNA Sequence Analyzer"
VERSION = "V12"

# IUPAC nucleotide symbols, including ambiguity codes.
IUPAC_DNA = set("ACGTRYSWKMBDHVN")
IUPAC_RNA = set("ACGURYSWKMBDHVN")

DNA_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVN",
    "TGCAYRSWMKVHDBN"
)

RNA_COMPLEMENT = str.maketrans(
    "ACGURYSWKMBDHVN",
    "UGCAYRSWMKVHDBN"
)

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

AA_WEIGHTS = {
    "A":89.09, "R":174.20, "N":132.12, "D":133.10,
    "C":121.15, "Q":146.15, "E":147.13, "G":75.07,
    "H":155.16, "I":131.17, "L":131.17, "K":146.19,
    "M":149.21, "F":165.19, "P":115.13, "S":105.09,
    "T":119.12, "W":204.23, "Y":181.19, "V":117.15
}


# ================================================================
# SEQUENCE UTILITIES
# ================================================================

class SequenceEngine:
    @staticmethod
    def clean(text):
        """Remove FASTA headers, whitespace, and digits."""
        lines = []
        for line in text.splitlines():
            if not line.strip().startswith(">"):
                lines.append(line)
        return re.sub(r"[\s0-9]", "", "".join(lines)).upper()

    @staticmethod
    def validate(sequence):
        if not sequence:
            return False, "Please enter or load a sequence."

        symbols = set(sequence)
        allowed = IUPAC_DNA | IUPAC_RNA
        invalid = sorted(symbols - allowed)

        if invalid:
            return False, "Unsupported character(s): " + ", ".join(invalid)

        if "T" in sequence and "U" in sequence:
            return False, "The sequence contains both T and U. Use DNA or RNA, not a mixture."

        return True, "Sequence accepted."

    @staticmethod
    def is_rna(sequence):
        return "U" in sequence and "T" not in sequence

    @staticmethod
    def reverse(sequence):
        return sequence[::-1]

    @staticmethod
    def complement(sequence):
        if SequenceEngine.is_rna(sequence):
            return sequence.translate(RNA_COMPLEMENT)
        return sequence.translate(DNA_COMPLEMENT)

    @staticmethod
    def reverse_complement(sequence):
        return SequenceEngine.complement(sequence)[::-1]

    @staticmethod
    def transcribe(sequence):
        if SequenceEngine.is_rna(sequence):
            return sequence
        return sequence.replace("T", "U")

    @staticmethod
    def statistics(sequence):
        counts = Counter(sequence)
        length = len(sequence)
        a, t, g, c = (counts[x] for x in "ATGC")
        gc = g + c
        at = a + t
        canonical = a + t + g + c

        entropy = 0.0
        if canonical:
            for base in "ATGC":
                if counts[base]:
                    p = counts[base] / canonical
                    entropy -= p * math.log2(p)

        return {
            "length": length,
            "A": a, "T": t, "G": g, "C": c,
            "N": counts["N"],
            "ambiguous": length - canonical,
            "GC": (gc / canonical * 100) if canonical else 0.0,
            "AT": (at / canonical * 100) if canonical else 0.0,
            "purines": a + g,
            "pyrimidines": c + t,
            "gc_skew": ((g-c)/gc) if gc else 0.0,
            "at_skew": ((a-t)/at) if at else 0.0,
            "entropy": entropy,
            "canonical": canonical
        }

    @staticmethod
    def translate(sequence, frame=0):
        # Translation is defined here for DNA coding sequence.
        dna = sequence.replace("U", "T")
        protein = []
        for pos in range(frame, len(dna)-2, 3):
            codon = dna[pos:pos+3]
            protein.append(CODON_TABLE.get(codon, "X"))
        return "".join(protein)

    @staticmethod
    def six_frames(sequence):
        dna = sequence.replace("U", "T")
        reverse = SequenceEngine.reverse_complement(dna)
        result = {}
        for frame in range(3):
            result[f"+{frame+1}"] = SequenceEngine.translate(dna, frame)
        for frame in range(3):
            result[f"-{frame+1}"] = SequenceEngine.translate(reverse, frame)
        return result

    @staticmethod
    def amino_acid_counts(protein):
        return Counter(aa for aa in protein if aa in AA_WEIGHTS)

    @staticmethod
    def molecular_weight(protein):
        residues = [aa for aa in protein if aa in AA_WEIGHTS]
        if not residues:
            return 0.0
        return sum(AA_WEIGHTS[aa] for aa in residues) - 18.015 * (len(residues)-1)

    @staticmethod
    def codon_counts(sequence):
        dna = sequence.replace("U", "T")
        counts = Counter()
        for i in range(0, len(dna)-2, 3):
            codon = dna[i:i+3]
            if len(codon) == 3 and set(codon) <= set("ACGT"):
                counts[codon] += 1
        return counts

    @staticmethod
    def sliding_gc(sequence, window, step):
        x, y = [], []
        for start in range(0, len(sequence)-window+1, step):
            part = sequence[start:start+window]
            valid = sum(part.count(b) for b in "ATGC")
            gc = part.count("G") + part.count("C")
            if valid:
                x.append(start + window/2)
                y.append(gc / valid * 100)
        return x, y

    @staticmethod
    def sliding_skew(sequence, window, step):
        x, gc_values, at_values = [], [], []
        for start in range(0, len(sequence)-window+1, step):
            part = sequence[start:start+window]
            g, c, a, t = (part.count(b) for b in "GCAT")
            if g+c or a+t:
                x.append(start + window/2)
                gc_values.append((g-c)/(g+c) if g+c else 0)
                at_values.append((a-t)/(a+t) if a+t else 0)
        return x, gc_values, at_values

    @staticmethod
    def sliding_entropy(sequence, window, step):
        x, y = [], []
        for start in range(0, len(sequence)-window+1, step):
            part = sequence[start:start+window]
            valid = sum(part.count(b) for b in "ATGC")
            entropy = 0.0
            if valid:
                for base in "ATGC":
                    count = part.count(base)
                    if count:
                        p = count / valid
                        entropy -= p * math.log2(p)
                x.append(start + window/2)
                y.append(entropy)
        return x, y


class ORFEngine:
    @staticmethod
    def find_orfs(sequence, minimum_aa=10):
        dna = sequence.replace("U", "T")
        strands = {
            "+": dna,
            "-": SequenceEngine.reverse_complement(dna)
        }
        results = []

        for strand, seq in strands.items():
            for frame in range(3):
                for start in range(frame, len(seq)-2, 3):
                    if seq[start:start+3] != "ATG":
                        continue

                    for end in range(start+3, len(seq)-2, 3):
                        if seq[end:end+3] in ("TAA", "TAG", "TGA"):
                            aa_length = (end-start)//3
                            if aa_length >= minimum_aa:
                                if strand == "+":
                                    genomic_start = start + 1
                                    genomic_end = end + 3
                                else:
                                    genomic_start = len(seq) - (end + 3) + 1
                                    genomic_end = len(seq) - start

                                results.append({
                                    "strand": strand,
                                    "frame": frame + 1,
                                    "start": genomic_start,
                                    "end": genomic_end,
                                    "length": aa_length,
                                    "protein": SequenceEngine.translate(seq[start:end])
                                })
                            break
        return sorted(results, key=lambda item: (item["start"], item["strand"]))


