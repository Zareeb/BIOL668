# Name: Zareeb Lorenzana
# GitHub: Zareeb
# Repo: Zareeb/BIOL668
# Date: 03/24/2026

import re


class Sequence:
    def __init__(self, sequence: str = "", gene_name=None, species=None, *args, **kwargs):
        self.sequence_id = None
        self.kmers = []
        self.sequence = sequence
        self.gene_name = gene_name
        self.species = species
    
    @property
    def sequence(self) -> str:
        return self._sequence

    @property
    def length(self) -> int:
        return len(self._sequence)
    
    @sequence.setter
    def sequence(self, sequence):
        if not isinstance(sequence, str):
            raise TypeError("Not a string input.")
        
        if sequence.strip().startswith(">"):            
            seq_id, parsed_seq = self.parse_fasta(sequence)
            self._sequence = self.validate_sequence(parsed_seq)
            
        else:
            self._sequence = self.validate_sequence(sequence)
                             
    @staticmethod
    def validate_sequence(sequence: str) -> str:
        
        sequence = sequence.upper()
        
        pattern = r"[A-Z>\-]+$"
        validated_sequence = re.sub(r"\s+", "", sequence)
        if not re.fullmatch(pattern, validated_sequence):
            raise ValueError(f"Unknown sequence format. Use FASTA with newline characters or raw sequence")
        
        return validated_sequence
    
    def __str__(self):
        
        return f"{self.gene_name} {self.species}: {self.sequence}"
    
    def print_record(self):
        # EFfectively similar to @property sequence getter
        
        return self.sequence
    
    def make_kmers(self, kmer_length: int = 3,):
        self.kmers = [self.sequence[i:i+kmer_length] for i in range(self.length - kmer_length + 1)]
        
        return self.kmers
    
    def parse_fasta(self, sequence):            
        re_match = re.findall(r">(.[^\n]*)\n([^>]*)", sequence, flags=re.MULTILINE)
        
        if len(re_match) > 1:
            raise ValueError("Multi-record FASTA not supported")
        
        for i in re_match:
            self.sequence_id = i[0]
            sequence = i[1]
        
        return self.sequence_id, sequence
            
    def fasta(self):
        # Returns a fasta formatted string
        if self.sequence_id is not None:
            return f">{self.sequence_id}\n{self.sequence}"
        
        elif (self.species is None) and (self.gene_name is not None):
            return f">{self.gene_name}\n{self.sequence}"
        
        elif (self.species is not None) and (self.gene_name is None):
            return f">{self.species}\n{self.sequence}"
        
        elif (self.gene_name is not None) and (self.species is not None):
            return f">{self.gene_name} {self.species}\n{self.sequence}"
        
        else:
            return f">{self.sequence}"

class DNA(Sequence):
    def __init__(self, sequence: str = "", gene_name = None, species = None, gene_id = None, *args, **kwargs):
        super().__init__(
            sequence,
            gene_name,
            species,
            *args,
            **kwargs)
        self.gene_id = gene_id
        
    def validate_sequence(self, sequence: str) -> str:
        sequence = super().validate_sequence(sequence)
        sequence = re.sub("[^ATGCU>]", 'N', sequence)
        
        return sequence
    
    def reverse_complement(self) -> str:
        complement_dictionary = str.maketrans("ATGC", "TACG")
        complement = self.sequence.translate(complement_dictionary)
        reverse_complement = complement[::-1]

        return reverse_complement
    
    def count_bases(self) -> int:
        counts_dictionary = {base: 0 for base in "ATGCU"}
        
        for base in self.sequence:
            if base in counts_dictionary.keys():
                counts_dictionary[base] += 1
                
        return counts_dictionary
    
    def analysis(self) -> float:
        if self.length == 0:
            return 0
        
        counts = self.count_bases()
        g_count = counts["G"]
        c_count = counts["C"]
        
        gc_content = (g_count + c_count)/self.length

        return f"{gc_content:.4f}"
    
    def six_frames(self):
        forward_sequence = self.sequence
        reverse_sequence = self.reverse_complement()
        self.six_frames = {}

        # Forward sequence
        for i in range(3):
            self.six_frames[f"F{i + 1}"] = forward_sequence[i:]

        # Reverse complement
        for i in range(3):
            self.six_frames[f"R{i + 1}"] = reverse_sequence[i:]

        return self.six_frames
    
    def print_info(self):
        headers = " ".join(item for item in (self.gene_id, self.species, self.gene_name) if item is not None)
        return f"{headers}: {self.sequence}"
        
class RNA(Sequence):
    def __init__(self, *args, **kwargs):
        pass


class Protein(Sequence):
    def __init__(self, *args, **kwargs):
        pass
    

class TableOfValues:
    standard_code = {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
        "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
        "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
        "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
        "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
        "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
        "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
        "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
        "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
        "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
        }
    
    kyte_doolittle = {
        'A':1.8,'C':2.5,'D':-3.5,'E':-3.5,'F':2.8,
        'G':-0.4,'H':-3.2,'I':4.5,'K':-3.9,'L':3.8, 
        'M':1.9,'N':-3.5,'P':-1.6,'Q':-3.5,'R':-4.5,
        'S':-0.8,'T':-0.7,'V':4.2,'W':-0.9,'X':0,'Y':-1.3
        }
    
    aa_mol_weights = {
        'A':89.09,'C':121.15,'D':133.1,'E':147.13,'F':165.19,
        'G':75.07,'H':155.16,'I':131.17,'K':146.19,'L':131.17,
        'M':149.21,'N':132.12,'P':115.13,'Q':146.15,'R':174.2,
        'S':105.09,'T':119.12,'V':117.15,'W':204.23,'X':0,'Y':181.19
        }
    
    @classmethod
    def get_standard_code(cls):
        return cls.standard_code
    
    @classmethod
    def get_kyte_doolittle(cls):
        return cls.kyte_doolittle
    
    @classmethod
    def get_aa_mol_weights(cls):
        return cls.aa_mol_weights
    

def main():
    query1 = \
    """
    >hg19_refGene_NM_007294 range=chr17:41196312-41277381 5'pad=0 3'pad=0 strand=- repeatMasking=none
    GCTGAGACTTCCTGGACGGGGGACAGGCTGTGGGGTTTCTCAGATAACTG
    GGCCCCTGCGCTCAGGAGGCCTTCACCCTCTGCTCTGGGTAAAGgtagta
    gagtcccgggaaagggacagggggcccaagtgatgctctggggtactggc
    """
    
    query2 = "GCTGAGACTTCCTG"
    
    gene_name = "BRCA1"
    species = "H. sapiens"
    gene_id = "AX5667.2"
    
    k = 3
    
    dna1 = DNA(query1)
    dna2 = DNA(query1, gene_name, species, gene_id)
    dna3 = DNA("GATCTC","my_dna","D.terebrans","AX5667.2")

    six_frames = dna1.six_frames()
    
    for frame, sequence in six_frames.items():
        print(f"{frame}: {sequence}")

    
if __name__ == '__main__':
    main()