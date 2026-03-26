# Name: Zareeb Lorenzana
# GitHub: Zareeb
# Repo: tbd
# Date: 03/24/2026

import re


class Sequence:
    def __init__(self, sequence: str = "", gene_id=None, species=None, *args, **kwargs):
        self.sequence_id = None
        self.kmers = []
        self.sequence = sequence
        self.gene_id = gene_id
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
        return f"{self.gene_id} {self.species}: {self.sequence}"
    
    def print_record(self):
        "Effectively the same as @property sequence"
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
        
        elif (self.species is None) and (self.gene_id is not None):
            return f">{self.gene_id}\n{self.sequence}"
        
        elif (self.gene_id is None) and (self.species is not None):
            return f">{self.species}\n{self.sequence}"
        
        elif (self.gene_id is not None) and (self.species is not None):
            return f">{self.gene_id} {self.species}\n{self.sequence}"
        
        else:
            return f">{self.sequence}"

class DNA(Sequence):
    def __init__(self, sequence: str = "", *args, **kwargs):
        super().__init__(sequence, *args, **kwargs)
        
    def validate_sequence(self, sequence: str) -> str:
        sequence = super().validate_sequence(sequence)
        sequence = re.sub("[^ATGCU>]", 'N', sequence)
        
        return sequence
    
    def reverse_complement(self) -> str:
        complement_dictionary = str.maketrans("ATGC", "TACG")
        complement = self.sequence.translate(complement_dictionary)
        reverse_complement = complement[::-1]

        return reverse_complement
    
    def print_info(self):
        return f"{self.gene_id} {self.species} {self.sequence}"
        
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
    query = \
    """
    >OX724082.1 Severe acute respiratory syndrome coronavirus 2 genome assembly, complete genome: monopartite
    AGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACG
    """
    
    query1 = "GATTACA"
    
    k = 3
    
    seq = Sequence(query1, "BRCA1", "H. Sapiens")
    
    print(seq.fasta())
    
if __name__ == '__main__':
    main()