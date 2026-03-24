# Name: Zareeb Lorenzana
# GitHub: Zareeb
# Repo: tbd
# Date: 03/24/2026

import re


class Sequence:
    def __init__(self, sequence: str = "", kmer_length: int = 3, *args, **kwargs):
        self.sequence = sequence
        self.kmer_length = kmer_length
        self.kmers = [] # Intializes to an empty list
    
    @property
    def sequence(self) -> str:
        return self._sequence

    @property
    def length(self) -> int:
        return len(self._sequence)
        
    @sequence.setter
    def sequence(self, sequence):
        self._sequence = self.validate_sequence(sequence)
        
    @staticmethod
    def validate_sequence(sequence: str) -> str:
        sequence = sequence.upper()
        
        pattern = r"^[A-Z]+$"
        validated_sequence = re.sub(r"\s+", "", sequence)

        if not re.fullmatch(pattern, validated_sequence):
            raise ValueError(f"You can only use standard DNA, RNA, or Protein sequences")

        return validated_sequence
    
    def make_kmers(self):
        # Makes overlapping kmers of a given length
        self.kmers = [self.sequence[i:i+self.kmer_length] for i in range(self.length - self.kmer_length + 1)]
        
        return self.kmers
    
    def fasta(self):
        # Returns a fasta formatted string
        pass

    

class DNA(Sequence):
    def __init__(self, *args, **kwargs):
        pass

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
    query = "To be or not to be that is the question"
    seq = Sequence(query)
    print(seq.sequence)
    print(seq.make_kmers())

if __name__ == '__main__':
    main()