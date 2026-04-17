# Name: Zareeb Lorenzana
# GitHub: Zareeb
# Repo: Zareeb/BIOL668
# Date: 03/26/2026

import re

class Sequence:
    def __init__(self, sequence: str = "", 
                 gene_name = None, 
                 species = None, 
                 *args, 
                 **kwargs
                 ):
        """Initializes a Sequence object.

        Assigns the input sequence to the internal sequence attribute using the property setter,
        which automatically validates and normalizes the sequence. Optionally stores metadata
        such as gene name and species.

        Args:
            sequence (str, optional): Input sequence string. Can be a raw sequence or FASTA-formatted string.
            gene_name (str, optional): Name of the gene associated with the sequence.
            species (str, optional): Species associated with the sequence.
        """
        self.sequence_id = None
        self.kmers = []
        self.sequence = sequence
        self.gene_name = gene_name
        self.species = species
    
    @property
    def sequence(self) -> str:
        """Returns the validated sequence.

        Accesses the internally stored sequence string after validation and normalization.

        Returns:
        str: The validated sequence.
        """
        return self._sequence

    @property
    def length(self) -> int:
        """Returns the length of the sequence.

        Calculates the number of characters in the validated sequence.

        Returns:
            int: Length of the sequence.
        """
        
        return len(self._sequence)
    
    @sequence.setter
    def sequence(self, sequence: str):
        """Sets and validates the sequence.

        Determines whether the input is a FASTA-formatted string or a raw sequence. If FASTA,
        parses the header and sequence, storing the sequence ID and validating the parsed sequence.
        Otherwise, directly validates the input sequence.

        Args:
            sequence (str): Input sequence string.

        Raises:
            TypeError: If the input is not a string.
        """
        if not isinstance(sequence, str):
            raise TypeError("Not a string input.")
        
        if sequence.strip().startswith(">"):            
            seq_id, parsed_seq = self.parse_fasta(sequence)
            self._sequence = self.validate_sequence(parsed_seq)
            
        else:
            self._sequence = self.validate_sequence(sequence)
                             
    @staticmethod
    def validate_sequence(sequence: str) -> str:
        """Validates and normalizes a sequence string.

        Converts the sequence to uppercase, removes all whitespace characters, and ensures
        that the sequence contains only valid characters (A–Z, '*', '>', or '-').

        Args:
            sequence (str): Input sequence string.

        Returns:
            str: Cleaned and validated sequence.

        Raises:
            ValueError: If the sequence contains invalid characters.
        """
        
        sequence = sequence.upper()
        
        pattern = r"[A-Z>?*\-]+$"
        validated_sequence = re.sub(r"\s+", "", sequence)
        if not re.fullmatch(pattern, validated_sequence):
            raise ValueError(f"Unknown sequence format. Use FASTA with newline characters or raw sequence")
        
        return validated_sequence
    
    def __str__(self) -> str:
        """Returns a string representation of the sequence.

        Combines available metadata (gene name and species) with the sequence.

        Returns:
            str: Formatted string representation.
        """
        headers = " ".join(item for item in (self.gene_name, self.species) if item is not None)
        
        return f"{headers}: {self.sequence}"
    
    def print_record(self) -> str:
        """Returns the sequence string.

        Provides access to the sequence similar to the sequence property getter.

        Returns:
            str: The sequence.
        """
        
        # EFfectively similar to @property sequence getter
        return self.sequence 
    
    def make_kmers(self, kmer_length: int = 3,) -> list:
        """Generates k-mers from the sequence.

        Splits the sequence into overlapping substrings of specified length.

        Args:
            kmer_length (int, optional): Length of each k-mer. Defaults to 3.

        Returns:
            list: List of k-mer substrings.
        """
        self.kmers = [self.sequence[i:i+kmer_length] for i in range(self.length - kmer_length + 1)]
        
        return self.kmers
    
    def parse_fasta(self, sequence: str) -> tuple:
        """Parses a FASTA-formatted string.

        Extracts the sequence ID from the header line and retrieves the associated sequence.
        Only single-record FASTA inputs are supported.

        Args:
            sequence (str): FASTA-formatted string.

        Returns:
            tuple: (sequence_id, sequence)

        Raises:
            ValueError: If multiple FASTA records are detected.
        """          
        re_match = re.findall(r">(.[^\n]*)\n([^>]*)", sequence, flags=re.MULTILINE)
        
        if len(re_match) > 1:
            raise ValueError("Multi-record FASTA not supported")
        
        for i in re_match:
            self.sequence_id = i[0]
            sequence = i[1]
        
        return self.sequence_id, sequence
            
    def fasta(self):
        """Returns the sequence in FASTA format.

        Constructs a FASTA-formatted string using available metadata in priority order:
        sequence ID, gene name, species, or the sequence itself.

        Returns:
            str: FASTA-formatted string.
        """
        
        # Returns a fasta formatted string
        if self.sequence_id is not None:
            header = self.sequence_id
        
        elif self.gene_name  and self.species:
            header = f"{self.gene_name} {self.species}"
            
        elif self.gene_name:
            header = self.gene_name
        
        elif self.species:
            header = self.species
            
        else:
            header = "sequence"

        return f">{header}\n{self.sequence}"

class DNA(Sequence):
    def __init__(self, sequence: str = "", 
                 gene_name = None, 
                 species = None, 
                 gene_id = None, 
                 *args, 
                 **kwargs):
        super().__init__(
            sequence,
            gene_name,
            species,
            *args,
            **kwargs)
        self.gene_id = gene_id
        self.frames = {}
        
    def validate_sequence(self, sequence: str) -> str:
        sequence = super().validate_sequence(sequence)
        validated_sequence = re.sub("[^ATGCU>]", 'N', sequence)
        
        return validated_sequence
    
    def reverse_complement(self) -> str:
        if "U" not in self.sequence:
            complement_dictionary = str.maketrans("ATGC", "TACG")
        
        else:          
            complement_dictionary = str.maketrans("AUGC", "UACG")
        
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

        # Forward sequence
        for i in range(3):
            self.frames[f"F{i + 1}"] = forward_sequence[i:]

        # Reverse complement
        for i in range(3):
            self.frames[f"R{i + 1}"] = reverse_sequence[i:]

        return self.frames
    
    def print_info(self):
        headers = " ".join(item for item in (self.gene_id, self.species, self.gene_name) if item is not None)
        return f"{headers}: {self.sequence}"
        
        
class RNA(DNA):
    def __init__(self, sequence = "", *args, **kwargs):
        super().__init__(
            sequence.replace("T", "U"),
            *args,
            **kwargs
            )
        self.codons = []
        
    def validate_sequence(self, sequence):
        sequence = super().validate_sequence(sequence)
        sequence = sequence.replace("T", "U")

        return sequence
        
    def make_codons(self, sequence: str = None):
        if sequence is None:
            sequence = self.sequence
            
        size = 3
        codons = []
        for i in range(0, len(sequence), size):
            if len(sequence[i:i + size]) % size == 0:
                codons.append(sequence[i:i + size])
                
        self.codons = codons
        
        return codons
    
    def translate(self, frames = None):
        standard_table = TableOfValues().get_standard_code()
        protein = ""
        stop_at_stop = False
        
        if frames is None:
            codons = self.make_codons()
            
        else:
            codons = frames
            
        for codon in codons:
            try:
                aa = standard_table[codon]

            except KeyError:
                aa = "?"

            if aa != "*":
                protein += aa

            else:
                protein += aa

                if stop_at_stop:
                    return protein

        return protein
    
    def find_orf(self):
        # Assigns sequences
        forward_sequence = self.sequence
        reverse_sequence = self.reverse_complement()

        # Stores codons for each reading frame in a dictionary
        frames = {}
        translated_frames = {}

        for i in range(3):
            frames[f"F{i + 1}"] = self.make_codons(sequence=forward_sequence[i:])
            translated_frames[f"F{i + 1}"] = self.translate(frames=frames[f"F{i + 1}"])
            
        for i in range(3):
            frames[f"R{i + 1}"] = self.make_codons(sequence=reverse_sequence[i:])
            translated_frames[f"R{i + 1}"] = self.translate(frames=frames[f"R{i + 1}"])

        open_reading_frames = {}
        
        for frame in translated_frames:
            current_sequence = "".join(translated_frames[frame])
            match = re.search(r"M[A-Z?]+\*", current_sequence)

            try:
                open_reading_frames[frame] = match.group()

            except AttributeError:
                open_reading_frames[frame] = ""

        return open_reading_frames


class Protein(Sequence):
    def __init__(self, sequence: str = "", *args, **kwargs):
        super().__init__(
            sequence,
            *args,
            **kwargs
            )
        
    def validate_sequence(self, sequence: str) -> str:
        sequence = super().validate_sequence(sequence)
        validated_sequence = re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "X", sequence)
        
        return validated_sequence
    
    def total_hydro(self, window: int = None):
        hscale = TableOfValues.get_kyte_doolittle()
        protein_sequence = self.sequence
        
        hscores_dict = {}
        
        if window is not None:
            
            j = window 
            
            for i in range(0, len(protein_sequence) - window + 1):
                window_score = 0 # Resets scores
                
                for aa in protein_sequence[i:j]:
                    
                    if aa in hscale.keys():
                        window_score += hscale[aa]
                    else:
                        window_score += 0
                
                window_sequence = protein_sequence[i:j]
                j += 1
                average_hscores = round(window_score/window, 2)
                
                hscores_dict[i+1] = (window_sequence, average_hscores)
            
            return hscores_dict
        
        else:
            total_hscore = 0
            for aa in protein_sequence:
                total_hscore += hscale[aa]
                
            return total_hscore

    def mol_weight(self):
        aa_mol_weights = TableOfValues.get_aa_mol_weights()
        molar_mass  = 0
        
        for aa in self.sequence:
            molar_mass += aa_mol_weights[aa]
        
        return molar_mass


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

    # six_frames = dna1.six_frames()
    
    # for frame, sequence in six_frames.items():
    #     print(f"{frame}: {sequence}")
        
    rna = RNA("CCCTATGAACCCATTGCTTAGGAGAACCGTATGATCCCTAGCTCATTAAGCGTAGAGTGAGGGTTCGAATGTGGAACTGATGCTTATCATTCCTCATCT")
    
    six_frames = rna.six_frames()
    
    # print(rna.sequence)
    # print(rna.reverse_complement())
    
    for frame, sequence in six_frames.items():
        print(f"{frame}: {RNA(sequence).translate()}")
    
    translated_protein = rna.translate()
    print(translated_protein)
    
    orf = rna.find_orf()
    for frame, sequence in orf.items():
        print(f"{frame}: {sequence}")
    
    prot = Protein("PMNPLLRRTV*SLAH*A*SEGSNVELMLIIPH")
    prot1 = DNA(query1)
    
    hydropathic_values = prot.total_hydro()
    hydropathic_values_window = prot.total_hydro(window=5)

    print("\nTotal hydrophobicity score = ", hydropathic_values)
    
    print("\nSliding window average hydrophobicity score:")
    for val in hydropathic_values_window.values():
        print(f"{val[0]}: {val[1]}")
        
    print(f"\nTotal molecular weight = {prot.mol_weight():.2f} g/mol")
    r=RNA(" AUGg?ATATUUTAAGGACctttaGGATCCACUAG ","my_rna","G.gallus","R5990999")
    print(r)
        
if __name__ == '__main__':
    main()