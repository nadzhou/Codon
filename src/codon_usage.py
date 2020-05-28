from collections import defaultdict
from seq_extract import seq_extract


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from Bio import SeqIO


genetic_table = { 
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',}

biasness_dict = {
        'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0, 'CTT': 0,
        'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0, 'ATC': 0,
        'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0, 'GTA': 0,
        'GTG': 0, 'TAT': 0, 'TAC': 0, 'TAA': 0, 'TAG': 0,
        'CAT': 0, 'CAC': 0, 'CAA': 0, 'CAG': 0, 'AAT': 0,
        'AAC': 0, 'AAA': 0, 'AAG': 0, 'GAT': 0, 'GAC': 0,
        'GAA': 0, 'GAG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
        'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
        'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
        'GCC': 0, 'GCA': 0, 'GCG': 0, 'TGT': 0, 'TGC': 0,
        'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0, 'CGA': 0,
        'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0, 'AGG': 0,
        'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}


class CodonUsage: 
    """ Determine the codon bias of the given DNA sequence.
    """

    def __init__(self, seq, genetic_table, biasness_dict): 
        self.seq = seq
        self.genetic_table = genetic_table
        self.biasness_dict = biasness_dict

        self.rscu = {}
        self.aa_codon_count = []
        self.codon_sums_dict = {}
        self.aa_sum_over_codons = {}
        self.codon_bias = None


    def calculate_rscu(self):
        """ Calculate the RSCU values given the codon and amino acid counts

            Returns: 
                rscu [dict]: Dictionary of codon with RWSCU values
        """
        for codon, aa in self.genetic_table.items():

            rscu_value = 0
            
            if aa in self.aa_sum_over_codons and codon in biasness_dict: 
                rscu_value = (1/self.codon_sums_dict[aa]) * self.aa_sum_over_codons[aa]
                rscu_value = self.biasness_dict[codon] / (rscu_value + 1) 

            self.rscu[codon] = rscu_value
            
        return self.rscu


    def sum_over_codons(self): 
        """ See which codons code for the same amino acid, then sum over 
            the counts to get an overall picture of codon bias. 
        """
        self.codon_count_dict = defaultdict(list)
        for x, y in self.genetic_table.items(): 
            if x in self.biasness_dict: 
                self.codon_count_dict[y].append(self.biasness_dict[x])


        for x, y in self.codon_count_dict.items(): 
            count = 0
            for i in y: 
                count += i
            self.aa_sum_over_codons[x] = count
        

    def make_aa_codon_count(self): 
        """ Given the codons, find the amino acids they code for and transfer these 
            amino acids with the counts. 
        """
        for i, j in self.genetic_table.items(): 
            if i in self.biasness_dict: 
                new_value = (self.genetic_table[i], 
                            self.aa_sum_over_codons[j], 
                            i, 
                            self.biasness_dict[i])
                
                self.aa_codon_count.append(new_value)
    
    
        for i, j in self.codon_count_dict.items(): 
            count = 0 
            for y in j: 
                count += 1
            self.codon_sums_dict[i] = count


    def count_codons(self): 
        """ Count occurence of each codon given the genetic code table
        """
        for i in range(0, len(self.seq), 3): 
            codon = self.seq[i : i + 3]
            if codon in self.biasness_dict: 
                self.biasness_dict[codon] += 1
           

    def call(self):
        """ Run the functions in this class
        """
        if self.codon_bias == None: 
            self.count_codons()
            self.sum_over_codons()
            self.make_aa_codon_count()
            self.codon_bias = self.calculate_rscu()

        return self.codon_bias

def main(): 

    index = np.arange(64)
    bar_width = 0.35

    rnase_record = seq_extract("rnase.fasta")
    rnase2_record = seq_extract("rnase2.fasta")

    coli_job = CodonUsage(rnase_record.seq, genetic_table, biasness_dict)
    result_1 = coli_job.call()

    enteceros_job = CodonUsage(rnase2_record.seq, genetic_table, biasness_dict)
    result_2 = enteceros_job.call()

    fig, ax = plt.subplots()

    coli = ax.bar(index, list(result_1.values()), bar_width, label="Coli")
    enteceros = ax.bar(index+bar_width, list(result_2.values()), bar_width, label="Enetecros")
    ax.set_title("Codon Usage for RNAse gene between two species")
    ax.set_xlabel("Codon")
    ax.set_ylabel("RSCU values")
    ax.xaxis.labelpad = 20
    ax.set_xticklabels(list(result_1.keys()))
    ax.legend()

    plt.show()

if __name__  == '__main__': 
    main()