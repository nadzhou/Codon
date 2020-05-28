from Bio import SeqIO



def seq_extract(file_name): 
    return SeqIO.read(file_name, "fasta")
    