from Bio import SeqIO

def ler_fasta(arquivo):
    sequencias = {}
    
    try:
        for record in SeqIO.parse(arquivo, "fasta"):
            header = ">" + record.description
            sequencias[header] = str(record.seq)
        
        print(f"Leitura concluída: {len(sequencias)} sequências encontradas no arquivo {arquivo}")
    except Exception as e:
        print(f"Erro ao ler o arquivo FASTA: {e}")
        
    return sequencias

def dna_para_rna(sequencia_dna):
    from Bio.Seq import Seq
    
    seq_dna = Seq(sequencia_dna.upper())
    seq_rna = seq_dna.transcribe()
    
    return str(seq_rna)

def criar_tabela_genetica():
    tabela = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    return tabela

def traduzir_rna(sequencia_rna, tabela, frame=0):
    proteina = ""
    for i in range(frame, len(sequencia_rna) - 2, 3):
        codon = sequencia_rna[i:i+3]
        if len(codon) == 3:
            aminoacido = tabela.get(codon, "X")
            proteina += aminoacido
    return proteina


def complemento_reverso(sequencia):
    from Bio.Seq import Seq
    
    seq_dna = Seq(sequencia.upper())
    seq_complementar = seq_dna.reverse_complement()
    
    return str(seq_complementar)

def salvar_fasta(sequencias, arquivo):
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    
    records = []
    for header, seq in sequencias.items():
        if header.startswith('>'):
            header = header[1:]
        
        record = SeqRecord(
            Seq(seq),
            id="",
            description=header
        )
        records.append(record)
    
    try:
        SeqIO.write(records, arquivo, "fasta")
        return True
    except Exception as e:
        print(f"Erro ao salvar o arquivo FASTA: {e}")
        return False

def main():
    try:
        caminho_fasta = "c:\\Users\\gabri\\Documents\\Unifesp\\Atividades\\AlgBio\\Atividade 4\\Ecoli_Sakai_cds_from_genomic.fna"
        sequencias_dna = ler_fasta(caminho_fasta)
        
        print("Transcrevendo DNA para RNA...")
        sequencias_rna = {}
        for header, seq in sequencias_dna.items():
            novo_header = header + " RNA"
            sequencias_rna[novo_header] = dna_para_rna(seq)
        
        print("Salvando arquivo de RNA...")
        caminho_rna = "c:\\Users\\gabri\\Documents\\Unifesp\\Atividades\\AlgBio\\Atividade 4\\Sakai_RNA.fasta"
        salvar_fasta(sequencias_rna, caminho_rna)
        print(f"Arquivo salvo: {caminho_rna}")
        
        tabela_genetica = criar_tabela_genetica()
        
        for frame in range(1, 7):
            print(f"Traduzindo sequências para o Frame {frame}...")
            sequencias_proteina = {}
            
            for header, seq_dna in sequencias_dna.items():
                novo_header = header + f" Proteína Frame {frame}"
                
                if frame <= 3:
                    seq_rna = dna_para_rna(seq_dna)
                    proteina = traduzir_rna(seq_rna, tabela_genetica, frame - 1)
                else:
                    seq_complementar = complemento_reverso(seq_dna)
                    seq_rna = dna_para_rna(seq_complementar)
                    proteina = traduzir_rna(seq_rna, tabela_genetica, frame - 4)
                
                sequencias_proteina[novo_header] = proteina
            
            caminho_proteina = f"c:\\Users\\gabri\\Documents\\Unifesp\\Atividades\\AlgBio\\Atividade 4\\Frame{frame}.fasta"
            salvar_fasta(sequencias_proteina, caminho_proteina)
            print(f"Arquivo salvo: {caminho_proteina}")
        
        print("\nProcessamento concluído com sucesso!")
        print("Arquivos gerados:")
        print("1. Sakai_RNA.fasta - Sequências de RNA")
        print("2-7. Frame1.fasta a Frame6.fasta - Sequências de proteínas nos 6 frames de leitura")
    
    except Exception as e:
        print(f"Erro durante o processamento: {e}")

if __name__ == "__main__":
    main()
