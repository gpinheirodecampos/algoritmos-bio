#!/usr/bin/env python
# Finder de proteínas - Atividade 5
# Identifica regiões do genoma da E.coli Sakai que potencialmente codificam proteínas

from Bio import SeqIO
from Bio.Seq import Seq
import os

NUMERO_RA = "156315"

# Dicionário com os códons de início e parada
GENETIC_CODE = {
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
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

# Códons de início e parada
START_CODON = 'ATG'  # Metionina (M)
STOP_CODONS = ['TAA', 'TAG', 'TGA']  # Codificados como '*'

def read_fasta(file_path):
    """Lê um arquivo FASTA e retorna o registro"""
    records = list(SeqIO.parse(file_path, "fasta"))
    if records:
        return records[0]
    return None

def translate_sequence(seq, frame):
    """
    Traduz a sequência de DNA para proteína no frame especificado
    frame: 1, 2, 3 para fita positiva e 4, 5, 6 para fita negativa
    """
    if frame <= 3:
        # Frames 1, 2, 3 (fita positiva)
        start_pos = frame - 1
        coding_seq = seq[start_pos:]
        protein_seq = coding_seq.translate()
    else:
        # Frames 4, 5, 6 (fita negativa / complementar)
        start_pos = abs(frame) - 4
        coding_seq = seq.reverse_complement()[start_pos:]
        protein_seq = coding_seq.translate()
    
    return protein_seq

def find_proteins(dna_seq, frame, min_length=50):
    """
    Encontra sequências de proteínas que começam com M (ATG) e terminam com um códon de parada.
    frame: 1-6 (1-3 fita +, 4-6 fita -)
    min_length: tamanho mínimo da proteína em aminoácidos
    Retorna lista de tuplas (proteína, posição_inicio, posição_fim)
    """
    proteins = []
    
    # Traduz a sequência no frame especificado
    if frame <= 3:
        # Frames 1, 2, 3 (fita positiva)
        start_pos = frame - 1
        dna_frame = dna_seq[start_pos:]
        protein_seq = str(dna_frame.translate())
        
        # Encontra todas as proteínas que começam com M e terminam com *
        m_positions = [i for i, aa in enumerate(protein_seq) if aa == 'M']
        
        for m_pos in m_positions:
            # Procura o próximo stop codon depois da metionina
            stop_pos = protein_seq.find('*', m_pos)
            if stop_pos > m_pos:
                protein = protein_seq[m_pos:stop_pos]
                
                # Verifica se a proteína tem o tamanho mínimo
                if len(protein) >= min_length:
                    # Calcula as posições no genoma original
                    genome_start = start_pos + (m_pos * 3)
                    genome_end = start_pos + (stop_pos * 3) + 2  # +2 para incluir o códon de parada
                    proteins.append((protein, genome_start, genome_end))
    else:
        # Frames 4, 5, 6 (fita negativa / complementar)
        start_pos = abs(frame) - 4
        dna_frame = dna_seq.reverse_complement()[start_pos:]
        protein_seq = str(dna_frame.translate())
        
        # Encontra todas as proteínas que começam com M e terminam com *
        m_positions = [i for i, aa in enumerate(protein_seq) if aa == 'M']
        
        for m_pos in m_positions:
            # Procura o próximo stop codon depois da metionina
            stop_pos = protein_seq.find('*', m_pos)
            if stop_pos > m_pos:
                protein = protein_seq[m_pos:stop_pos]
                
                # Verifica se a proteína tem o tamanho mínimo
                if len(protein) >= min_length:
                    # Calcula as posições no genoma original (invertido para a fita complementar)
                    genome_start = len(dna_seq) - (start_pos + (stop_pos * 3) + 2)
                    genome_end = len(dna_seq) - (start_pos + (m_pos * 3))
                    proteins.append((protein, genome_start, genome_end))
    
    return proteins

def save_to_fasta(proteins, record_id, frame, output_prefix, ra_number):
    """
    Salva as proteínas encontradas em um arquivo FASTA
    """
    if not proteins:
        print(f"Nenhuma proteína encontrada para {record_id} no frame {frame}")
        # Cria um arquivo vazio
        with open(f"{output_prefix}_frame{frame}_ativ5_{ra_number}.fasta", "w") as f:
            f.write("")
        return
    
    with open(f"{output_prefix}_frame{frame}_ativ5_{ra_number}.fasta", "w") as f:
        for i, (protein, start, end) in enumerate(proteins, 1):
            header = f">{record_id}, Frame {frame}, proteína {i}, [location={start+1}..{end+1}], {ra_number}"
            f.write(f"{header}\n{protein}\n")
    
    print(f"Salvou {len(proteins)} proteínas para {record_id} no frame {frame}")

def process_genome(file_path, prefix):
    """
    Processa um arquivo de genoma para todos os frames e salva os resultados
    """
    record = read_fasta(file_path)
    if not record:
        print(f"Erro ao ler o arquivo {file_path}")
        return
    
    for frame in range(1, 7):
        proteins = find_proteins(record.seq, frame)
        save_to_fasta(proteins, record.id, frame, prefix, NUMERO_RA)

def main():
    # Define os arquivos de entrada e prefixos de saída
    files_and_prefixes = [
        ("Ecoli Sakai sequence.fasta", "gc"),  # Genoma completo
        ("Ecoli Sakai Plasmid 1 sequence.fasta", "p1"),  # Plasmídeo 1
        ("Ecoli Sakai plasmid 2 sequence.fasta", "p2")  # Plasmídeo 2
    ]
    
    # Verifica se os arquivos existem
    dir_path = os.path.dirname(os.path.abspath(__file__))
    print(f"Diretório atual: {dir_path}")
    dir_files = os.listdir(dir_path)
    print(f"Arquivos no diretório: {dir_files}")
    
    # Processa cada arquivo
    for file_name, prefix in files_and_prefixes:
        file_path = os.path.join(dir_path, file_name)
        if not os.path.exists(file_path):
            print(f"ERRO: Arquivo {file_path} não encontrado!")
            continue
            
        print(f"Processando {file_name}...")
        record = read_fasta(file_path)
        if record:
            print(f"Registro lido: {record.id}, tamanho: {len(record.seq)} pb")
            process_genome(file_path, prefix)
        else:
            print(f"ERRO: Não foi possível ler o registro do arquivo {file_path}")
    
    print("\nProcessamento concluído!")
    print("Arquivos .fasta gerados com as potenciais proteínas em cada frame.")

if __name__ == "__main__":
    main()
