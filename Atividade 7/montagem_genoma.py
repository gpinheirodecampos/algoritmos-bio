"""
Algoritmo Guloso para Montagem de Genoma
=========================================

Este programa implementa um algoritmo guloso para montar um contig a partir de reads de DNA.
O algoritmo funciona encontrando as maiores sobreposições entre sequências e as fundindo
iterativamente até formar uma única sequência contígua.

Algoritmo Guloso Aplicado à Montagem de Contig:
-----------------------------------------------
1. LEITURA: Carrega todas as reads do arquivo FASTA usando BioPython
2. BUSCA DE SOBREPOSIÇÕES: Para cada par de reads, calcula a maior sobreposição possível
3. SELEÇÃO GULOSA: Escolhe o par com maior sobreposição (estratégia gulosa)
4. FUSÃO: Une as duas reads removendo a região sobreposta
5. ITERAÇÃO: Repete o processo até não haver mais sobreposições significativas

O algoritmo é chamado "guloso" porque a cada passo toma a decisão localmente ótima
(maior sobreposição disponível) sem considerar se essa escolha levará ao resultado
globalmente ótimo. Esta abordagem é eficiente e produz bons resultados na prática
para montagem de genomas.

Autor: Gabriel Pinheiro de Campos
RA: 156315
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def ler_fasta(arquivo):
    sequencias = {}
    
    # Usa BioPython para ler o arquivo FASTA
    for record in SeqIO.parse(arquivo, "fasta"):
        sequencias[record.id] = str(record.seq)
    
    return sequencias

def encontrar_sobreposicao(seq1, seq2, min_overlap=3):
    max_overlap = 0
    best_merge = ""
    best_type = ""
    
    # Testa seq1 seguida de seq2 (sufixo de seq1 sobrepõe com prefixo de seq2)
    for i in range(min_overlap, min(len(seq1), len(seq2)) + 1):
        if seq1[-i:] == seq2[:i]:  # Sufixo de seq1 == Prefixo de seq2
            if i > max_overlap:
                max_overlap = i
                best_merge = seq1 + seq2[i:]  # Remove a parte sobreposta de seq2
                best_type = "seq1_seq2"
    
    # Testa seq2 seguida de seq1 (sufixo de seq2 sobrepõe com prefixo de seq1)
    for i in range(min_overlap, min(len(seq1), len(seq2)) + 1):
        if seq2[-i:] == seq1[:i]:  # Sufixo de seq2 == Prefixo de seq1
            if i > max_overlap:
                max_overlap = i
                best_merge = seq2 + seq1[i:]  # Remove a parte sobreposta de seq1
                best_type = "seq2_seq1"
    
    return max_overlap, best_type, best_merge

def algoritmo_guloso_montagem(sequencias):
    # Converte para lista de sequências para facilitar manipulação
    reads = list(sequencias.values())
    identificadores = list(sequencias.keys())
    
    print(f"Iniciando montagem com {len(reads)} reads:")
    for i, (id_read, seq) in enumerate(zip(identificadores, reads)):
        print(f"  {id_read}: {seq}")
    print()
    
    passo = 1
    
    # Enquanto temos mais de uma sequência, procura a melhor fusão
    while len(reads) > 1:
        melhor_sobreposicao = 0
        melhor_fusao = ""
        indices_fusao = (-1, -1)
        
        print(f"Passo {passo}: Buscando melhor sobreposição...")
        
        # Compara todas as combinações de pares
        for i in range(len(reads)):
            for j in range(i + 1, len(reads)):
                overlap, tipo, fusao = encontrar_sobreposicao(reads[i], reads[j])
                
                print(f"  {reads[i]} + {reads[j]} -> sobreposição: {overlap}")
                
                if overlap > melhor_sobreposicao:
                    melhor_sobreposicao = overlap
                    melhor_fusao = fusao
                    indices_fusao = (i, j)
        
        # Se encontrou sobreposição, faz a fusão
        if melhor_sobreposicao > 0:
            i, j = indices_fusao
            print(f"  MELHOR: {reads[i]} + {reads[j]} (sobreposição: {melhor_sobreposicao})")
            print(f"  RESULTADO: {melhor_fusao}")
            
            # Remove as sequências originais (remove o índice maior primeiro)
            reads.pop(max(i, j))
            reads.pop(min(i, j))
            
            # Adiciona a nova sequência fundida
            reads.append(melhor_fusao)
            
            print(f"  Reads restantes: {len(reads)}")
            print()
        else:
            # Não há mais sobreposições, concatena as sequências restantes
            print("  Nenhuma sobreposição encontrada. Concatenando sequências restantes...")
            resultado = "".join(reads)
            break
        
        passo += 1
    
    # Se só resta uma sequência, ela é o contig final
    if len(reads) == 1:
        resultado = reads[0]
    
    return resultado

def salvar_contig_fasta(contig, arquivo_saida, ra_estudante="156315"):
    # Cria um objeto SeqRecord com a sequência do contig
    seq_record = SeqRecord(
        Seq(contig),
        id="contig_guloso",
        description=f"Contig montado utilizando Algoritmo Guloso {ra_estudante} - Sequencia gerada a partir do algoritmo guloso"
    )
    
    # Salva usando BioPython
    with open(arquivo_saida, 'w') as f:
        SeqIO.write(seq_record, f, "fasta")

def main():
    print("=" * 60)
    print("ALGORITMO GULOSO PARA MONTAGEM DE GENOMA")
    print("=" * 60)
    print()

    arquivo_entrada = "reads4.fasta"
    arquivo_saida = "contig.fasta"
    
    try:
        print("1. Lendo sequências do arquivo FASTA...")
        sequencias = ler_fasta(arquivo_entrada)
        print(f"   {len(sequencias)} sequências carregadas com sucesso!")
        print()
        
        print("2. Aplicando algoritmo guloso para montagem...")
        contig_final = algoritmo_guloso_montagem(sequencias)
        print()
        
        print("3. RESULTADO DA MONTAGEM:")
        print(f"   Contig final: {contig_final}")
        print(f"   Tamanho: {len(contig_final)} nucleotídeos")
        print()
        
        print("4. Salvando contig em arquivo FASTA...")
        salvar_contig_fasta(contig_final, arquivo_saida)
        print(f"   Contig salvo em: {arquivo_saida}")
        
        print()
        print("=" * 60)
        print("MONTAGEM CONCLUÍDA COM SUCESSO")
        print("=" * 60)
        
    except FileNotFoundError:
        print(f"ERRO: Arquivo '{arquivo_entrada}' não encontrado")
    except Exception as e:
        print(f"ERRO: {e}")

if __name__ == "__main__":
    main()
