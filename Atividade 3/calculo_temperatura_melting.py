#!/usr/bin/env python
# Nome: Gabriel Pinheiro de Campos
# RA: 156315
# Turma: NB
# Atividade 3: Cálculo de temperatura de melting e conteúdo GC de sequências de DNA

# Importando as bibliotecas necessárias
from Bio import SeqIO
import matplotlib.pyplot as plt
import csv
import math
import os

def calcular_temperatura_melting(conteudo_gc, comprimento, concentracao_sodio=100):
    """
    Calcula a temperatura de melting de uma sequência de DNA.
    
    Fórmula: Tm = 81.5 + 16.6 * log10([Na+]/1000) + 0.41 * (%GC) - 500/comprimento
    
    Parâmetros:
    - conteudo_gc: proporção de G+C na sequência (0 a 1)
    - comprimento: comprimento da sequência em nucleotídeos
    - concentracao_sodio: concentração de Na+ em mM (padrão: 100 mM)
    
    Retorna:
    - Temperatura de melting em °C
    """
    return 81.5 + 16.6 * math.log10(concentracao_sodio/1000) + 0.41 * (conteudo_gc*100) - 500/comprimento

def main():
    # Definindo os arquivos
    arquivo_entrada = "Ecoli_Sakai_cds_from_genomic.fna"
    arquivo_saida_dados = "Dados das sequencias.csv"
    arquivo_saida_gc = "Conteudo_GC.csv"
    arquivo_saida_temp = "Temperatura_x_GC.csv"
    arquivo_grafico = "GC_x_Temperatura.png"
    
    print(f"Processando arquivo: {arquivo_entrada}")
    
    # Verificar se o arquivo existe
    if not os.path.exists(arquivo_entrada):
        print(f"ERRO: O arquivo {arquivo_entrada} não foi encontrado!")
        return
    
    # Armazenar os dados de todas as sequências
    dados_sequencias = []
    
    # Ler o arquivo FASTA
    for registro in SeqIO.parse(arquivo_entrada, "fasta"):
        # Obter identificador e sequência
        identificador = registro.id
        sequencia = str(registro.seq).upper()
        
        # Contar nucleotídeos
        contador_A = sequencia.count("A")
        contador_T = sequencia.count("T")
        contador_C = sequencia.count("C")
        contador_G = sequencia.count("G")
        total = contador_A + contador_T + contador_C + contador_G
        
        # Calcular conteúdo GC
        conteudo_gc = (contador_G + contador_C) / total if total > 0 else 0
        
        # Calcular temperatura de melting
        temperatura_melting = calcular_temperatura_melting(conteudo_gc, total)
        
        # Armazenar os resultados
        dados_sequencias.append({
            "identificador": identificador,
            "contador_A": contador_A,
            "contador_T": contador_T,
            "contador_C": contador_C,
            "contador_G": contador_G,
            "total": total,
            "conteudo_gc": conteudo_gc,
            "temperatura_melting": temperatura_melting
        })
    
    print(f"Foram processadas {len(dados_sequencias)} sequências.")
    
    # Criar arquivo CSV com os dados das sequências
    with open(arquivo_saida_dados, 'w', newline='') as arquivo:
        escritor = csv.writer(arquivo)
        escritor.writerow(["Sequencia", "A", "T", "C", "G", "Total"])
        for dados in dados_sequencias:
            escritor.writerow([
                dados["identificador"],
                dados["contador_A"], 
                dados["contador_T"], 
                dados["contador_C"], 
                dados["contador_G"], 
                dados["total"]
            ])
    
    print(f"Arquivo criado: {arquivo_saida_dados}")
    
    # Criar arquivo CSV com o conteúdo GC
    with open(arquivo_saida_gc, 'w', newline='') as arquivo:
        escritor = csv.writer(arquivo)
        escritor.writerow(["Sequencia", "Conteudo GC"])
        for dados in dados_sequencias:
            escritor.writerow([
                dados["identificador"],
                f"{dados['conteudo_gc']:.3f}"  # 3 casas decimais
            ])
    
    print(f"Arquivo criado: {arquivo_saida_gc}")
    
    # Criar arquivo CSV com temperatura de melting e conteúdo GC
    with open(arquivo_saida_temp, 'w', newline='') as arquivo:
        escritor = csv.writer(arquivo)
        escritor.writerow(["Temperatura de Melting", "Conteudo GC (%)"])
        for dados in dados_sequencias:
            escritor.writerow([
                f"{dados['temperatura_melting']:.3f}",  # 3 casas decimais
                f"{dados['conteudo_gc']*100:.2f}"       # 2 casas decimais
            ])
    
    print(f"Arquivo criado: {arquivo_saida_temp}")
    
    # Criar gráfico de dispersão
    temperaturas = [dados["temperatura_melting"] for dados in dados_sequencias]
    conteudos_gc = [dados["conteudo_gc"] * 100 for dados in dados_sequencias]
    
    plt.figure(figsize=(10, 8))
    plt.scatter(temperaturas, conteudos_gc, alpha=0.5, s=10)
    plt.xlabel('Temperatura de Melting (°C)')
    plt.ylabel('Conteúdo GC (%)')
    plt.title('Relação entre Conteúdo GC e Temperatura de Melting\n(Gabriel Pinheiro de Campos)')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    # Salvar o gráfico como arquivo de imagem
    plt.savefig(arquivo_grafico, dpi=300, bbox_inches='tight')
    print(f"Gráfico salvo em: {arquivo_grafico}")
    
    # Mostrar o gráfico na tela
    plt.show()
    
if __name__ == "__main__":
    main()
