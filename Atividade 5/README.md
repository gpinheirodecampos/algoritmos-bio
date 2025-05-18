# Atividade 5 - Anotação Genômica

## Resposta Dissertativa

Esta atividade é parte de uma operação chamada de "Anotação" no sequenciamento de um genoma. A anotação genômica é o processo de identificar e rotular todos os elementos funcionais dentro de um genoma sequenciado. Esta atividade específica se concentra na identificação de regiões codificadoras de proteínas (ORFs - Open Reading Frames) nos seis frames possíveis de leitura.

Após esta operação de identificação de potenciais ORFs, a próxima operação a ser feita para finalizar a anotação seria a **caracterização funcional** das proteínas identificadas. Isso envolve:

1. **Análise de homologia**: Comparar as sequências de proteínas encontradas com bancos de dados de proteínas conhecidas (como o UniProt, NCBI nr, etc.) para identificar similaridades com proteínas já caracterizadas.

2. **Predição de domínios funcionais**: Utilizar ferramentas como Pfam, PROSITE, ou InterPro para identificar domínios conservados e motivos funcionais nas sequências.

3. **Análise de estrutura secundária e terciária**: Prever a estrutura da proteína para obter insights sobre sua função.

4. **Análise de contexto genômico**: Examinar a localização e o contexto das ORFs no genoma para identificar operons, ilhas genômicas, ou outras características funcionalmente relevantes.

5. **Validação experimental**: Em alguns casos, realizar experimentos para confirmar a funcionalidade das proteínas preditas.

A atividade que realizamos aqui (identificação de ORFs) se relaciona diretamente com esta próxima etapa porque fornece os "candidatos" - as sequências de proteínas potenciais que serão submetidas à caracterização funcional. Sem a identificação correta das ORFs, a anotação funcional não teria sequências para analisar.

Além disso, os dados que obtivemos (localização das ORFs no genoma, tamanho das proteínas, distribuição nos diferentes frames) já dão pistas importantes sobre características gerais do genoma, como densidade gênica, uso preferencial de códons, e organização genômica, que são informações valiosas para o processo completo de anotação.

## Código Python

O script `protein_finder.py` implementa a identificação de potenciais regiões codificadoras de proteínas nos 6 frames de leitura do genoma de E. coli Sakai, incluindo o genoma principal e os dois plasmídeos.

## Arquivos Gerados

Após a execução do script, serão gerados 18 arquivos FASTA:

- Genoma principal (BA000007.3):
  - gc_frame1_ativ5_XXXXX.fasta a gc_frame6_ativ5_XXXXX.fasta

- Plasmídeo 1 (AB011549.2):
  - p1_frame1_ativ5_XXXXX.fasta a p1_frame6_ativ5_XXXXX.fasta

- Plasmídeo 2 (AB011548.2):
  - p2_frame1_ativ5_XXXXX.fasta a p2_frame6_ativ5_XXXXX.fasta

Onde XXXXX deve ser substituído pelo número de RA do estudante.

## Como Executar

```bash
python protein_finder.py
```

Certifique-se de que as bibliotecas necessárias estão instaladas:

```bash
pip install biopython
```

## Observações

- O código considera apenas ORFs que começam com o códon ATG (Metionina) e terminam com um dos códons de parada (TAA, TAG, TGA).
- Apenas ORFs com pelo menos 50 aminoácidos são relatadas, para evitar muitos falsos positivos.
- A localização das ORFs é indicada no cabeçalho de cada sequência, conforme a convenção solicitada.
