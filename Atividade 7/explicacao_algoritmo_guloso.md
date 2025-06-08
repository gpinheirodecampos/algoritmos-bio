# Algoritmo Guloso para Montagem de Contig

## O que é um Algoritmo Guloso?

Um algoritmo guloso (greedy algorithm) é um tipo de algoritmo que resolve problemas fazendo sempre a escolha que parece melhor no momento atual, sem se preocupar com as consequências futuras. É como se fosse uma pessoa "gananciosa" que sempre pega a melhor opção disponível na hora, mesmo que isso não seja o melhor a longo prazo.

## Aplicação na Montagem de Genoma

Na montagem de genoma, o problema é juntar várias sequências pequenas de DNA (as reads) para formar uma sequência maior e contínua (o contig). É como montar um quebra-cabeça onde as peças são trechos de DNA que se sobrepõem.

O algoritmo guloso funciona assim:
1. **Olha todas as reads** que ainda não foram juntadas
2. **Procura quais reads têm a maior sobreposição** entre si
3. **Escolhe o par com maior sobreposição** (essa é a decisão "gulosa")
4. **Junta essas duas reads** em uma sequência maior
5. **Repete o processo** até sobrar apenas uma sequência final

## Como Funciona na Prática

Com as reads do nosso exercício:
- read_1: ATTAGACCTG
- read_2: CCTGCCGGAA  
- read_3: AGACCTGCCG
- read_4: GCCGGAATAC

**Passo 1:** O algoritmo compara todas as combinações e encontra que ATTAGACCTG e AGACCTGCCG têm 7 nucleotídeos iguais no final/início (AGACCTG). Como essa é a maior sobreposição, ele junta essas duas: ATTAGACCTGCCG

**Passo 2:** Agora compara as reads restantes e vê que CCTGCCGGAA e GCCGGAATAC também têm 7 nucleotídeos de sobreposição (GCCGGAA). Junta elas: CCTGCCGGAATAC

**Passo 3:** Finalmente, junta os dois pedaços maiores porque eles também se sobrepõem em 7 nucleotídeos (CCTGCCG), resultando no contig final: ATTAGACCTGCCGGAATAC

## Por que é "Guloso"?

O algoritmo é chamado de "guloso" porque a cada passo ele sempre escolhe a maior sobreposição disponível naquele momento, sem pensar se essa escolha vai atrapalhar depois. É uma estratégia simples e rápida, mas nem sempre garante o melhor resultado possível.

## Vantagens e Desvantagens

**Vantagens:**
- É fácil de entender e implementar
- Funciona rápido
- Dá bons resultados na maioria dos casos

**Desvantagens:**
- Pode não encontrar a melhor solução possível
- Se houver erros nas reads, pode montar errado
- Não consegue "voltar atrás" se fizer uma escolha ruim

## Conclusão

O algoritmo guloso é uma ferramenta útil para montagem de genoma porque é simples e eficiente. Mesmo não sendo perfeito, ele consegue resolver bem o problema de juntar reads para formar contigs, que é essencial no sequenciamento de DNA.
