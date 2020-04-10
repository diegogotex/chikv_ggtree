Filogenia CHIKV (em andamento)
================
Diego G. Teixeira
09/04/2020

## Introdução

Nesse estudo eu estou investigando a Filogenia do Vírus Chikungunya
(CHIKV). O CHIKV é um Alphavirus, da família Togaviridae, normalmente
definido por possuir um genoma de média de 11.5 kb, composto por duas
janelas de leitura (Open Read Frame - ORF). A primeira ORF comporta 4
proteínas não estruturais nsP1, nsP2, nsP3 e nsP4; A segunda ORF
comporta os genes responsáveis pela estrutura do virion, os genes C, E1,
E2, 6K/TF e E3.

<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/Fig_genoma.png" title="Estrutura genômica do CHIKV" alt="Estrutura genômica do CHIKV" width="90%" style="display: block; margin: auto;" />

Apesar de ser caracterizado como um vírus que possui 2 ORFs e o genoma
de referência, disponível no NCBI
([NC\_004162](https://www.ncbi.nlm.nih.gov/nuccore/NC_004162.2)),
corroborar com a afirmação, grande parte dos genomas completos de CHIKV
disponíveis no NCBI apresentam um códon de parada entre os genes
codificadores das proteínas nsP3 e nsP4. Inclusive, o genoma do vírus
isolados de pessoas diferentes em Natal, também mostra que há o códon de
parada nos isolados daqui.

<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/opal_genome.png" title="Códon de parada (destacado pela linha horizontal vermelha) em 3 isolados do Brasil" alt="Códon de parada (destacado pela linha horizontal vermelha) em 3 isolados do Brasil" width="90%" style="display: block; margin: auto;" />

A partir desse alinhamento, tive a ideia de observar a frequência do
códon de parada em todos os genomas completos de CHIKV já publicados.

## Obtenção das sequências completas de CHIKV

### NCBI

1.  O primeiro passo para observar a frequência do códom de parada nos
    genomas de CHIVK, é baixar os genomas no
    [NCBI](https://www.ncbi.nlm.nih.gov/);

2.  Em seguida, nós pesquisamos por “Chikungunya” no banco de dados de
    “nucleotide”;

3.  selecionamos o “Chikungunya virus” no campo “Results by taxon”;
    
    3.1. O campo de busca vai ficar dessa forma: “(chikungunya)
    AND”Chikungunya virus“\[porgn:\_\_txid37124\]”;

4.  O passo seguinte é selecionar apenas as sequências que tenham pelo
    menos 10.5kb de comprimento. Eu selecionei 13kb como comprimento
    máximo, mas isso foi apenas por causa do NCBI exigir que esse campo
    seja preenchido.
    
    <img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/ncbi_length.png" width="90%" style="display: block; margin: auto;" />

5.  por fim, vamos baixar o arquivo GenBank selecionando a opção Send to
    \> Format \> GenBank (full) e clicando em “Create File”.

Ao final dessa etapa, nós teremos um arquivo .gb com 958 entradas. O
arquivo desse trabalho foi baixado no dia 1/03/2020, é provável que esse
número aumente no decorrer do tempo.

### Formatando o arguivo

depois baixar o arguivo .gb, eu realizei uma formatação do dados,
transformando o arguivo .gb em um arquivo tabular com algumas
informações que eu achei necessárias para o estudo. Para formatar essa
tabela, eu utilizei um código em Python (v.3.7):

``` python
# -*- coding: utf-8 -*-

#carregando as bibliotecas
from Bio import SeqIO

#carregando o arquivo do genbank
#o ACC MK518395 tem um erro e a informação Assembly-Data tem que ser removida 
gbFile = "/Users/diegogotex/Dropbox/Doutorado/CHIKV_PHYLO/CHIKV_10_5kb_to_13kb_01_mar_20.gb"
#02/03/2020
#loop de teste para saber se deu certo carregar o arquivo
count = 0
for seq_record in SeqIO.parse(gbFile, "genbank"):
    #print(seq_record.id)
    count +=1
    
print(count)

#criando um arquivo de escrita
file = open("/Users/diegogotex/Dropbox/Doutorado/CHIKV_PHYLO/CHIKV_>10.5kb.tab", "w+")

#loop2 
for seq_record in SeqIO.parse(gbFile, "genbank"): 
    source = seq_record.features[0] #puxando as features de cada ACCESSION
    country = str(source.qualifiers.get('country')) #
    country = country.replace(" ", "_") #removendo os espaços
    country = country.strip("[']") #removendo a aspa simples
    date = str(source.qualifiers.get('collection_date')) #puxando a data de coleta
    date = date.strip("[']") #removendo a aspa simples
    #printando as variáveis separadas por tabulação
    #print("%s\t%s\t%s\t%s" % (seq_record.id, country, date, seq_record.seq))
    file.write("%s\t%s\t%s\t%s" % (seq_record.id, country, date, seq_record.seq))
    file.write("\n")
    
file.close
```

**OBS:** *inicialmente eu tive alguns erros durante o processo de
formatação, no passo do loop2, por causa de um erro nas informações da
sequência MK518395. por algum motivo, quando o loop chegava nela dava
erro. inicialmente eu removi entrada e rodei o script, aí deu certo. mas
em seguida, eu readicionei a entrada e saí removendo cada uma das
features até identificar que o que estava causando a interrupção na
execução do scritp estava na informação Assembly-Data.*

Após o processo de transformação do arquivo .gb (se quiser baixar o
mesmo arquivo que eu trabalhei é esse aqui
[CHIKV\_10\_5kb\_to\_13kb\_01\_mar\_20.gb](https://github.com/diegogotex/chikv_ggtree/blob/master/files/CHIKV_10_5kb_to_13kb_01_mar_20.gb))
em um arquivo tabular, eu realizei a etapa restante de formatação das
informações de forma manual. Inseri a tabela no LibreOffice Calc (TALEZ
o excell sirva…) e removi as sequências que não continham informação de
País na qual a amopstra foi coletada, na entrada do genbank, passando de
958 entradas para 880. Eu filtrei essas sequências, pois a minha
intenção nessa análise é identificar a presença do códon OPAL nos
genomas e procurar alguma correlação com o local de infecção.

<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/tab_file.png" width="90%" style="display: block; margin: auto;" />

Eu gosto de montar um arquivo com múltiplas planilhas, assim eu mantenho
todo o fluxo de trabalho e sei quais sequências saíram e o motivo delas
terem saído. Claro, você pode fazer o mesmo por meio de algum script,
mas eu prefiro trabalhar dessa forma porquê eu tenho um contato maior
com as sequências.

O passo segunte é transformar esta tabela em um arquivo fasta. No meu
arquivo, eu mantive apenas o ACC, País e o códom na posição OPAL dos
genomas (ACC|País|Códon), como cabeçalho. Ex.: \>LC259094.1|Angola|TGA.

## Alinhamento

Para alinha meu arquivo fasta, eu utilizei a ferramenta
[MAFFT](https://mafft.cbrc.jp/alignment/software/).

``` bat
mafft --reorder --auto CHIKV_Country_OPAL.fasta > CHIKV_Country_OPAL_MAFFT.fasta
```

uma vez alinhado, abri o alinhamento utilizando o
[Aliview](https://ormbunkar.se/aliview/). **OBS:**\_Eu costumo utilizar
o Aliview porque é um software leve e com compatibilidade tanto com
Linux, quanto com o macOS e o Windwos.\_

Uma vez aberto no aliview, o alinhamento vai ficar dessa forma:

Ali1

Agora eu vou “tratar” esse alinhamento. O primeiro passo é transformar
todas as bases mal sequenciadas (geralmente são traocadasas bases
normais \[A,T,C,G\] por o caractére “n”):

``` bat
cat CHIKV_Country_OPAL_MAFFT.fasta | sed -e '^[^>]/s/[^ATCGatcg-]/-/g' > CHIKV_Country_OPAL_MAFFT_noN.fasta
```

em seguida eu vou remover as regiões 5’e 3’UTR, veja na imagem que o
alinhamento em ambas as extremidades é ruim por apresentar uma baixa
cobertura e muitas regiões com supostos indels (inserções ou deleções),
as quais, no fim, não irão adicionar muita coisa na inferência da
filogenia. O passo seguinte é remover as regiões de baixa cobertura
dentro da região de interesse no alinhamento. *Resuzir o alinhamento
dessa forma faz o programa de inferência rodar mais rápido,além de
reduzir a quantidade de sítios sem informação de parcimônia (sítios que
contenham pelo menos dois tipos de nucleotídeos diferentes, e que pelo
menos dois deles ocorram com uma freqência mínima de duas vezes).*

## Inferindo a filogenia

## Admixture

## Filogenia + Admixture

## teste de Recombinação
