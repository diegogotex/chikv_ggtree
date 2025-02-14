Filogenia CHIKV (em andamento)
================
Diego G. Teixeira
10/04/2020

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

**OBS:** *é importante que no cabeçalho não contenha **espaço** e/ou
**tabs** e nenhum dos seguintes caracteres: “:”, “,”, “)”, “(”, “;”,
“\]”, “\[”, “’”. o RAxML, um dos programas que vamos utilziar pra
inferir a filogenia, é cheio de frescura com o cabeçalho.*

## Alinhamento

Para alinha meu arquivo fasta, eu utilizei a ferramenta
[MAFFT](https://mafft.cbrc.jp/alignment/software/).

``` bat
mafft --reorder --auto CHIKV_Country_OPAL.fasta > CHIKV_Country_OPAL_MAFFT.fasta
```

uma vez alinhado, abri o alinhamento utilizando o
[Aliview](https://ormbunkar.se/aliview/). **OBS:** *Eu costumo utilizar
o Aliview porque é um software leve e com compatibilidade tanto com
Linux, quanto com o macOS e o Windows.*

Uma vez aberto no aliview, o alinhamento vai ficar dessa forma:

<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/ali1.png" width="90%" style="display: block; margin: auto;" />

Agora eu vou “tratar” esse alinhamento. O primeiro passo é transformar
todas as bases mal sequenciadas (geralmente são traocadasas bases
normais \[A,T,C,G\] por o caractére “n”):

```shell
cat CHIKV_Country_OPAL_MAFFT.fasta | sed -e '^[^>]/s/[^ATCGatcg-]/-/g' > CHIKV_Country_OPAL_MAFFT_noN.fasta
```

em seguida eu vou remover as regiões 5’e 3’UTR, veja na imagem que o
alinhamento em ambas as extremidades é ruim por apresentar uma baixa
cobertura e muitas regiões com supostos indels (inserções ou deleções),
as quais, no fim, não irão adicionar muita coisa na inferência da
filogenia. O passo seguinte é remover as regiões de baixa cobertura
dentro da região de interesse no alinhamento. *Resuzir o alinhamento
dessa forma faz o programa de inferência rodar mais rápido, além de
reduzir a quantidade de sítios sem informação de parcimônia (sítios que
contenham pelo menos dois tipos de nucleotídeos diferentes, e que pelo
menos dois deles ocorram com uma freqência mínima de duas vezes). Um
exemplo do que eu estou falando é mostrado a seguir:*

<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/ali2.jpg" title="A seta vermelha indica um exemplo de sítio que não irá alterar a análise filogenética caso seja removido." alt="A seta vermelha indica um exemplo de sítio que não irá alterar a análise filogenética caso seja removido." width="90%" style="display: block; margin: auto;" />

Por fim, irei remover da análise todas as seqências com cumprimento
menor do que 75% do comprimento total do alinhamento. No caso, o nosso
alinhamento apresenta 11.337 sítios, 75% disso corresponde a 8.427.
desse modo, eu irei remover todas as sequências com comprimento inferior
a 8.427 bases nucleotídicas, onde mais de 25% é composto por gaps. Para
realizar essa etapa da análise, utilizarei uma ferramenta chamada
[CIAlign](https://pypi.org/project/cialign/).

```shell
CIAlign --infile CHIKV_Country_OPAL_MAFFT_noN_TRIM.fasta --outfile_stem CHIKV_Country_OPAL_MAFFT_noN_TRIM --remove_short --remove_min_length 8427 --plot_coverage_input --plot_coverage_output
```

Abrindo o alinhamento final no Aliview, nos temos ele desta forma:

<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/ali3.png" width="90%" style="display: block; margin: auto;" />

No final, o dataset ficará com um total de 826 sequências. Para a estapa
de inferência filogenética, exporte o alinhamento em .fasta, .phy e
.nexus.

## Inferindo a filogenia

Em uma análise inicial de inferência filogenética, nós iremos utilizar o
programa [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/). eu
configurei o programa para rodar uma inferência rápida com bootstrap,
com a opção “-f a”, modelo ge substituição GTR e com árvore final
avaliada sob variação gama “-m GTRCAT”. A opção “-x 123” determina a
seed. a “-p 456” número aleatório da seed para a inferência por
parcimônia. “-T 15” é o número de threads que eu utilizei e “-N 500” o
número de replicatas de bootstrap.

```shell
raxmlHPC-PTHREADS -m GTRCAT -f a -x 123 -p 456 -T 15 -N 500 -s CHIKV_Country_OPAL_MAFFT_noN_TRIM_Coverage_cleaned.phy -n CHIKV_100BS
```

Após terminar de processar todas as replicatas, o RAxML vai retornar
diversos arquivos, dentre esses iremos utilizar o RAxML\_bipartitions,
que contém a árvore final com os valores de bootstrap associados às
dicotomias. *Caso você queira trabalhar a topologia da árvore no R, dá
pra importar tanto o arquivo RAxML\_bipartitions quanto o
RAxML\_bootstrap, que a funções dos pacotes ggtree e app conseguem
realizar a compilação de todas as árvores em uma só.* Iremos abrir o
arquivo RAxML\_bibartitions com o programa
[FigTree](http://tree.bio.ed.ac.uk/software/figtree/). Como os genótipos
de CHIKV já foram bem descritos antes, eu irei colorir os ramos de
acordo com os genótipos ECSA, Asian e West African. Dentro do genótipo
ECSA eu vou destacar as amostras isoladas no Brasil e as amostras IOL
(um “subgenótipo” dentro do ECSA).

<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/arvre1.png" width="90%" style="display: block; margin: auto;" />

Na próxima iamgem, eu adicionei o nome das sequências, mas omiti os
valores de suporte de bootstrap nas dicotomias, para evitar uma poluição
visual.

<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/arvre2.png" width="90%" style="display: block; margin: auto;" />

## Filogenia - ggtree

O pacote ggtree apresenta uma das documentações mais mal escritas dentre os pacotes do R que eu tive/quis que trabalhar, mas apesar desse defento, o pacote é muito útil e eficiente no que se propões em gerar e editar topologias de árvores filogenéticas. 

Uma vez que você resolver utilizar, é importante que você saiba identificar o número de cada dicotomia na topologia da árvore, já que parte do processo de anotação da árvore se dá a partir dessa informação. Então, pensando em realizar essa tarefa nós vamos executar os seguintes comandos:

```R
#puxando logo todas as bibliotecas que eu vou usar
library(ggtree)
library(treeio)
library(tidytree)
library(ape)
library(tidyr)
library(ggplot2)
library(phangorn) 
library(ggstance)


#carregando a árvore
arvre <- read.newick("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/RAXML_out/20_03/all/RAxML_bipartitions.CHIKV_100BS", 
                     node.label = "support" #classificando o bootstrap como label
                     )

# uma vez carregada a árvore, a melhor forma de lidar com o dado é transformando-o em objeto "phylo" do pacote ape
# o pacote tydytree tem a função as_tibble para converter o objeto phylo em data frame.
#


#convertendo o objeto treedata em dataframe
arvre_tbl <- as_tibble(arvre)
arvre_tbl[is.na(arvre_tbl$label),]$label <- paste("n",c(1:nrow(arvre_tbl[is.na(arvre_tbl$label),])), sep = "")

#separando uma tabela para o support
arvre_tbl_sup <- as.data.frame(arvre_tbl[,c(4,5)])
colnames(arvre_tbl_sup) <- c("name","support")

#o dataframe permite uma manipulação melhor dos dados. 
#com isso, vamos adicionar a informação a respeito do codom OPAL
head <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/head_OPAL.txt", header = T, stringsAsFactors = F)


#ligando alguma infomação evolutiva à arvore
#a função midpoit seta o midpoint root da árvore
arvre_teste <- midpoint(as.phylo(arvre_tbl))


#vendo o número dos nós
ggtree(arvre_teste) %<+% arvre_tbl + 
  geom_label(aes(x=branch, label=node), size=2.5)

``` 

o output fica assim:
<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/tree_nodes.png" width="100%" style="display: block; margin: auto;" />

Agora que eu sei que eu quero os nós 841 (genótipo ECSA); 1223 (Asian); 828 (African); 904 (IOL); 859 (Brasil), eu vou adicionar essas informações no objeto da árvore e em seguida eu vou adicionar pontos coloridos nas folhas com diferentes codons na posição de transição entre a nsp3 e nsp4.

Vou primeiro plotar a árvore com o codom CGA que é o mais frente, além do códom de parada.
```R

#highlight nos nós para definir os genótipos
g <- ggtree(arvre_teste) + 
  geom_hilight(node = 1223, fill = "tomato3", alpha = .5)+ #Asian genotype
  geom_hilight(node = 841, fill = "springgreen4", alpha = .5)+ #ECSA genotype
  geom_hilight(node = 828, fill = "steelblue", alpha = .5)+ #African
  geom_balance(node = 904, fill = "springgreen3", color = NA)+ #IOL
  geom_balance(node = 859, fill = "yellowgreen", color = NA)+ #Brasil
  geom_treescale(x=0.0001,y=250,offset = 10)

#associando o stop codon opal
#somente o codon mais frequente, o CGA
g <- g %<+% head + 
  geom_tippoint(aes(color=Opal), size=1)+
  scale_color_manual(values = alpha(c("red","black","blue","white", "yellow2","violetred4"),
                                    c(0,1,0,0,0,0)))
g

```
<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/tree_CGA.png" width="100%" style="display: block; margin: auto;" />


Agora colocando todos os codons diferentes.
```R

g1 <- g %<+% head + 
  geom_tippoint(aes(color=Opal), size=1)+
  scale_color_manual(values = alpha(c("red","black","blue","white", "yellow2","violetred4"),
                                    c(1,1,1,0,1,1)))

g1

```
<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/tree_annot.png" width="100%" style="display: block; margin: auto;" />


## Admixture

A realização do teste de admixture é um pouco laborosa, quando o dado de partida é o alinhamento múltiplo de sequências, mas é factível. Para essa análise, o ponto de partida é o alinhamento multiplo de sequências já tratado, após o passo do CIAlign.

O alinhamento vai servir como entrada no software [snp-sites](https://github.com/sanger-pathogens/snp-sites) o qual irá produzor um arquivo VCF.

```shell 
#transformar MSA em vcf
#snp-sites -v -o CHIKV_all_COUNTRY_MAFFT_TRIM_removeN_noREF.vcf CHIKV_all_COUNTRY_MAFFT_TRIM_removeN_noREF.fasta

```

em seguida, precisaremos formatar o arquivo .vcf no formato do [plink](https://www.cog-genomics.org/plink/2.0/), o qual nos vai fornocer o arquivo .ped

```shell
#VCF em .ped
#~/Programs/plink_mac_20200219/plink --vcf CHIKV_all_COUNTRY_MAFFT_TRIM_removeN_noREF.vcf --double-id  --recode12  --out CHIKV
```

uma vez que eu consegui chegar ao arquivo .ped, eu estou apto a executar o programa [ADMIXTURE] (http://dalexander.github.io/admixture/download.html). O ADMIXTURE é um programa que estima a ancestralidade de indivíduos não relacionados, como é descrito no manuscrito do aritgo do software [Alexander2009]. 

Como sugerido no manual do programa, uma etapa importante para a análise é identificar a quantidade de populações ancestreais, postulada pela letra K, realizando um processo chamado de cross-validação. Eu irei testar a validação para um número de 1 até 10 populações diferentes. 

```shell
#testar o melhor K no admixture
for K in 1 2 3 4 5 6 7 8 9 10; do ~/Programs/admixture_macosx-1.3.0/admixture --cv CHIKV.ped $K | tee log${K}.out; done

#pegando os valores de CV
grep -h CV log*.out
```
```bash
CV error (K=1): 0.34922
CV error (K=2): 0.14131
CV error (K=3): 0.09570
CV error (K=4): 0.07196
CV error (K=5): 0.06459
CV error (K=6): 0.06036
CV error (K=7): 0.05869
CV error (K=8): 0.04946
CV error (K=9): 0.04873
CV error (K=10): 0.04782
```

após a execução pro programa, eu vou importar os resultados no R e gerar um gráfico para o CV. **OBS**: esse passo não é necessário, uma vez que você já tem os valores e consegue decidir o valor de K somente por eles, mas ainda assim eu gosto de realizá-lo.

```R

CV_error <- as.data.frame(matrix(nrow = 10, ncol = 2))
colnames(CV_error) <- c("K","CV")
CV_error$K <- c(1:10)
CV_error$CV <- c(0.34922, 0.14131, 0.09570, 0.07196, 0.06459, 0.06036, 0.05869, 0.04946, 0.04873, 0.04782)

plot(CV_error$K, 
     CV_error$CV, 
     type = "b", #both, linha e ponto
     xlab = "K",
     ylab = "Cross-Validation error",
     xaxt = "n")
axis(1, at = c(1:10))

```
<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/admixture_cv.png" width="80%" style="display: block; margin: auto;" />


Agora que eu já conheço o CV para cada um das "quantidades de ancestrais" eu posso executar o admixture para cada um dos __K__ que eu achar representativo para a minha amostra. 


```R
#puxando o arquivo
admixK3 <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/K3/CHIKV.3.Q")
#puxando o nome das amostras
admixK3$sample <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/CHIKV.head.txt", header = F)$V1
#dando nome às colunas relacionado ao genótipo 
colnames(admixK3) <- c("Asian","ECSA","West African","id")
#reordenando as colunas
admixK3<- admixK3[,c(4,1:3)]

#plotando
barplot(t(as.matrix(admixK3[,2:4])),
        col= rainbow(3),
        xlab="Virus", 
        ylab = "Ancestry",
        space = 0,
        border = NA
)
```
K3</br>
<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/admix_K3.png" width="80%" style="display: block; margin: auto;" />


```R
#Já expliquei o código acima, aqui vou só jogar o código dos outros plots mesmo.

admixK4 <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/K4/CHIKV.4.Q")
admixK4$sample <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/CHIKV.head.txt", header = F)$V1
colnames(admixK4) <- c("Asian","ECSA","West African","Mix","id")
admixK4<- admixK4[,c(5,1:4)]


barplot(t(as.matrix(admixK4[,2:5])),
        col= rainbow(4),
        xlab="Virus", 
        ylab = "Ancestry",
        space = 0,
        border = NA
        )



admixK5 <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/K5/CHIKV.5.Q")
admixK5$sample <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/CHIKV.head.txt", header = F)$V1
colnames(admixK5) <- c("Asian","ECSA","West African","Mix","Mix2","id")
admixK5<- admixK5[,c(6,1:5)]


barplot(t(as.matrix(admixK5[,2:6])),
        col= rainbow(4),
        xlab="Virus", 
        ylab = "Ancestry",
        space = 0,
        border = NA
)


admixK6 <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/K6/CHIKV.6.Q")
admixK6$sample <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/CHIKV.head.txt", header = F)$V1
colnames(admixK6) <- c("Asian","ECSA","West African","Mix","Mix2","Mix3","id")
admixK6<- admixK6[,c(7,1:6)]


barplot(t(as.matrix(admixK6[,2:7])),
        col= rainbow(6),
        xlab="Virus", 
        ylab = "Ancestry",
        space = 0,
        border = NA
)



admixK7 <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/K7/CHIKV.7.Q")
admixK7$sample <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/CHIKV.head.txt", header = F)$V1
colnames(admixK7) <- c("Asian","ECSA","West African","Mix","Mix2","Mix3","Mix4","id")
admixK7<- admixK7[,c(8,1:7)]


barplot(t(as.matrix(admixK7[,2:8])),
        col= rainbow(7),
        xlab="Virus", 
        ylab = "Ancestry",
        space = 0,
        border = NA
)



admixK8 <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/K8/CHIKV.8.Q")
admixK8$sample <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/CHIKV.head.txt", header = F)$V1
colnames(admixK8) <- c("Asian","ECSA","West African","Mix","Mix2","Mix3","Mix4","Mix5","id")
admixK8<- admixK8[,c(9,1:8)]


barplot(t(as.matrix(admixK8[,2:9])),
        col= rainbow(8),
        xlab="Virus", 
        ylab = "Ancestry",
        space = 0,
        border = NA
)

admixK9 <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/K9/CHIKV.9.Q")
admixK9$sample <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/CHIKV.head.txt", header = F)$V1
colnames(admixK9) <- c("Asian","ECSA","West African","Mix","Mix2","Mix3","Mix4","Mix5","Mix6","id")
admixK9<- admixK9[,c(10,1:9)]


barplot(t(as.matrix(admixK9[,2:10])),
        col= rainbow(9),
        xlab="Virus", 
        ylab = "Ancestry",
        space = 0,
        border = NA
)


admixK10 <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/K10/CHIKV.10.Q")
admixK10$sample <- read.table("~/Dropbox/Doutorado/CHIKV_PHYLO/Country/CHIKV_Country_OPAL_MAFFT_noN_TRIM/vcf/admixture/CHIKV.head.txt", header = F)$V1
colnames(admixK10) <- c("Asian","ECSA","West African","Mix","Mix2","Mix3","Mix4","Mix5","Mix6","Mix7","id")
admixK10<- admixK10[,c(11,1:10)]


barplot(t(as.matrix(admixK10[,2:11])),
        col= rainbow(10),
        xlab="Virus", 
        ylab = "Ancestry",
        space = 0,
        border = NA
)

```
K4</br>
<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/admix_K4.png" width="80%" style="display: block; margin: auto;" />

K5</br>
<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/admix_K5.png" width="80%" style="display: block; margin: auto;" />

K6</br>
<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/admix_K6.png" width="80%" style="display: block; margin: auto;" />

K7</br>
<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/admix_K7.png" width="80%" style="display: block; margin: auto;" />

K8</br>
<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/admix_K8.png" width="80%" style="display: block; margin: auto;" />

K9</br>
<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/admix_K9.png" width="80%" style="display: block; margin: auto;" />

K10</br>
<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/admix_K10.png" width="80%" style="display: block; margin: auto;" />


Essa visualização padrão tanto do Admixture quanto do Structure são legais, mas infelizmente a gente não consegue distinguir direito a qual sequência perntece a barra, e adicionar o nome à cada uma deixa o gráfico ainda mais poluído. Como você pode observar, a medida que vamos aumentando o número de K, também vai aumentando a complexidade de que interpretar o gráfico. Por causa disso, eu vou elencar somente 3 gráficos para continuar a análise:

**1. K3 (quantidade de gentótipos diferentes identificados para CHIKV);**</br>
**2. K4 (Valor cd CV baixo, próximo ao de K5, K6 e K7);**</br>
**3. K8 (Valor de CV baixo, similar ao K9 e K10);**

## Filogenia + Admixture

Para ajudar na visualização do Admixture, eu resolvi plotar os gráficos juntos a filogenia. A melhor forma de fazer isso é por meio do ggtree, já que o Figtree não permite esse tipo de visualização e eu não consigo imaginar o tamanho do trabalho que vai dar pra alinhar as colunas às folhas da topologia de forma manual. 

a primeira coisa que eu devo fazer é transformar a tabela do admixture de modo que seja lida como um objeto do ggplot (pacote gráfico associado ao ggtree), para que eu associe o valores às folhas da ávore por meio no nome da sequência.

```R
#convertendo as informacões do dataframe de um modo que o ggplot2 compreenda
admix3_re <- reshape2::melt(admixK3)
admix4_re <- reshape2::melt(admixK4)
admix8_re <- reshape2::melt(admixK8)
```

Uma vez que eu transformei o dataframe, eu posso associar ao objeto da árvore no ggtree:

```R
#plotando com o admixture K3
g <- facet_plot(g, panel = "Admixture K3", data = admix3_re,
                geom = geom_barh,
                mapping = aes(x = value, fill=as.factor(variable)),
                stat = 'identity',
                width=1)

#plotando com o admixture K4
g <- facet_plot(g, panel = "Admixture K4", data = admix4_re,
                geom = geom_barh,
                mapping = aes(x = value, fill=as.factor(variable)),
                stat = 'identity',
                width=1)

g <- facet_plot(g, panel = "Admixture K8", data = admixK8_re,
                geom = geom_barh,
                mapping = aes(x = value, fill=as.factor(variable)),
                stat = 'identity',
                width=1)

g
```

<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/arvre_admix.png" width="80%" style="display: block; margin: auto;" />



## teste de Recombinação
