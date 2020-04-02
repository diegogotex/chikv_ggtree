Filogenia CHIKV (em andamento)
================
Diego G. Teixeira
02/04/2020

Introdução
----------

Nesse estudo eu estou investigando a Filogenia do Vírus Chikungunya (CHIKV). O CHIKV é um Alphavirus, da família Togaviridae, normalmente definido por possuir um genoma de média de 11.5 kb, composto por duas janelas de leitura (Open Read Frame - ORF). A primeira ORF comporta 4 proteínas não estruturais nsP1, nsP2, nsP3 e nsP4; A segunda ORF comporta os genes responsáveis pela estrutura do virion, os genes C, E1, E2, 6K/TF e E3.

<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/Fig_genoma.png" alt="Estrutura genômica do CHIKV" width="90%" />
<p class="caption">
Estrutura genômica do CHIKV
</p>

Apesar de ser caracterizado como um vírus que possui 2 ORFs e o genoma de referência, disponível no NCBI ([NC\_004162](https://www.ncbi.nlm.nih.gov/nuccore/NC_004162.2)), corroborar com a afirmação, grande parte dos genomas completos de CHIKV disponíveis no NCBI apresentam um códon de parada entre os genes codificadores das proteínas nsP3 e nsP4. Inclusive, o genoma do vírus isolados de pessoas diferentes em Natal, também mostra que há o códon de parada nos isolados daqui.

<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/opal_genome.png" alt="Códon de parada (destacado pela linha horizontal vermelha) em 3 isolados do Brasil" width="90%" />
<p class="caption">
Códon de parada (destacado pela linha horizontal vermelha) em 3 isolados do Brasil
</p>

A partir desse alinhamento, tive a ideia de observar a frequência do códon de parada em todos os genomas completos de CHIKV já publicados.

Obtenção das sequências completas de CHIKV
------------------------------------------

### NCBI

1.  O primeiro passa para observar a frequência do códom de parada nos genomas de CHIVK, é baixar os genomas no [NCBI](https://www.ncbi.nlm.nih.gov/);
2.  Em seguida, nós pesquisamos por "Chikungunya" no banco de dados de "nucleotide";
3.  selecionamos o "Chikungunya virus" no campo "Results by taxon";

    3.1. O campo de busca vai ficar dessa forma: "(chikungunya) AND "Chikungunya virus"\[porgn:\_\_txid37124\]";
4.  O passo seguinte é selecionar apenas as sequências que tenham pelo menos 10.5kb de comprimento. Eu selecionei 13kb como comprimento máximo, mas isso foi apenas por causa do NCBI exigir que esse campo seja preenchido.

    <img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/ncbi_length.png" width="90%" style="display: block; margin: auto;" />

5.  por fim, vamos baixar o arquivo GenBank selecionando a opção Send to &gt; Format &gt; GenBank (full) e clicando em "Create File".

Ao final dessa etapa, nós teremos um arquivo .gb com 958 entradas. O arquivo desse trabalho foi baixado no dia 1/03/2020, é provável que esse número aumente no decorrer do tempo.

### Filtrando o arguivo

depois baixar o arguivo .gb, eu realizei uma formatação do dados, transformando o arguivo .gb em um arquivo tabular com algumas informações que eu achei necessárias para o estudo. Para formatar essa tabela, eu utilizei um código em Python (v.3.7):

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

**OBS:** *inicialmente eu tive alguns erros durante o processo de formatação, no passo do loop2, por causa de um erro nas informações da sequência MK518395. por algum motivo, quando o loop chegava nela dava erro. inicialmente eu removi entrada e rodei o script, aí deu certo. mas em seguida, eu readicionei a entrada e saí removendo cada uma das features até identificar que o que estava causando a interrupção na execução do scritp estava na informação Assembly-Data.*
