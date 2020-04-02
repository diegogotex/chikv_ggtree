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
4.  O passo seguinte é selecionar apenas as sequências que tenham pelo menos 11kb de comprimento. Eu selecionei 13kb como comprimento máximo, mas isso foi apenas por causa do NCBI exigir que esse campo seja preenchido.

<img src="https://github.com/diegogotex/chikv_ggtree/blob/master/Figs/ncbi_length.png" alt="Códon de parada (destacado pela linha horizontal vermelha) em 3 isolados do Brasil" width="90%" />
<p class="caption">
Códon de parada (destacado pela linha horizontal vermelha) em 3 isolados do Brasil
</p>
