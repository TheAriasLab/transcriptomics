# Transcriptomics resources for wild and cultivate cacao: Finding genes responsible for domestication

A continuación, se describe el análisis de datos realizado para la identificación de genes responsables en la domesticación de Cacao en comparación con cultivos domésticos y salvajes. 

Los datos se corrieron un clúster mediante la administración de código a través de SLURM. 

##	Descarga de las secuencias
Las secuencias pueden ser descargadas desde la librería de SRA del NCBI mediante la herramienta sratoolkit https://github.com/ncbi/sra-tools  
Las secuencias usadas fueron las siguientes; 

**Cultivado**
SRR3217276, SRR3217277, SRR3217298, SRR3217299, SRR3217301, SRR3217304, SRR3217315, SRR3217317, SRR3217318, SRR3217319

**Silvestre**
SRR3217278, SRR3217279, SRR3217280, SRR3217281, SRR3217282, SRR3217283, SRR3217284, SRR3217292, SRR3217294, SRR3217297

```
module load software/bioinformatics/sratoolkit/2.8.1
prefetch SRR3217276
```
Los archivos descargados por medio del SRA deben ser separadas en archivos pair-end, además este archivo (el header especificamente) debe ser compatible con procedimientos posteriores como Trinity. Por tal usamos la función fastq-dump perteneciente también al SRA tools del NCBI. 

```
module load software/bioinformatics/sratoolkit/2.8.1
fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files *.sra -o /BIOS-Share/home/areyesb/Secuencias
#/BIOS-Share/home/areyesb/Secuencias hace referencia a la carpeta donde están todas las secuencias con extensión .sra
#Se usa *.sra, para que el procedimiento se haga a todas las secuencias. 
```

##	Calidad de las secuencias
Luego del procedimiento con fastq-dump, las secuencias quedan con formato fastq, así que usamos la herramienta fastqc (http://www.bioinformatics.babraham.ac.uk/projects/fastqc/ ) para revisar la calidad de las secuencias. 

```
module load software/bioinformatics/fastqc/0.11.4
fastqc *.fastq
```

Posteriormente al analizar la calidad de las secuencias, se decidió elegir las siguientes características para todas las secuencias, para la limpieza de la mala calidad proveniente del proceso de secuenciación. 
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3
Para esta limpieza se usó la herramienta Trimmomatic (http://www.usadellab.org/cms/index.php?page=trimmomatic ) que permite realizar todas las estas ediciones directamente. 

```
module load software/bioinformatics/trimmomatic/0.36
TrimmomaticPE /BIOS-Share/home/areyesb/Secuencias/SRR3217278_1.fastq /BIOS-Share/home/areyesb/Secuencias/SRR3217278_2.fastq SRR3217278_1_Paired.fastq SRR3217278_1_UNPaired.fastq SRR3217278_2_Paired.fastq SRR3217278_2_UNPaired.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:28 HEADCROP:10 MINLEN:3
```
Las secuencias se volvieron a revisar usando la herramienta de fastqc, para verificar que el proceso de Trimmomatic haya resultado exitoso. 

##	Descontaminación de las lecturas

Se decide descontaminar las lecturas para eliminar datos provenientes de otros organismos, principalmente de bacterias, que puedan interferir en el ensamblaje o en la creación de transcritos incorrectos. 
Para esta tarea se decide usar la herramienta de kaiju (http://kaiju.binf.ku.dk/ ). Primero se debe crear la base de datos de kaiju dentro del clúster debido a la cantidad de datos que se deben analizar referentes a las lecturas. En este caso se usó la base de datos más completa para realizar la mejor clasificación taxonómica posible, nr_euk (Subset of NCBI BLAST nr database containing all proteins belonging to Archaea, Bacteria and Viruses, and additionally include proteins from fungi and microbial eukaryotes). 

```
module load devtools/bioconda/bioconda3
source activate kaiju_env
kaiju-makedb -s nr_euk
```

Posteriormente, se realiza la clasificación taxonómica de kaiju, el cual tiene la siguiente sintaxis; 
kaiju -t nodes.dmp -f kaiju_db_*.fmi -i firstfile.fastq and -j secondfile.fastq
Esto se realiza para todas las secuencias. 

```
module load devtools/bioconda/bioconda3
source activate kaiju_env
kaiju -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -f /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nr_euk/kaiju_db_nr_euk.fmi -i /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217276_1_Paired.fastq -j /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217276_2_Paired.fastq -o kaiju_cultivado1.out
```

Luego se debe clasificar los datos en crudo de kaiju, asignando una taxonomía a cada una de las lecturas mapeadas, para este caso se eligió usar la taxonomía al nivel de -superkingdom, para clasificar de una manera mas adecuada y grupal las lecturas que mapearían con bacterias.
Esto se realiza para todos los outputs del proceso anterior de kaiju.

```
module load devtools/bioconda/bioconda3
source activate kaiju_env
kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_cultivado_SRR3217276.out -o kaiju_names_cult_superkingdom1.out
```

Debido a que se deben extraer las lecturas que se van a usar en los análisis posteriores, se decide hacer una extracción de estas lecturas de los outputs de kaiju usando comandos externos. Esto se realiza para todos los outputs de kaiju.

**a)**	Se extraen las lecturas que no mapean con ninguna taxonomía (Lo que influye en que estas lecturas son presentan ninguna contaminación, para esto se usa el output en bruto del primer proceso de kaiju).

awk '/U/ {print $2"/1"}' kaiju_silvestre_SRR3217278.out > lista_Uclasi_silv_SRR3217278_1.txt
awk '/U/ {print $2"/2"}' kaiju_silvestre_SRR3217278.out > lista_Uclasi_silv_SRR3217278_2.txt

**b)**	Se extraen las lecturas que mapean como Eukaryota (Se usa el archivo output de la taxonomía de kaiju).
awk '/Eukaryota/ {print $2"/1"}' kaiju_names_silv_superkingdom1.out > lista_silv_Eukaryota_SRR3217278_1.txt
awk '/Eukaryota/ {print $2"/1"}' kaiju_names_silv_superkingdom1.out > lista_silv_Eukaryota_SRR3217278_1.txt

**c)**	Se unen los dos archivos generados en el paso a y b, para generar una única lista de lecturas a extraer a las lecturas, y que serán usadas en los análisis posteriores. 

```
cat /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/lista_Uclasi_silv_SRR3217278_1.txt /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju-clasificado-eukaryota/lista_silv_Eukaryota_SRR3217278_1.txt > /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas-eukaryota/Lista_filtro_silv-Eukaryota_SRR3217278_1.txt
```
_Nota: Hay que recordar que cada uno de estos comandos debe correrse por cada una de las secuencias, pero además de eso, cada comando debe generar una lista duplicada, debido a que se necesita una lista por cada archivo pair-end, en este caso cada archivo pair-end tiene un header 2/1 y 2/2, haciendo referencia al pair-end de cada lectura, por tal cada uno debe ser extraído en listas diferentes para poder realizar su extracción correctamente. 

Las lecturas deben ser extraídas con la lista de lecturas que acabamos de generar, para esto usamos la herramienta seqtk (https://github.com/lh3/seqtk ), que nos permite a través de una lista que contenga el header de la lectura extraer en un nuevo archivo fastq. 

```
module load software/bioinformatics/seqtk/1.2
seqtk subseq /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado/SRR3217276_1_Paired.fastq /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/Archivos-listas/Lista_filtro_culti_SRR3217276_1.txt > /BIOS-Share/home/areyesb/TrimmomaticResult/Cultivado_filtrado/SRR3217276_1_Paired_filtrado.fastq
```

Este último comando ya nos permitirá tener en una carpeta todas nuestras lecturas ya descontaminadas.



##	Ensamblaje de las lecturas usando Trinity
Para el ensamblaje de las lecturas se usó Trinity, el cual es un ensamblador transcriptómico de novo de fácil uso (https://github.com/trinityrnaseq/trinityrnaseq/wiki ).

```
module load software/bioinformatics/trinityrnaseq/2.8.4
Trinity --seqType fq --max_memory 150G --left /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217278_1_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217279_1_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217280_1_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217281_1_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217282_1_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217283_1_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217284_1_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217292_1_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217294_1_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217297_1_Paired.fastq --right /BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217278_2_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217279_2_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217280_2_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217281_2_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217282_2_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217283_2_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217284_2_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217292_2_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217294_2_Paired.fastq,/BIOS-Share/home/areyesb/TrimmomaticResult/Silvestre/SRR3217297_2_Paired.fastq --CPU 32
```


**Estadisticas basicas de Trinity**
Algunas estadisticas basicas del ensamblaje de Trinity puede obtenerse fácilmente, así aseguramos que el ensamblaje haya sido correcto. 

```
module load software/bioinformatics/trinityrnaseq/2.8.4
TrinityStats.pl Trinity.fasta
```

Este nos dará información básica del ensamblaje, como el siguiente ejemplo; 

```
Counts of transcripts, etc.

Total trinity 'genes':  63526
Total trinity transcripts:      96620
Percent GC: 40.75

Stats based on ALL transcript contigs:

        Contig N10: 3937
        Contig N20: 3133
        Contig N30: 2617
        Contig N40: 2240
        Contig N50: 1918
        Median contig length: 599.5
        Average contig: 1073.57
        Total assembled bases: 103728470

Stats based on ONLY LONGEST ISOFORM per 'GENE':

        Contig N10: 3742
        Contig N20: 2933
        Contig N30: 2432
        Contig N40: 2050
        Contig N50: 1682
        Median contig length: 400
        Average contig: 833.70
        Total assembled bases: 52961586
```


## Validación del ensamblaje
Una forma de validar las lecturas es usando herramientas que nos permitan calcular si las lecturas son suficientes pata nuestro análisis, como, por ejemplo, si el muestreo de secuenciación fue suficiente o si las lecturas tienen los genes suficientes para nuestros análisis posteriores. 
Para este caso usaremos busco, que nos permite hallar el porcentaje de genes ortólogos que se encuentran en nuestras lecturas con base en un set de datos más cercano filogenéticamente a nuestro organismo de estudio, para este caso se usó la base de datos embryophyta.

```
module load devtools/bioconda/bioconda3
source activate busco_env
busco -i /BIOS-Share/home/areyesb/kaiju/Filtrar_secuencias/Filter_script/TrinityCultivadofiltrado.fasta -o busco.out -l embryophyta_odb10 -m tran -c 32 -f
```

Un resultado de busco aparece como el siguiente;

```
  Results
	
  C:81.0%[S:42.9%,D:38.1%],F:11.1%,M:7.9%,n:1614	   
	1308	Complete BUSCOs (C)			   
	693	Complete and single-copy BUSCOs (S)	   
	615	Complete and duplicated BUSCOs (D)	   
	179	Fragmented BUSCOs (F)			   
	127	Missing BUSCOs (M)			   
	1614	Total BUSCO groups searched
```



## Anotación de novo del ensamblaje
Después de realizar la evaluación de calidad del ensamblaje se prosigue con la anotación del ensamblaje, usualmente la anotación de novo se realiza con el objetivo de encontrar la mayor cantidad de genes en genomas desconocidos, en este caso se decidió hacer una anotación de novo para el transcriptoma del cacao silvestre, y una anotación normal para el transcriptoma cultivado.

Para esto se deicidio usar la herramienta Trinotate (https://github.com/Trinotate/Trinotate.github.io/wiki), la cual involucra una serie de pasos para mapear y anotar todos los genes posibles por medio de diferentes herramientas.

Instalación de las bases de datos requeridas
Se debe instalar una base de datos denominada Trinotate.sqlite, que será la encargada de recopilar todos los datos de anotación del transcriptoma. Esta base de datos se eligió descargarse una ubicación diferente al del programa. 

```
module load devtools/bioconda/bioconda3
source activate trinotate_env
Build_Trinotate_Boilerplate_SQLite_db.pl /BIOS-Share/home/areyesb/Anotacion/Trinotate-db
```

Luego se requiere descargar dos bases de datos, una perteneciente a la base de datos básica de uni-prot, para reconocer proteínas, y otra denominada Pfam.

```
module load devtools/bioconda/bioconda3
source activate trinotate_env
makeblastdb -in uniprot_sprot.pep -dbtype prot
```

Descomprimir la base de datos de Pfam. 
```
gunzip Pfam-A.hmm.gz
```

### Uso de Transdecoder
Se requiere el uso de Transdecoder (https://github.com/TransDecoder/TransDecoder/wiki), el cual requiere el ensamblaje y un archivo denominado Trinity.fasta.transdecoder.pep, para generar este archivo se deben seguir los siguientes pasos; 

**Paso 1: Extraer marcos de lectura abiertos**
```
module load devtools/bioconda/bioconda3
source activate transdecoder_env
TransDecoder.LongOrfs -t /BIOS-Share/home/areyesb/Trinity/Cultivado-Ensam/trinity_out_dir/TrinityCultivado.fasta
```

**Paso 2. Incluir búsquedas de homología como criterios de retención de ORF**
Aquí se incluye una base de datos personalizada de swissprot, se incluyó la de uniprot de cacao (https://www.uniprot.org/uniprot/?query=theobroma%20cacao&fil=proteome%3AUP000026915+AND+organism%3A%22Theobroma+cacao+%28Cacao%29+%28Cocoa%29+%5B3641%5D%22&sort=score ).
```
module load software/bioinformatics/ncbi-blast/2.7.1
blastp -query /BIOS-Share/home/areyesb/Anotacion/Transdecoder/TrinitySilvestre.fasta.transdecoder_dir/longest_orfs.pep -db /BIOS-Share/home/areyesb/ProteCacao_Uniprot/uniprot-theobroma_cacao-filtered.fasta  -max_target_seqs 1 -outfmt 6 -evalue 1e-5 -num_threads 10 > blastp.outfmt6
```
_Usar la base de datos Pfam_
```
module load devtools/bioconda/bioconda3
source activate hmmer_env
hmmpress /BIOS-Share/home/areyesb/Anotacion/Transdecoder/Pfam-db/Pfam-A.hmm
```
_Creacion archive .pep_
```
module load devtools/bioconda/bioconda3
source activate hmmer_env
hmmscan --cpu 32 --domtblout pfam.domtblout /BIOS-Share/home/areyesb/Anotacion/Transdecoder/Pfam-db/Pfam-A.hmm /BIOS-Share/home/areyesb/Anotacion/Transdecoder/TrinitySilvestre.fasta.transdecoder_dir/longest_orfs.pep
```
**Paso 3. Integrar los resultados de búsqueda de Blast y Pfam en la selección de la región de codificación**
```
module load devtools/bioconda/bioconda3
source activate transdecoder_env
TransDecoder.Predict -t /BIOS-Share/home/areyesb/Trinity/Silvestre-Ensam/trinity_out_dir/TrinitySilvestre.fasta --retain_pfam_hits /BIOS-Share/home/areyesb/Anotacion/Transdecoder/pfam.domtblout --retain_blastp_hits /BIOS-Share/home/areyesb/Anotacion/Transdecoder/blastp.outfmt6
```

**Paso 4. Capturar homologías BLAST**
En este paso se involucran dieversas bases de datos, e incouso bases de datos personalizadas. 
* Blastx busca transcripciones de Trinity
* Buscar proteínas predichas por transdecodificador
* Inclusión de base de datos de referencia general de proteínas.
```
module load software/bioinformatics/ncbi-blast/2.7.1
blastx -query /BIOS-Share/home/areyesb/Trinity/Silvestre-Ensam/trinity_out_dir/TrinitySilvestre.fasta -db /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/uniprot_sprot.pep -num_threads 32 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx1.outfmt6
```
```
module load software/bioinformatics/ncbi-blast/2.7.1
blastp -query /BIOS-Share/home/areyesb/Anotacion/Transdecoder/TrinitySilvestre.fasta.transdecoder.pep -db /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/uniprot_sprot.pep -num_threads 32 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp2.outfmt6
```

**Inclusión de base de datos de referencia de cacao**
```
module load software/bioinformatics/ncbi-blast/2.7.1
blastx -query /BIOS-Share/home/areyesb/Trinity/Silvestre-Ensam/trinity_out_dir/TrinitySilvestre.fasta -db /BIOS-Share/home/areyesb/ProteCacao_Uniprot/uniprot-theobroma_cacao-filtered.fasta -num_threads 32 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastx3.outfmt6
```
```
module load software/bioinformatics/ncbi-blast/2.7.1
blastp -query /BIOS-Share/home/areyesb/Anotacion/Transdecoder/TrinitySilvestre.fasta.transdecoder.pep -db /BIOS-Share/home/areyesb/ProteCacao_Uniprot/uniprot-theobroma_cacao-filtered.fasta -num_threads 32 -max_target_seqs 1 -outfmt 6 -evalue 1e-3 > blastp4.outfmt6
```

**Ejecución de HMMER para identificar dominios de proteínas**
```
module load devtools/bioconda/bioconda3
source activate hmmer_env
hmmscan --cpu 32 --domtblout /BIOS-Share/home/areyesb/Anotacion/Transdecoder/pfam.domtblout /BIOS-Share/home/areyesb/Anotacion/Transdecoder/Pfam-db/Pfam-A.hmm /BIOS-Share/home/areyesb/Anotacion/Transdecoder/TrinitySilvestre.fasta.transdecoder.pep > pfam.log
```

**Ejecución de signalP para predecir péptidos señal**
```
module load software/bioinformatics/signalp/4.1
signalp -f short -n signalp.out /BIOS-Share/home/areyesb/Anotacion/Transdecoder/TrinitySilvestre.fasta.transdecoder.pep
```

**Ejecución de tmhmm para predecir regiones transmembrana**
```
module load software/bioinformatics/tmhmm/2.0c
tmhmm --short < /BIOS-Share/home/areyesb/Anotacion/Transdecoder/TrinitySilvestre.fasta.transdecoder.pep > tmhmm.out
```

**Ejecución de tmhmm para predecir regiones transmembrana**```

```
module load software/bioinformatics/tmhmm/2.0c
tmhmm --short < /BIOS-Share/home/areyesb/Anotacion/Transdecoder/TrinitySilvestre.fasta.transdecoder.pep > tmhmm.out
```

**Ejecución de RNAMMER para identificar las transcripciones de ARNr**
```
/BIOS-Share/Software/devtools/bioconda/miniconda3/python3.7/envs/trinotate_env/bin/RnammerTranscriptome.pl --transcriptome /BIOS-Share/home/areyesb/Trinity/Silvestre-Ensam/trinity_out_dir/TrinitySilvestre.fasta --path_to_rnammer /BIOS-Share/Software/bioinformatics/rnammer/1.2
```

### Cargar los resultados generados en una base de datos Trinotate SQLite

**1. Cargar transcripciones y regiones de codificación**
```
module load software/bioinformatics/trinityrnaseq/2.8.4
/BIOS-Share/Software/bioinformatics/trinityrnaseq/2.8.4/util/support_scripts/get_Trinity_gene_to_trans_map.pl /BIOS-Share/home/areyesb/Trinity/Silvestre-Ensam/trinity_out_dir/TrinitySilvestre.fasta > Trinity.fasta.Silv.gene_trans_map
```

_Cargar la información en la base de datos Trinotate sqlite_
```
module load devtools/bioconda/bioconda3
source activate trinotate_env
Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite init --gene_trans_map /BIOS-Share/home/areyesb/Anotacion/sqlite-db/Trinity.fasta.Silv.gene_trans_map --transcript_fasta /BIOS-Share/home/areyesb/Trinity/Silvestre-Ensam/trinity_out_dir/TrinitySilvestre.fasta --transdecoder_pep /BIOS-Share/home/areyesb/Anotacion/Transdecoder/TrinitySilvestre.fasta.transdecoder.pep
```

**2. Carga de homologías BLAST**
```
module load devtools/bioconda/bioconda3
source activate trinotate_env
Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite LOAD_swissprot_blastp /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Blastp/blastp2.outfmt6
```
```
module load devtools/bioconda/bioconda3
source activate trinotate_env
Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite LOAD_swissprot_blastx /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Blastp/blastx1.outfmt6
```

_Cargar bases de datos blast de referencia personalizadas_
_Se incluyo la base de datos de referencia de cacao_
```
module load devtools/bioconda/bioconda3
source activate trinotate_env
Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite LOAD_custom_blast --outfmt6 /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Blastp/blastx3.outfmt6 --prog blastx --dbtype blastx_cacao
```
```
module load devtools/bioconda/bioconda3
source activate trinotate_env
Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite LOAD_custom_blast --outfmt6 /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Blastp/blastp4.outfmt6 --prog blastp --dbtype blastp_cacao
```

**3. Cargar entradas de dominio Pfam**
```
module load devtools/bioconda/bioconda3
source activate trinotate_env
Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite LOAD_pfam /BIOS-Share/home/areyesb/Anotacion/Transdecoder/pfam.log
```

**4. Cargar dominios transmembrane**
```
module load devtools/bioconda/bioconda3
source activate trinotate_env
Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite LOAD_tmhmm /BIOS-Share/home/areyesb/Anotacion/tmhmm/tmhmm.out
```

**5. Cargar predicciones de péptidos señal**
```
module load devtools/bioconda/bioconda3
source activate trinotate_env
Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite LOAD_signalp /BIOS-Share/home/areyesb/Anotacion/signal/signalp.out
```

## Trinotate: Generar un informe de anotación
```
module load devtools/bioconda/bioconda3
source activate trinotate_env
Trinotate /BIOS-Share/home/areyesb/Anotacion/Trinotate-db/Trinotate.sqlite report [opts] > trinotate_annotation_report.xls
```
