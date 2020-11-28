#!/bin/bash
#SBATCH -J kaiju_superkindom_cult
#SBATCH -o kaiju_superkindom_cult.out
#SBATCH -e kaiju_superkindom_cult.err
#SBATCH -n 32
#SBATCH --partition=normal

module load devtools/bioconda/bioconda3
source activate kaiju_env

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_cultivado_SRR3217276.out -o kaiju_names_cult_superkingdom1.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_cultivado_SRR3217277.out -o kaiju_names_cult_superkingdom2.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_cultivado_SRR3217298.out -o kaiju_names_cult_superkingdom3.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_cultivado_SRR3217299.out -o kaiju_names_cult_superkingdom4.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_cultivado_SRR3217301.out -o kaiju_names_cult_superkingdom5.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_cultivado_SRR3217304.out -o kaiju_names_cult_superkingdom6.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_cultivado_SRR3217315.out -o kaiju_names_cult_superkingdom7.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_cultivado_SRR3217317.out -o kaiju_names_cult_superkingdom8.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_cultivado_SRR3217318.out -o kaiju_names_cult_superkingdom9.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_cultivado_SRR3217319.out -o kaiju_names_cult_superkingdom10.out