#!/bin/bash
#SBATCH -J kaiju_eukaryota_silv
#SBATCH -o kaiju_eukaryota_silv.out
#SBATCH -e kaiju_eukaryota_silv.err
#SBATCH -n 32
#SBATCH --partition=normal

module load devtools/bioconda/bioconda3
source activate kaiju_env

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_silvestre_SRR3217278.out -o kaiju_names_silv_superkingdom1.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_silvestre_SRR3217279.out -o kaiju_names_silv_superkingdom2.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_silvestre_SRR3217280.out -o kaiju_names_silv_superkingdom3.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_silvestre_SRR3217281.out -o kaiju_names_silv_superkingdom4.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_silvestre_SRR3217282.out -o kaiju_names_silv_superkingdom5.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_silvestre_SRR3217283.out -o kaiju_names_silv_superkingdom6.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_silvestre_SRR3217284.out -o kaiju_names_silv_superkingdom7.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_silvestre_SRR3217292.out -o kaiju_names_silv_superkingdom8.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_silvestre_SRR3217294.out -o kaiju_names_silv_superkingdom9.out

kaiju-addTaxonNames -t /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/nodes.dmp -n /BIOS-Share/home/areyesb/kaiju/Kaiju-NR_EUK/names.dmp -r superkingdom -i /BIOS-Share/home/areyesb/kaiju/kaiju-secuencias/kaiju-grupos/kaiju_silvestre_SRR3217297.out -o kaiju_names_silv_superkingdom10.out