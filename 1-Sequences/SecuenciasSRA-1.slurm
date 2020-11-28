#!/bin/bash
#SBATCH -J sec-download
#SBATCH -o sec-download.log
#SBATCH -e sec-download.log
#SBATCH -n 16
#SBATCH --partition=normal

module load software/bioinformatics/sratoolkit/2.8.1

#descargar secuencias

module load software/bioinformatics/sratoolkit/2.8.1

#Cultivado
prefetch SRR3217276 SRR3217277 SRR3217298 SRR3217299 SRR3217301 SRR3217304 SRR3217315 SRR3217317 SRR3217318 SRR3217319

#Silvestre
prefetch SRR3217278 SRR3217279 SRR3217280 SRR3217281 SRR3217282 SRR3217283 SRR3217284 SRR3217292 SRR3217294 SRR3217297

