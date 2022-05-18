#! /usr/bin/env bash
set -v

# Draft pull from vipr, feel free to edit, replace, put into snakemake rule, etc

# === output
output_fasta="vipr_download.fasta"

# === Split out into parameters for generalizability
FAMILY=pneumoviridae
VIRUS="Respiratory%20syncytial%20virus"
MINYEAR=2000
MINLEN=5000 #<= maybe
METADATA="genbank,strainname,segment,date,host,country,genotype,species"

URL="https://www.viprbrc.org/brc/api/sequence?datatype=genome&family=${FAMILY}&${VIRUS}&fromyear=${MINYEAR}&minlength=${MINLEN}&metadata=${METADATA}&output=fasta"
echo ${URL}
curl ${URL} \
    | tr '-' '_' \
    | tr ' ' '_' \
    | sed 's:N/A:NA:g' \
    > ${output_fasta}

# # === Fix dates
# python3 format_downloaded_genomes.py --in_fasta vipr_download.fasta --out_fasta rsv.fasta
# 
# # === Parse into sequences.fasta and metadata.tsv
# augur parse --sequences rsv.fasta \
#       --output-sequences sequences.fasta \
#       --output-metadata metadata.tsv \
#       --fields strain strain_name segment date host country subtype virus
#     
