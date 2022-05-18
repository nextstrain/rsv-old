#! /usr/bin/env python
# FROM: https://github.com/blab/rsv_adaptive_evo/blob/master/processing_scripts/format_downloaded_genomes.ipynb

# === pkgs
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys
import argparse

# === input
def parse_args():
    parser = argparse.ArgumentParser(
        description = "Takes an vipr downloaded RSV fasta, formats dates in the 4th header column"
    )
    parser.add_argument(
        "--in_fasta",
        help = "The vipr downloaded RSV fasta.",
        required = True
    )
    parser.add_argument(
        "--out_fasta",
        help = "The expected output fasta",
        required = False
    )
    return parser.parse_args()

# === main
def fix_rsv_dates(data_file, out_file, min_year):
  with open(data_file, 'r') as handle:
      
      edited_records = []
      
      for virus in SeqIO.parse(handle, 'fasta'):
          date = virus.id.split('|')[3]
          # exclude sequences with no date
          if date!= 'NA' and date != 'May_2016/Dec_2017':
              if int(date[0:4])<min_year:
                  min_year = int(date[0:4])
              # if date only has year, add -XX-XX for month and day
              if len(date)==4:
                  formatted_date = date+'-XX-XX'
              else:
                  formatted_date = date.replace('_', '-')
                  # if date only has month, add -XX for day
                  if len(formatted_date)!=10:
                      formatted_date = formatted_date+'-XX'
              virus.id = virus.id.replace(date, formatted_date)
              virus.description = virus.id
                      
              edited_records.append(SeqRecord(virus.seq, id=virus.id, description=virus.description))
              
      SeqIO.write(edited_records, out_file, "fasta")

def main():
  args = parse_args()
  data_file = args.in_fasta # 'data/vipr_download.fasta'
  out_file = args.out_fasta # 'data/rsv.fasta'
  min_year = 2000
  fix_rsv_dates(data_file, out_file, min_year)

if __name__ == "__main__":
  main()
