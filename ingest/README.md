# ingest

Feel free to modify

```
conda activate nextstrain # or some env with augur

bash ingest_vipr.sh

ls -ltr vipr_download.fasta

python3 format_downloaded_genomes.py \
  --in_fasta vipr_download.fasta \
  --out_fasta rsv.fasta
  
ls -ltr rsv.fasta

augur parse --sequences rsv.fasta \
      --output-sequences sequences.fasta \
      --output-metadata metadata.tsv \
      --fields strain strain_name segment date host country subtype virus

ls -ltr sequences.fasta metadata.tsv
#>   51M May 18 13:01 sequences.fasta
#>  306K May 18 13:01 metadata.tsv

grep -c ">" sequences.fasta
#> 3550
```

* [ ] Sanity check, is there supposed to be 3550 sequences? Or does REST have a default maximum.
