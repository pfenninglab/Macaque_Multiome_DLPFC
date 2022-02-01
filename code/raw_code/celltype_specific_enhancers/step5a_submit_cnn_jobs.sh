#!/bin/bash

ls ../../../data/tidy_data/celltype_specific_enhancers/fasta/*vs* | awk '{split($0, a, "/"); split(a[8], b, "_fold1"); split(b[1], c, "10_"); print c[2]}'  | sort | uniq | while read line; do
  echo $line
  sbatch step5b_run_cnn.sb $line
done
