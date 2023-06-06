#!/bin/sh

samtools view Primary.bam | python bam_to_ml_csv.py 1> ml_stats.csv 2> ml_stats.err

