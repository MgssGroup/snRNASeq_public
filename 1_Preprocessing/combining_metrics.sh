#!/bin/bash

cat female_counts/*_count/_cmdline > female_cmdlines_May2021.txt
awk 'NR == 1 {print "Filenames," $0 }' female_counts/F1_count/outs/metrics_summary.csv > female_metrics_May2021.csv
awk 'NR % 2 ==0 { print FILENAME "," $0 }' female_counts/*_count/outs/metrics_summary.csv >> female_metrics_May2021.csv
cd female_counts
for f in *_count
do
    cp "$f"/outs/web_summary.html  ../web_summaries/"${f%_*}"_web_summary.html
done
cd ..

cat male_counts/*_count/_cmdline > male_cmdlines_May2021.txt
awk 'NR == 1 {print "Filenames," $0 }' male_counts/M1_count/outs/metrics_summary.csv > male_metrics_May2021.csv
awk 'NR % 2 ==0 { print FILENAME "," $0 }' male_counts/*_count/outs/metrics_summary.csv >> male_metrics_May2021.csv
cd male_counts
for f in *_count
do
    cp "$f"/outs/web_summary.html  ../web_summaries/"${f%_*}"_web_summary.html
done
cd ..
zip -r web_summaries.zip web_summaries
