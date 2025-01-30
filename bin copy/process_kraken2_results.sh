#!/bin/bash

awk '$4 == "S" || $4 == "U"' $1 | sed 's|   *||g' | sed 's| |_|g' | awk '{print $6,$1,"kraken2"}' | sed 's| _|\t|g' | sed 's| |\t|g' | sed 's|_| |g' | sort -t$'\t' -k 2,2nr