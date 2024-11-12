#!/bin/bash

awk '$1 != "AnalysisLevel"' $1 | awk '$3 != "readsLongEnough"' | awk '$3 != "totalReads"' | sed 's| |_|g' | awk '{print $3,$6,"metamaps"}' | sed 's| |\t|g' | sed 's|_| |g' | sed 's|Unclassified|unclassified|g'