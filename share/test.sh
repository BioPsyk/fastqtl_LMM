#!/bin/bash

if [[ $( ./bin/alfonso --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 100000,4 --Sinv data/Sinv --eigen data/eigen  | sha1sum | awk {'print toupper($1)'}) == "80745CCF282B5069B75B99470F26DB8466B15ABE" ]]; then
    echo "Passed: based on Sinv and eigen."
else
    echo "Failed: based on Sinv and eigen."
    exit 1
fi

if [[ $( ./bin/alfonso --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 100000,4 --kinship data/kinship  | sha1sum | awk {'print toupper($1)'}) == "0E50C4075C347F5E596EE120A68B7154C64625C2" ]]; then
    echo "Passed: based on kinship matrix."
else
    echo "Failed: based on kinship matrix."
    exit 1
fi

if [[ $( ./bin/alfonso --bed data/phenotype.bed --job-number 1 --genes 10 --vcf data/genotype.vcf.gz --perm 100000,4 --kinship data/kinship --calc-h2 | sha1sum | awk {'print toupper($1)'}) == "F2C6522C50D9084D26ABD865C8AA3285DF9A1304" ]]; then
    echo "Passed: test with heritability calculation."
else
    echo "Failed: test with heritability calculation."
    exit 1
fi
