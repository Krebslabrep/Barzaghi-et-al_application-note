#!/bin/bash

# this writes a sampleFile.txt in the current working directory using all the bam files contained in the path passed as $1 argument  

echo "FileName  SampleName" > sampleFile.txt
for file in $(ls $1/*.bam); do echo $file $(echo $file | sed 's/\//\t/g' | rev | cut -f 1 | rev | sed 's/.bam//g'); done | sed 's/ /\t/g' >> sampleFile.txt

