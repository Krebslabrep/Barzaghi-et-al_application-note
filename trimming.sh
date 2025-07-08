java -jar ./Trimmomatic-0.36/trimmomatic-0.36.jar \
  PE \
  -threads 4 \
  forward.fq.gz reverse.fq.gz \
  forward_paired.fq.gz \
  forward_unpaired.fq.gz \
  reverse_paired.fq.gz \
  reverse_unpaired.fq.gz \
  ILLUMINACLIP:./Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 \
  TRAILING:3 \
  SLIDINGWINDOW:4:18