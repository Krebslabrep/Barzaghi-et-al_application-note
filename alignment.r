library(QuasR)
library(BSgenome.Mmusculus.UCSC.mm10, lib.loc = "/g/krebs/barzaghi/R/x86_64-pc-linux-gnu-library/4.1/")

sample_file = "/path/to/sample/file" 
# this contains paths to the trimmed fastq pairs
# For more details on how to specify this file refer to the QuasR documentation at https://www.rdocumentation.org/packages/QuasR/versions/1.12.0/topics/qAlign

cl = makeCluster(40)
prj <- qAlign(sampleFile=sample_file,
              genome="BSgenome.Mmusculus.UCSC.mm10",
              aligner = "Rbowtie",
              paired="fr",
              bisulfite="undir",
              projectName="prj",
              alignmentParameter = "-e 70 -X 1000 -k 2 --best -strata",
              alignmentsDir="./",
              snpFile = SNPs,
              clObj = cl, 
              cacheDir = "/scratch/barzaghi/tmp")
stopCluster(cl)