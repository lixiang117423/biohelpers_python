hisat2 -x example/acuce.genome.hisat2.index -1 example/R24129662_1.fq.gz -2 example/R24129662_2.fq.gz -p 60 | samtools sort -@ 60 -O BAM -o mapping/R24129662.sorted.bam -
hisat2 -x example/acuce.genome.hisat2.index -1 example/R24129663_1.fq.gz -2 example/R24129663_2.fq.gz -p 60 | samtools sort -@ 60 -O BAM -o mapping/R24129663.sorted.bam -
hisat2 -x example/acuce.genome.hisat2.index -1 example/R24129662_1.fq.gz -2 example/R24129662_2.fq.gz -p 60 | samtools sort -@ 60 -O BAM -o mapping/R24129662.sorted.bam -
hisat2 -x example/acuce.genome.hisat2.index -1 example/R24129663_1.fq.gz -2 example/R24129663_2.fq.gz -p 60 | samtools sort -@ 60 -O BAM -o mapping/R24129663.sorted.bam -
