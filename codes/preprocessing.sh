##Index
STAR --runMode genomeGenerate --genomeDir mouse/index --genomeFastaFiles database/mouse/GRCm38.p6.genome.fa --sjdbGTFfile database/mouse/gencode.vM23.annotation.gtf --sjdbOverhang 100

##Alignment
ls *_1_paired.fastq.gz | while read filename
do 
samplename=`basename ${filename} _1_paired.fastq.gz`
STAR --runThreadN 10 --genomeDir mouse/index --readFilesIn ${samplename}_1_paired.fastq.gz ${samplename}_2_paired.fastq.gz --outSAMstrandField intronMotif --quantMode TranscriptomeSAM --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${samplename} --outFilterType BySJout && rm *.toTranscriptome.out.bam && java -jar software/picard-tools-1.119/AddOrReplaceReadGroups.jar I=${samplename}Aligned.sortedByCoord.out.bam O=${samplename}Aligned.sortedByCoord.out_Pit.bam PL=illumine ID=${samplename}_sequence LB=${samplename}_sequence SM=Si PU=HWI-ST303 && java -jar software/picard-tools-1.119/ReorderSam.jar I=${samplename}Aligned.sortedByCoord.out_Pit.bam O=${samplename}Aligned.sortedByCoord.out_Pit_reorder.bam R=database/mouse/GRCm38.p6.genome.fa CREATE_INDEX=TRUE && java -jar software/picard-tools-1.119/MarkDuplicates.jar INPUT=${samplename}Aligned.sortedByCoord.out_Pit_reorder.bam OUTPUT=${samplename}Aligned.sortedByCoord.out_Pit_re_marked.bam REMOVE_DUPLICATES=false METRICS_FILE=${samplename}.txt ASSUME_SORTED=true && samtools index ${samplename}Aligned.sortedByCoord.out_Pit_re_marked.bam
done

##RSEM-calculate-expression
##single end
ls *.after.fastq.gz | while read filename
do 
samplename=`basename ${filename} .after.fastq.gz`
mkdir -p ${samplename}
rsem-calculate-expression --star --star-path ~/software/STAR_2.6.0a/bin/Linux_x86_64 --star-gzipped-read-file -p 8 ${samplename}.after.fastq.gz mouse/index ${samplename} 
done

##paired end
ls *_1_paired.fastq.gz | while read filename
do 
samplename=`basename ${filename} _1_paired.fastq.gz`
mkdir -p ${samplename}
rsem-calculate-expression --paired-end --star --star-path ~/software/STAR_2.6.0a/bin/Linux_x86_64 --star-gzipped-read-file -p 8 ${samplename}_1_paired.fastq.gz ${samplename}_2_paired.fastq.gz mouse/index ${samplename} 
done

##Get the expression matrix
ls *.genes.results > samplename
cat samplename |while read samplename
do 
cat ${samplename} | awk 'NR==1{print "'$samplename'"}NR>=2{print $5}'>${samplename}.count
cat ${samplename} | awk 'NR==1{print "'$samplename'"}NR>=2{print $6}'>${samplename}.tpm; 
cat ${samplename} | awk 'NR==1{print "'$samplename'"}NR>=2{print $7}'>${samplename}.fpkm; 
done

cat SRR9161643.genes.results | awk '{print $1}' > gene_ID
paste gene_ID *.count >raw_count
paste gene_ID *.tpm >raw_tpm
paste gene_ID *.fpkm >raw_fpkm

