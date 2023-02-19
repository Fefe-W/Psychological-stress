##Converting bams to juncs
for bamfile in `ls stress/GSE89692/vHIP/CSDS/*.bam`; do
    echo Converting $bamfile to $bamfile.junc
    samtools index $bamfile
    regtools junctions extract -a 8 -m 50 -M 500000 $bamfile -o $bamfile.junc
    echo $bamfile.junc >> vHIP_CSDS_juncfiles.txt
done

##Intron clustering
python ../clustering/leafcutter_cluster.py -j vHIP_CSDS_juncfiles.txt -m 50 -o vHIP_CSDS -l 500000

##psi quantification (used to select covariates)
../scripts/leafcutter_quantify_psi.R vHIP_CSDS_perind_numers.counts -o vHIP_CSDS_raw_psi.txt.gz -p 8

##Differential intron excision analysis
../scripts/leafcutter_ds.R vHIP_CSDS_perind_numers.counts.gz groups_file.txt -o vHIP_CSDS -e ../leafcutter/data/mouse.txt.gz -p 8

