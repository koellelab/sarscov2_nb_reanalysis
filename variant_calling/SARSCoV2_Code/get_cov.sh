for i in ./data/seq/*/*_sorted_viral_rh_clip_viral.bam
do 
echo $i
samtools depth -d 0 -aa $i > ${i%.bam}.cov
done
