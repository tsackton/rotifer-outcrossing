#code to extract haplotypes for mapping

bioawk -c fastx '{print FILENAME,$name,$seq}' */*.fasta > regions.tab

#fiter to keep only those in the final list, removing overlaps
#need to fix filters to add N to line

sed -i 's/$/N/' MACR_filter.txt
sed -i 's/^/\//' MACR_filter.txt
grep "^MA_CR" regions.tab | grep -F -f MACR_filter.txt > MACR_regions.tab

sed -i 's/$/N/' MAMM_filter.txt
sed -i 's/^/\//' MAMM_filter.txt
grep "^MA_MM" regions.tab | grep -F -f MAMM_filter.txt > MAMM_regions.tab

sed -i 's/$/N/' MMCR_filter.txt
sed -i 's/^/\//' MMCR_filter.txt
grep "^MM_CR" regions.tab | grep -F -f MMCR_filter.txt > MMCR_regions.tab

sed -i 's/$/N/' three_way_filter.txt
sed -i 's/^/\//' three_way_filter.txt
grep "^three-way" regions.tab | grep -F -f three_way_filter.txt > three_way_regions.tab

#merge
cat MACR_regions.tab MMCR_regions.tab MAMM_regions.tab three_way_regions.tab > nr_regions.tab

#decide to map to CR1, MA1, MM1
#FILENAME
#MA_CR/100Nsremoved_gb_rename_nospace_srt_1kbtrim.fasta

grep "CR1" nr_regions.tab | perl -p  -e 's/^(\w+[-_]\w+)\/(\d+)\S+/>$1.$2/' | awk '{print $1"_"$2; print $3}' > CR_regions.fasta
grep "MA1" nr_regions.tab | perl -p  -e 's/^(\w+[-_]\w+)\/(\d+)\S+/>$1.$2/' | awk '{print $1"_"$2; print $3}' > MA_regions.fasta
grep "MM1" nr_regions.tab | perl -p  -e 's/^(\w+[-_]\w+)\/(\d+)\S+/>$1.$2/' | awk '{print $1"_"$2; print $3}' > MM_regions.fasta

####### mapping -make index ####

for FILE in $(ls *.fasta); do bwa index $FILE; done
for FILE in $(ls *.fa); do bwa index $FILE; done

#### mapping - run mapping, use SLURM scripts ##
