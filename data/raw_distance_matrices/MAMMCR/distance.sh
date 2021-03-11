for f in /n/home04/vlaine/combined/combined/09Gblocks_SNPs/06problemregions/02reorder/*.fasta; do
 T=${f##*/}
 T=${T%.fasta}_dist.tab
snp-dists "$f"  -b > "$T"
done
