for f in /n/home04/vlaine/combined/MA-MM/MA-MM/05gblocks_trimming/trim/*.fasta; do
 T=${f##*/}
 T=${T%.fasta}_dist.tab
snp-dists "$f"  -b > "$T"
done
