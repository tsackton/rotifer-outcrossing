for f in /n/home04/vlaine/combined/combined/09Gblocks_SNPs/06problemregions/*.fasta; do
 T=${f##*/}
 T=${T%.fasta}_srt.fasta
seqkit sort "$f"  > "$T"
done
