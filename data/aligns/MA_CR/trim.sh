for f in /n/home04/vlaine/combined/MA-CR/MA-CR/05gblocks_trimming/reorder/*.fasta; do
 T=${f##*/}
 T=${T%.fasta}_1kbtrim.fasta
seqret -sequence "$f"  -outseq "$T" -sbegin 1000 -send -1000
done
