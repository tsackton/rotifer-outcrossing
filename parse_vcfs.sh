vcftools --vcf CR.vcf --remove-indels -c --max-alleles 2 --extract-FORMAT-info AD > CR.ad.out
vcftools --vcf CR.vcf -c --get-INFO AC > CR.ac.out
vcftools --vcf MA.vcf --remove-indels -c --max-alleles 2 --extract-FORMAT-info AD > MA.ad.out
vcftools --vcf MA.vcf -c --get-INFO AC > MA.ac.out
vcftools --vcf MM.vcf --remove-indels -c --max-alleles 2 --extract-FORMAT-info AD > MM.ad.out
vcftools --vcf MM.vcf -c --get-INFO AC > MM.ac.out


##notes on triallelic sites

wc -l *.ac.out
374390 CR.ac.out
318466 MA.ac.out
461617 MM.ac.out

grep -c "," *.ac.out
CR.ac.out:25222
MA.ac.out:21815
MM.ac.out:30081

wc -l *.ad.out
298153 CR.ad.out
253915 MA.ad.out
369429 MM.ad.out

egrep -c "\s+[ACTG],[ACTG]\s+" *.ac.out
CR.ac.out:3593
MA.ac.out:3054
MM.ac.out:4463

triallelic:
CR: 1.19%
MA: 1.19%
MM: 1.19%
