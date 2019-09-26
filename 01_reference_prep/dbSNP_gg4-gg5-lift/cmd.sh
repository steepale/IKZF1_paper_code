grep -v "#" Gallus_gallus.vcf |grep "TSA=SNV" |grep -v "," |awk '{print $1"\t"$2-1"\t"$2}' >bedfile4Control-freec.bed 
