#!/usr/bin/env zsh
set -euo pipefail
# 需要注意的是vcf文件的染色体号为1,2,3,4,5,Mt,Pt

perl ~/biosoft/bin/table_annovar.pl 2_Region.vcf ~/Documents/database/TAIR10_annovar/ \
--vcfinput --buildver AT --remove --protocol refGene --operation g --nastring . --outfile 3_Annovar
