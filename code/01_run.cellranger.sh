#! /bin/bash

samplelist=(YF1 YF2 YF3 YM1 YM2 YM3 OF1 OF2 OF3 OF4 OM1 OM2 OM3)

for sample in ${samplelist[*]}
do
/software/CellRanger/cellranger-6.0.2/cellranger count \
--id=$sample \
--transcriptome=/software/Reference_genome/Macaca_fascicularis_5.0.102/Macaca_fascicularis \
--fastqs=/object/08_monkey_cochlea/snRNA_seq/Rawdata/$sample \
--sample=$sample \
--include-introns \
--localcores=30 \
--chemistry=SC3Pv3
done