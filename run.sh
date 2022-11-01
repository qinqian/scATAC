#!/bin/bash -ex


#bs download run -i 219537323 -o .

echo "Lane,Sample,Index" > sample_sheet1.csv
echo "1-4,MSK93202,SI-NA-A1" >> sample_sheet1.csv
echo "1-4,MSK74711,SI-NA-A2" >> sample_sheet1.csv

__conda_setup="$('/data/pinello/SHARED_SOFTWARE/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/data/pinello/SHARED_SOFTWARE/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/data/pinello/SHARED_SOFTWARE/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/data/pinello/SHARED_SOFTWARE/anaconda3/bin:$PATH"
    fi  
fi
unset __conda_setup
conda activate sc-tutorial       ## load with vim, emacs

#./cellranger-atac-2.0.0/cellranger-atac mkfastq --run data/ --csv sample_sheet1.csv

#./cellranger-atac-2.0.0/cellranger-atac count --id=MSK74711 \
#                      --reference=refdata-cellranger-atac-GRCh38-and-mm10-2020-A-2.0.0/ \
#                      --fastqs=HLJJKBGXJ/outs/fastq_path/HLJJKBGXJ/MSK74711/ \
#                      --sample=MSK74711 \
#                      --localcores=8 \
#                      --localmem=64

#./cellranger-atac-2.0.0/cellranger-atac count --id=MSK93202 \
#                      --reference=refdata-cellranger-atac-GRCh38-and-mm10-2020-A-2.0.0/ \
#                      --fastqs=HLJJKBGXJ/outs/fastq_path/HLJJKBGXJ/MSK93202/ \
#                      --sample=MSK93202 \
#                      --localcores=20 \
#                      --localmem=64
