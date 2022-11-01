#!/bin/bash -ex

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

./cellranger-atac-2.0.0/cellranger-atac count --id=MSK74711_hg38 \
                      --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
                      --fastqs=HLJJKBGXJ/outs/fastq_path/HLJJKBGXJ/MSK74711/ \
                      --sample=MSK74711 \
                      --localcores=20 \
                      --localmem=64
