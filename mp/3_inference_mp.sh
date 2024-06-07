#!/bin/bash
#SBATCH -J py_gyc
#SBATCH --partition gpu
#SBATCH --gpus t4:1
#SBATCH -N 1
#SBATCH -n 3
#SBATCH -o stderr_out/stdout.%j
#SBATCH -e stderr_out/stderr.%j         

module load conda
conda activate uni_core

data_path="/home/gaoyuchen/mpp/1_mpbp_20230905/7_knowledge_new/1_mp/10_knowledge/12/20"  # replace to your data path
results_path="3_infer_mp"  # replace to your results path
weight_path="2_save_finetune_mp/checkpoint_best.pt"  # replace to your ckpt path
batch_size=32
task_name='mp' # data folder name 
task_num=1
loss_func='finetune_smooth_mae'
dict_name='dict.txt'
conf_size=11
only_polar=0

CUDA_VISIBLE_DEVICES="0" 
python ./unimol/infer.py --user-dir ./unimol $data_path --task-name $task_name --valid-subset test \
       --results-path $results_path \
       --num-workers 1 --ddp-backend=c10d --batch-size $batch_size \
       --task mol_finetune --loss $loss_func --arch unimol_base \
       --classification-head-name $task_name --num-classes $task_num \
       --dict-name $dict_name --conf-size $conf_size \
       --only-polar $only_polar \
       --path $weight_path \
       --fp16 --fp16-init-scale 4 --fp16-scale-window 256 \
       --log-interval 50 --log-format simple \
