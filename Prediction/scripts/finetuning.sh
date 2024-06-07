data_path="your_path"  # replace to your data path
save_dir="save_finetune_xxx"  # replace to your save path
n_gpu=1
MASTER_PORT=10086
dict_name="dict.txt"
weight_path="ckp/mol_pre_no_h_220816.pt"  # replace to your ckpt path
task_name="property"  # molecular property prediction task name 
task_num=1
loss_func="finetune_smooth_mae"
lr=1e-4
batch_size=32
epoch=500
dropout=0
warmup=0.06
local_batch_size=32
only_polar=0
conf_size=11
seed=0

if [ "$task_name" == "qm7dft" ] || [ "$task_name" == "qm8dft" ] || [ "$task_name" == "qm9dft" ] || [ "$task_name" == "ep" ] || [ "$task_name" == "lipo" ] || [ "$task_name" == "mp" ] || [ "$task_name" == "bp" ] || [ "$task_name" == "fp" ]; then
	metric="valid_agg_mae"
elif [ "$task_name" == "esol" ] || [ "$task_name" == "freesolv" ]; then
    metric="valid_agg_rmse"
else 
    metric="valid_agg_auc"
fi

export NCCL_ASYNC_ERROR_HANDLING=1
export OMP_NUM_THREADS=1
update_freq=`expr $batch_size / $local_batch_size`
python -m torch.distributed.launch --nproc_per_node=$n_gpu --master_port=$MASTER_PORT $(which unicore-train) $data_path --task-name $task_name --user-dir ./unimol --train-subset train --valid-subset valid \
       --conf-size $conf_size \
       --num-workers 8 --ddp-backend=c10d \
       --dict-name $dict_name \
       --task mol_finetune --loss $loss_func --arch unimol_base  \
       --classification-head-name $task_name --num-classes $task_num \
       --optimizer adam --adam-betas "(0.9, 0.99)" --adam-eps 1e-6 --clip-norm 1.0 \
       --lr-scheduler polynomial_decay --lr $lr --max-epoch $epoch --batch-size $local_batch_size --pooler-dropout $dropout\
       --update-freq $update_freq --seed $seed \
       --fp16 --fp16-init-scale 4 --fp16-scale-window 256 \
       --log-interval 100 --log-format simple \
       --validate-interval 1 \
       --finetune-from-model $weight_path \
       --best-checkpoint-metric $metric --patience 20 \
       --save-dir $save_dir --only-polar $only_polar \
       --reg

# --reg, for regression task
# --maximize-best-checkpoint-metric, for classification task