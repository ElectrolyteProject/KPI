import os
import pickle
import pandas as pd
import numpy as np
from tqdm import tqdm
import time

def get_csv_results(predict_path, csv_path):
    predict = pd.read_pickle(predict_path)
    # print(predict)
    # ep_id = []
    smiles = []
    predict_list0 = []
    target_list0 = []

    for batch in predict:
        sz = batch["bsz"]
        for i in range(sz):
            # print(i)
            # ep_id.append(batch["smi_name"][i])
            smiles.append(batch["smi_name"][i])
            predict_list0.append(batch["predict"][i][0].cpu().tolist())
            target_list0.append(batch["target"][i][0].cpu().tolist())
            
    predict_df = pd.DataFrame({"SMILES": smiles,\
                                "MP_pred": predict_list0,
                                "MP_targ": target_list0,})


    predict_df = predict_df.groupby("SMILES")[["MP_pred", "MP_targ"]].mean().reset_index()
    predict_df.to_csv(csv_path,index=False)
    # print(predict_df)
    return predict_df


print(time.strftime('%Y.%m.%d %H:%M:%S', time.localtime(time.time())))

# test
start_time = time.time()
predict_path='./3_infer_mp/2_save_finetune_mp_test.out.pkl'  # replace to your results path
csv_path='./3_infer_mp/3_infer_mp_20230906.csv'

# # total
# start_time = time.time()
# predict_path='./infer_total/save_finetune_data.out.pkl'  # replace to your results path
# csv_path='./infer_total/total_difficult_single_results_epoch_20_20221211.csv'

predict_df = get_csv_results(predict_path, csv_path)
predict_df.info(), predict_df.head()

print("--- Count: %s seconds ---" % (time.time() - start_time))