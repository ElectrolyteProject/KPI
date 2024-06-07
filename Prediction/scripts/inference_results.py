import pandas as pd
import time

def get_csv_results(predict_path, csv_path):
    """
    Reads prediction data from a pickle file, processes it, and saves the results to a CSV file.

    Parameters:
    - predict_path: str, path to the pickle file containing prediction data.
    - csv_path: str, path to the output CSV file.

    Returns:
    - predict_df: DataFrame, a DataFrame containing the processed results.
    """
    # Read the pickle file
    predict = pd.read_pickle(predict_path)
    smiles = []
    predict_list0 = []
    target_list0 = []

    # Process each batch in the prediction data
    for batch in predict:
        sz = batch["bsz"]
        for i in range(sz):
            smiles.append(batch["smi_name"][i])
            predict_list0.append(batch["predict"][i][0].cpu().tolist())
            target_list0.append(batch["target"][i][0].cpu().tolist())

    # Create a DataFrame from the processed data
    predict_df = pd.DataFrame({
        "SMILES": smiles,
        "pred": predict_list0,
        "targ": target_list0
    })

    # Group by SMILES and calculate the mean of predictions and targets
    predict_df = predict_df.groupby("SMILES")[["pred", "targ"]].mean().reset_index()
    
    # Save the DataFrame to a CSV file
    predict_df.to_csv(csv_path, index=False)
    
    return predict_df

if __name__ == "__main__":
    # Print the current time
    print(time.strftime('%Y.%m.%d %H:%M:%S', time.localtime(time.time())))

    # Start timing
    start_time = time.time()
    
    # Define file paths
    predict_path = './infer_xxx/save_finetune_xxx_test.out.pkl'  # Replace with your results path
    csv_path = './infer_xxx/infer_fp_xxx.csv'  # Replace with your desired output path

    # Process the predictions and save to CSV
    predict_df = get_csv_results(predict_path, csv_path)
    
    # Display DataFrame info and head
    print(predict_df.info())
    print(predict_df.head())

    # Print the elapsed time
    print("--- Count: %s seconds ---" % (time.time() - start_time))
