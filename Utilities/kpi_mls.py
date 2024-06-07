import os
import csv
import joblib
import matplotlib.pyplot as plt
import deepchem as dc
import kpi_features as kf
from math import sqrt
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

def features_standard(data, props):
    """
    Standardize features extracted from SMILES strings and return the dataset.
    
    Parameters:
    data (DataFrame): DataFrame containing SMILES strings and properties.
    props (str): The name of the property column to be used as the target variable.
    
    Returns:
    list: A list containing the standardized features and the dataset in DeepChem's NumpyDataset format.
    """
    smiles = data.loc[data["SMILES"].notnull(), "SMILES"].values
    props_list = data.loc[data[props].notnull(), props].values
    X = [kf.extract_features(smi) for smi in smiles]
    y = list(props_list)
    
    stdsc = StandardScaler()
    X_dataset = stdsc.fit_transform(X)
    standard_dataset = dc.data.NumpyDataset(X_dataset, y)
    
    return [X_dataset, standard_dataset]

def try_different_method(train_dataset, test_dataset, props, model, method):
    """
    Train the model, evaluate it on the test dataset, and save the results and metrics.
    
    Parameters:
    train_dataset (NumpyDataset): The training dataset.
    test_dataset (NumpyDataset): The testing dataset.
    props (str): The name of the property being predicted.
    model (sklearn estimator): The machine learning model to be trained.
    method (str): The name of the method used for training.
    """
    model.fit(train_dataset.X, train_dataset.y)
    score = model.score(test_dataset.X, test_dataset.y)
    
    result_train = model.predict(train_dataset.X)
    result_test = model.predict(test_dataset.X)
    
    folder_name = "2_results"
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        
    # Save the model
    joblib.dump(model, f'{folder_name}/{method}.model')
    
    # Save the predictions
    with open(f"{folder_name}/{method}.csv", 'w') as f:
        writer = csv.writer(f)
        for actual, predicted in zip(train_dataset.y, result_train):
            writer.writerow([actual, predicted, 1])
        for actual, predicted in zip(test_dataset.y, result_test):
            writer.writerow([actual, predicted, 2])
    
    print(f"Method: {method} --- Score: {score:.3f}")
    
    # Calculate and print performance metrics
    metrics = {
        "MAE": mean_absolute_error(test_dataset.y, result_test),
        "MSE": mean_squared_error(test_dataset.y, result_test),
        "RMSE": sqrt(mean_squared_error(test_dataset.y, result_test)),
        "R2": r2_score(test_dataset.y, result_test)
    }
    
    print(f"Train MAE: {mean_absolute_error(train_dataset.y, result_train):.3f}")
    print(f"Train MSE: {mean_squared_error(train_dataset.y, result_train):.3f}")
    print(f"Train RMSE: {sqrt(mean_squared_error(train_dataset.y, result_train)):.3f}")
    print(f"Train R2: {r2_score(train_dataset.y, result_train):.3f}")
    print(f"Test MAE: {metrics['MAE']:.3f}")
    print(f"Test MSE: {metrics['MSE']:.3f}")
    print(f"Test RMSE: {metrics['RMSE']:.3f}")
    print(f"Test R2: {metrics['R2']:.3f}")
    
    # Save metrics to a CSV file
    with open(f"1_{props}_metrics.csv", 'a') as f:
        writer = csv.writer(f)
        writer.writerow([method, metrics["MAE"], metrics["MSE"], metrics["RMSE"], metrics["R2"]])
    
    # Plot the results
    plt.figure(dpi=300)
    ax = plt.gca()
    
    plt.scatter(train_dataset.y, result_train, c='#5E9AEC', s=3, label="Train dataset")
    plt.scatter(test_dataset.y, result_test, c='#F47575', s=3, label="Test dataset")
    
    font = {'family': 'Arial', 'weight': 'normal', 'size': 12}
    plt.xlabel('Exp. boiling point (℃)', font)
    plt.ylabel('Pred. boiling point (℃)', font)
    
    for spine in ax.spines.values():
        spine.set_linewidth(2)
    
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    
    plt.legend(prop=font, frameon=False)
    x = [min(result_train) - 1, max(result_train) + 1]
    y = x
    plt.plot(x, y, ls='--', c='k', alpha=0.5)
    plt.show()
