import deepchem as dc
import kpi_plots as kp
import kpi_features as kf
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw.MolDrawing import DrawingOptions

def plot_cluster(input_path, folder_name1, folder_name2, _featurizer, _color, _property, _perplexity=35):
    """
    Generates a cluster plot using the specified parameters and saves the output and image.

    Parameters:
    - input_path: str, path to the input file.
    - folder_name1: str, name of the folder to save the output CSV file.
    - folder_name2: str, name of the folder to save the output image file.
    - _featurizer: str, the featurizer method to use.
    - _color: str, the color scheme to use in the plot.
    - _property: str, the property to visualize.
    - _perplexity: int, the perplexity parameter for the clustering algorithm (default is 35).
    """
    output_path = f"./{folder_name1}/cluster_featurizer_{_featurizer}_perplexity_{_perplexity}.csv"
    image_path = f"./{folder_name2}/cluster_featurizer_{_featurizer}_perplexity_{_perplexity}.jpg"
    kp.plot_cluster(input_path, output_path, image_path, _featurizer, _perplexity, _color, _property)
    print(f"Cluster plot saved to {image_path} and data saved to {output_path}")

def find_furthest_molecules(csv_file_path):
    """
    Finds the maximum Euclidean distance between any two molecules based on their t-SNE coordinates.

    Parameters:
    - csv_file_path: str, path to the CSV file containing t-SNE coordinates.

    Returns:
    - float, the maximum distance between any two molecules.
    """
    # Read the CSV file
    df = pd.read_csv(csv_file_path)
    
    # Extract t-SNE coordinates
    coordinates = df[['TSNE_x', 'TSNE_y']].values
    
    # Initialize maximum distance
    max_distance = 0
    
    # Calculate the distance between each pair of molecules
    for i in range(len(coordinates)):
        for j in range(i + 1, len(coordinates)):
            # Calculate Euclidean distance
            distance = np.linalg.norm(coordinates[i] - coordinates[j])
            # Update maximum distance
            if distance > max_distance:
                max_distance = distance

    # Return the maximum distance
    return max_distance


def find_similar_molecules(original_csv, new_csv, task_name, target_smiles, rou1, rou2, props):
    """
    Finds molecules similar to a target molecule based on t-SNE coordinates and a specified property.

    Parameters:
    - original_csv: str, path to the original CSV file containing molecule data.
    - new_csv: str, path to the output CSV file to save similar molecules.
    - task_name: str, the property name to compare.
    - target_smiles: str, the SMILES string of the target molecule.
    - rou1: float, the minimum distance threshold.
    - rou2: float, the maximum distance threshold.
    - props: float, the maximum allowed difference in the specified property.

    Returns:
    - None
    """
    # Read the CSV file
    df = pd.read_csv(original_csv)
    
    # Find the target molecule's row
    target_row = df[df['SMILES'] == target_smiles]
    
    if target_row.empty:
        print("No matching molecule found.")
        return
    
    # Print target molecule data
    print("Target Molecule Data:")
    print(target_row)
    
    # Extract target molecule's coordinates and property
    target_tsne_x = target_row['TSNE_x'].values[0]
    target_tsne_y = target_row['TSNE_y'].values[0]
    target_bp = target_row[task_name].values[0]
    
    # Define a function to check similarity
    def is_similar(row):
        distance = np.sqrt((row['TSNE_x'] - target_tsne_x)**2 + (row['TSNE_y'] - target_tsne_y)**2)
        prop_diff = abs(row[task_name] - target_bp)
        return rou1 < distance < rou2 and prop_diff <= props
    
    # Apply the function to filter similar molecules
    similar_molecules = df.apply(is_similar, axis=1)
    result_df = df[similar_molecules]
    
    # Save results to a new CSV file
    result_df.to_csv(new_csv, index=False)
    print(f"Results saved to {new_csv}")


