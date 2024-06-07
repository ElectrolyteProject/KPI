import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors

def process_smiles_and_save(input_csv, output_csv, property_name='mp'):
    """
    Reads a CSV file containing SMILES strings and a specified property (e.g., melting point),
    calculates molecular weight and heavy atom count, and saves the processed data to a new CSV file.

    Parameters:
    - input_csv: str, path to the input CSV file containing SMILES strings and the specified property.
    - output_csv: str, path to the output CSV file to save the processed data.
    - property_name: str, name of the property to include in the output (default is 'mp' for melting point).
    """
    # Read the CSV file
    df = pd.read_csv(input_csv)
    
    # Prepare output data dictionary
    output_data = {'SMILES': [], property_name: [], 'MolWt': [], '#Heavy': []}
    
    # Iterate through the SMILES column
    for index, row in df.iterrows():
        smiles = row['SMILES']
        property_value = row[property_name]
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # Calculate molecular weight and heavy atom count
            molecular_weight = Descriptors.MolWt(mol)
            atom_count = mol.GetNumHeavyAtoms()  # Heavy atom count, excluding hydrogen
            
            # Update output data
            output_data['SMILES'].append(smiles)
            output_data[property_name].append(property_value)
            output_data['MolWt'].append(molecular_weight)
            output_data['#Heavy'].append(atom_count)
    
    # Save the processed data to a new CSV file
    output_df = pd.DataFrame(output_data)
    output_df.to_csv(output_csv, index=False)
    
    # Print the number of SMILES processed
    print(f"Total SMILES processed: {len(output_data['SMILES'])}")

def plot_histograms(input_csv, variables, num_bins=50, facecolor='blue', alpha=0.5):
    """
    Reads a CSV file and plots histograms for specified variables.

    Parameters:
    - input_csv: str, path to the input CSV file.
    - variables: list of str, list of variables (column names) to plot histograms for.
    - num_bins: int, number of bins for the histograms (default is 50).
    - facecolor: str, color of the bars in the histograms (default is 'blue').
    - alpha: float, transparency level of the bars (default is 0.5).
    """
    # Read the CSV file
    df = pd.read_csv(input_csv)
    
    # Plot histograms for each specified variable
    for var in variables:
        plt.figure()  # Create a new figure
        plt.hist(df[var], num_bins, facecolor=facecolor, alpha=alpha)
        plt.title(f'Histogram of {var}')
        plt.xlabel(var)
        plt.ylabel('Density')
        plt.show()

def calculate_statistics(input_csv, columns_to_analyze):
    """
    Reads a CSV file and calculates statistics (mean, variance, max, min) for specified columns.

    Parameters:
    - input_csv: str, path to the input CSV file.
    - columns_to_analyze: list of str, list of columns to calculate statistics for.
    
    Returns:
    - statistics: dict, a dictionary containing the statistics for each column.
    """
    # Load the CSV file
    df = pd.read_csv(input_csv)
    
    # Initialize a dictionary to store the statistics
    statistics = {}
    
    # Calculate statistics for each specified column
    for column in columns_to_analyze:
        statistics[column] = {
            'mean': df[column].mean(),
            'variance': df[column].var(),
            'max': df[column].max(),
            'min': df[column].min()
        }
    
    # Print the results
    for column, stats in statistics.items():
        print(f"Statistics for {column}:")
        for stat, value in stats.items():
            print(f"  {stat}: {value}")
        print()  # Add a blank line for better readability
    
    return statistics

def compute_histogram_data(input_csv, output_csv, variables, num_bins=50):
    """
    Reads a CSV file, computes histogram data for specified variables, normalizes the index, and saves the results to a new CSV file.

    Parameters:
    - input_csv: str, path to the input CSV file.
    - output_csv: str, path to the output CSV file to save the histogram data.
    - variables: list of str, list of variables (column names) to compute histograms for.
    - num_bins: int, number of bins for the histograms (default is 50).
    """
    # Read the CSV file
    df = pd.read_csv(input_csv)

    # Initialize a DataFrame to save the histogram data
    hist_data = pd.DataFrame()

    # Compute histogram data for each specified variable
    for var in variables:
        # Use numpy.histogram to compute the histogram
        n, bins = np.histogram(df[var], bins=num_bins)
        
        # Add the histogram data to the DataFrame
        hist_data[f'{var}_n'] = n

    # Add an index column
    hist_data['Index'] = range(1, len(hist_data) + 1)

    # Normalize the index column
    hist_data['Normalized_Index'] = (hist_data['Index'] - hist_data['Index'].min()) / (hist_data['Index'].max() - hist_data['Index'].min())

    # Reorder the columns to have the index columns at the front
    cols = hist_data.columns.tolist()
    cols = cols[-2:] + cols[:-2]  # Move the last two columns to the front
    hist_data = hist_data[cols]

    # Save the histogram data to a CSV file
    hist_data.to_csv(output_csv, index=False)

    print(f"Histogram data saved to '{output_csv}'.")
    