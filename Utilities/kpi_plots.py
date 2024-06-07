import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Descriptors
import deepchem as dc
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from matplotlib.patches import Circle

def plot_histogram(data, color='blue', label='Molwt', 
                   x_major=0, y_major=0, x_minor=0, y_minor=0,
                   xlim_min=0, xlim_max=0, ylim_min=0, ylim_max=0):
    """
    Plot a histogram with the given data and settings.

    Parameters:
    data (list): Data to plot in the histogram.
    color (str): Color of the bars in the histogram.
    label (str): Label for the x-axis.
    x_major (int): Major ticks for the x-axis.
    y_major (int): Major ticks for the y-axis.
    x_minor (int): Minor ticks for the x-axis.
    y_minor (int): Minor ticks for the y-axis.
    xlim_min (int): Minimum limit for the x-axis.
    xlim_max (int): Maximum limit for the x-axis.
    ylim_min (int): Minimum limit for the y-axis.
    ylim_max (int): Maximum limit for the y-axis.
    """

    font_settings = {'family': 'arial', 'weight': 'normal', 'size': 16}

    plt.xlabel(f'{label}', font_settings)
    plt.ylabel('Count', font_settings)
    plt.tick_params(axis='both', which='major', length=6, width=2, direction='in', labelsize=14)
    plt.tick_params(axis='both', which='minor', length=3, width=1, direction='in', labelsize=14)

    ax = plt.gca()
    plt.hist(data, bins=50, color=color, edgecolor='black')

    for spine in ax.spines.values():
        spine.set_linewidth(2)

    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'

    if x_major:
        ax.xaxis.set_major_locator(plt.MultipleLocator(x_major))
    if y_major:
        ax.yaxis.set_major_locator(plt.MultipleLocator(y_major))
    if x_minor:
        ax.xaxis.set_minor_locator(plt.MultipleLocator(x_minor))
    if y_minor:
        ax.yaxis.set_minor_locator(plt.MultipleLocator(y_minor))

    if xlim_min or xlim_max:
        plt.xlim(xlim_min, xlim_max)
    if ylim_min or ylim_max:
        plt.ylim(ylim_min, ylim_max)

    plt.subplots_adjust(left=0.15)
    plt.show()
    
def analyze_data_distribution(df, property_column):
    """
    Analyze the distribution of molecular properties and plot histograms.

    Parameters:
    df (DataFrame): DataFrame containing SMILES strings and properties.
    property_column (str): Column name for the property to analyze.
    """

    elements = ['C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl', 'Br', 'I']
    element_count = {element: 0 for element in elements}
    molecular_weights = []
    atom_counts = []
    property_values = []

    # Iterate through each row in the DataFrame
    for _, row in df.iterrows():
        smiles = row['SMILES']
        prop = row[property_column]
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            # Count unique elements in the molecule
            unique_elements = set(atom.GetSymbol() for atom in mol.GetAtoms() if atom.GetSymbol() in element_count)
            for element in unique_elements:
                element_count[element] += 1

            # Compute molecular weight and atom count
            molecular_weights.append(Descriptors.MolWt(mol))
            atom_counts.append(mol.GetNumAtoms())
            # Record property value
            property_values.append(prop)

    # Print element distribution statistics
    print("Element distribution statistics:")
    for element, count in element_count.items():
        print(f"{element}: {count} molecules contain")

    # Plot element distribution
    plt.bar(element_count.keys(), element_count.values(), color='blue')
    plt.xlabel('Element')
    plt.ylabel('Number of Molecules Containing Element')
    plt.title('Distribution of Elements in Molecules')
    plt.show()

    # Plot property distribution
    plot_histogram(property_values, color='green', label=f'Molecular Property ({property_column})')

    # Calculate and print statistical metrics
    max_value = np.max(property_values)
    min_value = np.min(property_values)
    mean_value = np.mean(property_values)
    std_dev = np.std(property_values)

    print(f"Max value: {max_value}")
    print(f"Min value: {min_value}")
    print(f"Mean value: {mean_value}")
    print(f"Standard deviation: {std_dev:.2f}")

    # Plot distribution of the number of atoms in molecules
    plot_histogram(atom_counts, color='pink', label='Number of Atoms in Molecule')

    # Plot distribution of molecular weights
    plot_histogram(molecular_weights, color='purple', label='Molecular Weight')

def plot_cluster(input_path, output_path, image_path, _featurizer, _perplexity, _color, _property):
    """
    Generates a cluster plot using t-SNE based on the specified parameters and saves the output and image.

    Parameters:
    - input_path: str, path to the input CSV file.
    - output_path: str, path to the output CSV file.
    - image_path: str, path to the output image file.
    - _featurizer: str, the featurizer method to use ('ECFP' or 'MACCS').
    - _perplexity: int, the perplexity parameter for the t-SNE algorithm.
    - _color: str, the color scheme to use in the plot.
    - _property: str, the property to visualize.
    """
    all_smiles = []
    all_props = []

    # Read the input CSV file
    with open(input_path, 'r') as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            all_smiles.append(row[0])
            all_props.append(eval(row[1]))

    print(f"Total SMILES: {len(all_smiles)}")

    # Choose the featurizer
    if _featurizer == "ECFP":
        featurizer = dc.feat.CircularFingerprint()
    else:
        featurizer = dc.feat.MACCSKeysFingerprint()

    # Featurize the SMILES
    features = featurizer.featurize(all_smiles)
    dataset = dc.data.NumpyDataset(features)
    
    # Apply t-SNE for dimensionality reduction
    X_tsne = TSNE(n_components=2, init='pca', perplexity=_perplexity, random_state=0).fit_transform(dataset.X)

    # Write the t-SNE results to the output CSV file
    with open(output_path, 'w') as f:
        writer = csv.writer(f)
        writer.writerow(["SMILES", "TSNE_x", "TSNE_y", _property])
        for smiles, x, y, prop in zip(all_smiles, X_tsne[:, 0], X_tsne[:, 1], all_props):
            writer.writerow([smiles, x, y, prop])

    # Plot the t-SNE results
    plt.figure(figsize=(6, 5))
    plt.scatter(X_tsne[:, 0], X_tsne[:, 1], s=5, c=all_props, cmap=_color, label="t-SNE")

    # Customize plot appearance
    ax = plt.gca()
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)
    ax.spines['top'].set_linewidth(2)

    plt.xticks([])
    plt.yticks([])
    plt.axis('off')

    # Save and show the plot
    plt.savefig(image_path, dpi=120)
    plt.show()

def plot_molecules_with_circles(csv_file_path, smiles_list, rou2, rou4):
    """
    Plots t-SNE coordinates of molecules from a CSV file, marks specified molecules with different markers,
    and draws circles around them with specified radii.

    Parameters:
    - csv_file_path: str, path to the CSV file containing t-SNE coordinates and molecular data.
    - smiles_list: list of str, list of SMILES strings to mark and draw circles around.
    - rou2: float, radius of the first circle to draw around specified molecules.
    - rou4: float, radius of the second circle to draw around specified molecules.

    Returns:
    - None
    """
    # Read the CSV file
    df = pd.read_csv(csv_file_path)
    
    # Plot scatter plot colored by melting point (mp)
    scatter = plt.scatter(df['TSNE_x'], df['TSNE_y'], c=df['mp'], cmap='viridis', alpha=0.5)
    plt.colorbar(scatter, label='property')
    plt.xlabel('TSNE_x')
    plt.ylabel('TSNE_y')
    plt.title('Molecular Visualization with t-SNE Coordinates')

    # Define markers for specified SMILES
    markers = ['*', '^']  # Assume up to two molecules need marking
    
    # Mark specified SMILES and draw circles
    for index, smiles in enumerate(smiles_list):
        # Find the corresponding molecule
        target_row = df[df['SMILES'] == smiles]
        if not target_row.empty:
            x, y = target_row['TSNE_x'].values[0], target_row['TSNE_y'].values[0]
            # Use different markers for specified molecules
            plt.scatter(x, y, color='red', s=50, marker=markers[index % len(markers)], edgecolors='black', label=f'SMILES: {smiles}')
            # Draw the first circle
            circle1 = Circle((x, y), rou2, color='blue', fill=False, linewidth=1.5, linestyle='-')
            plt.gca().add_artist(circle1)
            # Draw the second circle
            circle2 = Circle((x, y), rou4, color='red', fill=False, linewidth=1.5, linestyle='-')
            plt.gca().add_artist(circle2)

    plt.legend()
    plt.grid(True)
    plt.show()
