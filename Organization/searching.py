import csv
import pubchempy as pcp
import requests
import json
from chemspipy import ChemSpider

def get_compound_data(smiles, cs):
    """
    Searches for compound data on ChemSpider using a SMILES string and returns the melting point, boiling point, and flash point.

    Parameters:
    - smiles: str, the SMILES string of the compound.
    - cs: ChemSpider, an instance of the ChemSpider class initialized with an API key.

    Returns:
    - tuple: (melting_point, boiling_point, flash_point) if data is found, otherwise (None, None, None).
    """
    results = cs.search(smiles)
    if results:
        compound = results[0]  # Assuming the first result is the desired one
        return compound.melting_point, compound.boiling_point, compound.flash_point
    else:
        return None, None, None

def ChemSpider_get_data(input_file, output_file, api_key):
    """
    Reads a CSV file containing SMILES strings, retrieves compound data from ChemSpider, and saves the results to a new CSV file.

    Parameters:
    - input_file: str, path to the input CSV file containing SMILES strings.
    - output_file: str, path to the output CSV file to save the retrieved compound data.
    - api_key: str, ChemSpider API key for accessing the ChemSpider database.
    """
    cs = ChemSpider(api_key)

    with open(input_file, newline='') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile)
        fieldnames = ['SMILES', 'Melting Point', 'Boiling Point', 'Flash Point']
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()

        for row in reader:
            smiles = row['SMILES']
            melting_point, boiling_point, flash_point = get_compound_data(smiles, cs)
            writer.writerow({
                'SMILES': smiles,
                'Melting Point': melting_point,
                'Boiling Point': boiling_point,
                'Flash Point': flash_point
            })

    print("Data processing complete.")

def read_smiles_from_csv(file_path):
    """
    Reads SMILES strings from a CSV file.

    Parameters:
    - file_path: str, path to the CSV file.

    Returns:
    - smiles: list of str, list of SMILES strings.
    """
    smiles = []
    with open(file_path, mode='r', encoding='utf-8') as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader)  # Skip header row if there is one
        for row in csvreader:
            smiles.append(row[0])
    return smiles

def get_properties(cid):
    """
    Retrieves chemical and physical properties for a given compound CID from PubChem.

    Parameters:
    - cid: int, PubChem Compound ID.

    Returns:
    - properties: dict, dictionary containing melting point, boiling point, and flash point.
    """
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON/?heading=Chemical+and+Physical+Properties'
    req = requests.get(url)
    proper_json = json.loads(req.text)
    
    properties = {'Melting Point': '', 'Boiling Point': '', 'Flash Point': ''}
    try:
        sections = proper_json['Record']['Section'][0]['Section'][1]['Section']
        for section in sections:
            if section['TOCHeading'] in properties:
                properties[section['TOCHeading']] = section['Information'][0]['Value']['StringWithMarkup'][0]['String']
    except KeyError:
        pass  # No data found for this compound

    return properties

def get_cids_and_properties_from_smiles(smi, results):
    """
    Retrieves PubChem CID and properties for a given SMILES string and appends the results to a list.

    Parameters:
    - smi: str, SMILES string.
    - results: list, list to append the results.
    """
    try:
        compounds = pcp.get_compounds(smi, 'smiles')
        if compounds:
            cid = compounds[0].cid
            properties = get_properties(cid)
            results.append([smi, cid, properties['Melting Point'], properties['Boiling Point'], properties['Flash Point']])
        else:
            results.append([smi, None, '', '', ''])
    except Exception as e:
        print(f"Error processing {smi}: {e}")
        results.append([smi, None, '', '', ''])

def write_to_csv(data, filename):
    """
    Writes data to a CSV file.

    Parameters:
    - data: list of lists, data to write to the CSV file.
    - filename: str, path to the output CSV file.
    """
    with open(filename, 'w', newline='', encoding='utf-8') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['SMILES', 'CID', 'Melting Point', 'Boiling Point', 'Flash Point'])
        for row in data:
            csvwriter.writerow(row)

def PubChem_get_data(input_csv, output_csv):
    """
    Main function to read SMILES strings, retrieve properties, and write the results to a CSV file.

    Parameters:
    - input_csv: str, path to the input CSV file containing SMILES strings.
    - output_csv: str, path to the output CSV file.
    """
    smiles_list = read_smiles_from_csv(input_csv)
    
    results = []
    for i, smi in enumerate(smiles_list, 1):
        get_cids_and_properties_from_smiles(smi, results)
        print(f"{i}: {smi}")
    
    write_to_csv(results, output_csv)
