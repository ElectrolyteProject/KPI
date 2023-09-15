import json
import os
import argparse
import pandas as pd
from rdkit.Chem import Draw
from rdkit import Chem
from rdkit.Chem.Draw.MolDrawing import DrawingOptions  # Only needed if modifying defaults


def print_graph():
    filename = args.f
    with open(args.f) as a:
        dict = json.load(a)
    suppl = [Chem.MolFromSmiles(value) for key, value in dict.items()]
    molnamelist = [(key, value) for key, value in dict.items()]
    print("文件{}中共包含{}个分子".format(filename, len(suppl)))

    img = Draw.MolsToGridImage(
        suppl,
        molsPerRow=5,
        subImgSize=(300, 300),
        #legends=[f'{x} {y}' for (x, y) in molnamelist]
        legends=[f'{x}' for (x, y) in molnamelist],
    )
    img.save(f'{filename.split(".")[0]}.jpg')


def print_graph_single():
    with open(args.f) as a:
        dict = json.load(a)
    index = 0
    for key, value in dict.items():
        index+=1
        mol = Chem.MolFromSmiles(value)
        img = Draw.MolsToGridImage(
            [mol],
            molsPerRow=1,
            subImgSize=(300, 300),
            #legends=[f'{key}']
        )
        img.save(f'{index}-{key}.jpg')


def print_graph_new():
    with open(args.f) as a:
        dict = json.load(a)
    opts = DrawingOptions()
    opts.includeAtomNumbers = True
    opts.bondLineWidth = 2.8
    opts.atomLabelFontSize = 25
    #opts.atomLabelFontFace = 'Arial'
    #IPythonConsole.drawOptions.fontFile = r'C:\Windows\Fonts\Arial.ttf'
    for key, value in dict.items():
        m = Chem.MolFromSmiles(value)
        draw = Draw.MolToImage(m, options=opts)
        draw.save(f'{key}.png')
        #break


def main():
    #print_graph()
    print_graph_single()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', help='File path of the .json which includes epid-smiles.')
    args = parser.parse_args()
    main()

