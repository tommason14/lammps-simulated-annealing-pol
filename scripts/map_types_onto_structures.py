from glob import glob
import os
import re

"""
Takes a directory of xyz files extracted from an MD trajectory using TRAVIS, and 
uses a labelled set of reference molecules to add atom types to the unlabelled xyzs

Manually add atom types for each distinct molecule in the cluster and change the molecules
dictionary to reflect this.
"""

molecules = {'C6H11N2': 'ref_c2mim.xyz', 'C2H3O2': 'ref_ac.xyz'}

def molecules_in_order(xyz):
    """
    Reads the second line of an xyz file created by TRAVIS,
    and returns only molecular formulae of each molecule in the file
    """
    with open(xyz) as f:
        next(f)
        info = next(f)
    splitline = re.split(r'=|,\s|\[|\]', info)
    return [i for i in splitline if re.match('^[0-9A-Z]+$', i) and not re.match('^[0-9]+$|^RM$', i)]

def types_from_xyz(xyz):
    """
    Return a list of atomic symbols (the first column of coords) of an xyz file
    """
    with open(xyz) as f:
        return [i.split()[0] for i in f.readlines()[2:]]

def coords_from_xyz(xyz):
    """
    Returns a nested list of coordinates in an xyz file
    """
    with open(xyz) as f:
        return [i.split()[1:] for i in f.readlines()[2:]]

newdir = os.path.join(os.getcwd(), 'annotated')
if not os.path.isdir(newdir):
    os.mkdir(newdir)

xyzs = [x for x in glob('*xyz') if 'ref' not in x] 
for xyz in xyzs:
    mols = molecules_in_order(xyz)
    newtypes = [types_from_xyz(molecules[mol]) for mol in mols]
    newtypes = [t for types in newtypes for t in types] # unnest
    coords = coords_from_xyz(xyz)
    with open(f"{newdir}/{xyz}", "w") as new:
        new.write(f"{len(newtypes)}\n\n")
        for symbol, coord in zip(newtypes, coords):
            coord = " ".join(coord)
            new.write(f"{symbol} {coord}\n")
