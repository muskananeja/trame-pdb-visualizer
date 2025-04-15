#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import vtk
import numpy as np

def extract_molecule_data(pdb_reader):
    """
    Extract metadata from the VTK PDB reader.
    
    Args:
        pdb_reader: vtkPDBReader instance with a loaded PDB file
        
    Returns:
        dict: Dictionary containing metadata about the molecule
    """
    # Get the output data from the reader
    output = pdb_reader.GetOutput()
    
    # Basic molecule statistics
    num_atoms = output.GetNumberOfPoints()
    num_bonds = output.GetNumberOfCells()
    
    # Extract array metadata (similar to the molecule.py reference)
    arrays_info = []
    
    # Process point data (atom data)
    pd = output.GetPointData()
    for i in range(pd.GetNumberOfArrays()):
        array = pd.GetArray(i)
        array_name = array.GetName()
        
        # Extract range info for scalar data
        if array.GetNumberOfComponents() == 1:
            scalar_range = array.GetRange()
        else:
            scalar_range = None
        
        arrays_info.append({
            "name": array_name,
            "association": "points",
            "size": array.GetNumberOfTuples(),
            "range": scalar_range
        })
    
    # Process cell data (bond data)
    cd = output.GetCellData()
    for i in range(cd.GetNumberOfArrays()):
        array = cd.GetArray(i)
        array_name = array.GetName()
        
        # Extract range info for scalar data
        if array.GetNumberOfComponents() == 1:
            scalar_range = array.GetRange()
        else:
            scalar_range = None
        
        arrays_info.append({
            "name": array_name,
            "association": "cells",
            "size": array.GetNumberOfTuples(),
            "range": scalar_range
        })
    
    # Compile all information
    molecule_info = {
        "num_atoms": num_atoms,
        "num_bonds": num_bonds,
        "arrays": arrays_info
    }
    
    # Extract global bounds
    bounds = output.GetBounds()
    molecule_info["bounds"] = bounds
    
    # Calculate molecule dimensions
    x_size = bounds[1] - bounds[0]
    y_size = bounds[3] - bounds[2]
    z_size = bounds[5] - bounds[4]
    molecule_info["dimensions"] = [x_size, y_size, z_size]
    
    # Get unique residue types
    unique_residues = set()
    residue_array = pd.GetArray("residue")
    if residue_array:
        for i in range(residue_array.GetNumberOfTuples()):
            unique_residues.add(int(residue_array.GetValue(i)))
    molecule_info["unique_residues"] = sorted(list(unique_residues))
    
    # Get unique atom types
    unique_atoms = set()
    atom_types_array = pd.GetArray("atom_types")
    if atom_types_array:
        for i in range(atom_types_array.GetNumberOfTuples()):
            unique_atoms.add(atom_types_array.GetValue(i))
    molecule_info["unique_atoms"] = sorted(list(unique_atoms))
    
    return molecule_info

def get_statistics_text(molecule_info):
    """Generate a formatted text summary of the molecule statistics"""
    if not molecule_info:
        return "No molecule loaded"
    
    text = []
    text.append(f"Atoms: {molecule_info.get('num_atoms', 'N/A')}")
    text.append(f"Bonds: {molecule_info.get('num_bonds', 'N/A')}")
    
    if 'dimensions' in molecule_info:
        dim = molecule_info['dimensions']
        text.append(f"Size: {dim[0]:.2f} × {dim[1]:.2f} × {dim[2]:.2f} Å")
    
    if 'unique_residues' in molecule_info:
        text.append(f"Residues: {len(molecule_info['unique_residues'])}")
    
    if 'unique_atoms' in molecule_info:
        text.append(f"Atom types: {len(molecule_info['unique_atoms'])}")
    
    return "\n".join(text)
