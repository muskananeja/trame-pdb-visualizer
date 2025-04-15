#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import vtk
import numpy as np

# Color mapping constants
ATOM_COLORS = {
    'H': [1.0, 1.0, 1.0],  # White
    'C': [0.5, 0.5, 0.5],  # Grey
    'N': [0.0, 0.0, 1.0],  # Blue
    'O': [1.0, 0.0, 0.0],  # Red
    'S': [1.0, 1.0, 0.0],  # Yellow
    'P': [1.0, 0.5, 0.0],  # Orange
    'Cl': [0.0, 1.0, 0.0], # Green
    'Ca': [0.5, 0.5, 0.5], # Grey
    'Fe': [0.7, 0.5, 0.0], # Brown
    'Na': [0.0, 0.0, 1.0], # Blue
    'K': [0.8, 0.6, 1.0],  # Purple
    'Zn': [0.5, 0.5, 0.5], # Grey
    'Mg': [0.0, 1.0, 0.0], # Green
    # Default for other atoms
    'default': [0.9, 0.9, 0.9] # Light Grey
}

RESIDUE_COLORS = {
    # Hydrophobic
    'ALA': [0.8, 0.8, 0.8],  # Light Grey
    'VAL': [0.8, 0.8, 0.8],  # Light Grey
    'LEU': [0.8, 0.8, 0.8],  # Light Grey
    'ILE': [0.8, 0.8, 0.8],  # Light Grey
    'PRO': [0.8, 0.8, 0.8],  # Light Grey
    'PHE': [0.8, 0.8, 0.8],  # Light Grey
    'MET': [0.8, 0.8, 0.8],  # Light Grey
    'TRP': [0.8, 0.8, 0.8],  # Light Grey
    
    # Polar
    'GLY': [0.0, 1.0, 0.0],  # Green
    'SER': [0.0, 1.0, 0.0],  # Green
    'THR': [0.0, 1.0, 0.0],  # Green
    'CYS': [0.0, 1.0, 0.0],  # Green
    'TYR': [0.0, 1.0, 0.0],  # Green
    'ASN': [0.0, 1.0, 0.0],  # Green
    'GLN': [0.0, 1.0, 0.0],  # Green
    
    # Acidic
    'ASP': [1.0, 0.0, 0.0],  # Red
    'GLU': [1.0, 0.0, 0.0],  # Red
    
    # Basic
    'LYS': [0.0, 0.0, 1.0],  # Blue
    'ARG': [0.0, 0.0, 1.0],  # Blue
    'HIS': [0.0, 0.0, 1.0],  # Blue
    
    # Default
    'default': [0.9, 0.9, 0.9]  # Light Grey
}

def create_ball_stick_visualization(pdb_reader):
    """
    Create a Ball and Stick visualization for the molecule.
    
    Args:
        pdb_reader: vtkPDBReader instance with a loaded PDB file
        
    Returns:
        vtkActor: Actor for rendering
    """
    # Create the ball and stick representation
    ball_stick = vtk.vtkMoleculeMapper()
    ball_stick.SetInputConnection(pdb_reader.GetOutputPort())
    ball_stick.UseBallAndStickSettings()
    
    # Set rendering parameters
    ball_stick.SetAtomicRadiusScaleFactor(0.3)  # Smaller atoms for clearer visualization
    ball_stick.SetBondRadius(0.08)
    
    # Create an actor for the ball and stick representation
    actor = vtk.vtkActor()
    actor.SetMapper(ball_stick)
    
    return actor

def create_ribbon_visualization(pdb_reader):
    """
    Create a Protein Ribbon visualization for the molecule.
    
    Args:
        pdb_reader: vtkPDBReader instance with a loaded PDB file
        
    Returns:
        vtkActor: Actor for rendering
    """
    # Create the ribbon representation using ProteinRibbonFilter
    ribbon_filter = vtk.vtkProteinRibbonFilter()
    ribbon_filter.SetInputConnection(pdb_reader.GetOutputPort())
    ribbon_filter.SetDrawSmallMolecules(True)
    
    # Mapper for the ribbon
    ribbon_mapper = vtk.vtkPolyDataMapper()
    ribbon_mapper.SetInputConnection(ribbon_filter.GetOutputPort())
    
    # Create an actor for the ribbon
    actor = vtk.vtkActor()
    actor.SetMapper(ribbon_mapper)
    
    return actor

def apply_color_mapping(actor, pdb_reader, color_mapping):
    """
    Apply a specific color mapping to the molecule visualization.
    
    Args:
        actor: vtkActor to apply color mapping to
        pdb_reader: vtkPDBReader instance with a loaded PDB file
        color_mapping: String indicating the color mapping to use ("atom", "bfactor", or "residue")
    """
    if color_mapping == "atom":
        # In Ball and Stick mode, the atoms are colored by element by default
        # This is handled by the vtkMoleculeMapper
        pass
    
    elif color_mapping == "bfactor":
        # Create a lookup table for B-factor coloring
        lut = vtk.vtkLookupTable()
        lut.SetHueRange(0.667, 0.0)  # Blue to red
        lut.SetTableRange(0.0, 100.0)  # Typical B-factor range
        lut.Build()
        
        # Get the mapper from the actor
        mapper = actor.GetMapper()
        
        # If a simple mapper is used (vtkPolyDataMapper)
        if isinstance(mapper, vtk.vtkPolyDataMapper):
            mapper.ScalarVisibilityOn()
            mapper.SetScalarModeToUsePointFieldData()
            mapper.SelectColorArray("bfactor")
            mapper.SetLookupTable(lut)
        # If a molecule mapper is used (vtkMoleculeMapper)
        else:
            # We need to convert to a standard mapper to use B-factor coloring
            output = pdb_reader.GetOutput()
            polydata = vtk.vtkPolyData()
            polydata.ShallowCopy(output)
            
            new_mapper = vtk.vtkPolyDataMapper()
            new_mapper.SetInputData(polydata)
            new_mapper.ScalarVisibilityOn()
            new_mapper.SetScalarModeToUsePointFieldData()
            new_mapper.SelectColorArray("bfactor")
            new_mapper.SetLookupTable(lut)
            
            actor.SetMapper(new_mapper)
    
    elif color_mapping == "residue":
        # Create a lookup table for residue coloring
        lut = vtk.vtkLookupTable()
        lut.SetNumberOfTableValues(20)  # 20 common amino acids
        lut.Build()
        
        # Get the mapper from the actor
        mapper = actor.GetMapper()
        
        # If a simple mapper is used (vtkPolyDataMapper)
        if isinstance(mapper, vtk.vtkPolyDataMapper):
            mapper.ScalarVisibilityOn()
            mapper.SetScalarModeToUsePointFieldData()
            mapper.SelectColorArray("residue")
            mapper.SetLookupTable(lut)
        # If a molecule mapper is used (vtkMoleculeMapper)
        else:
            # We need to convert to a standard mapper to use residue coloring
            output = pdb_reader.GetOutput()
            polydata = vtk.vtkPolyData()
            polydata.ShallowCopy(output)
            
            new_mapper = vtk.vtkPolyDataMapper()
            new_mapper.SetInputData(polydata)
            new_mapper.ScalarVisibilityOn()
            new_mapper.SetScalarModeToUsePointFieldData()
            new_mapper.SelectColorArray("residue")
            new_mapper.SetLookupTable(lut)
            
            actor.SetMapper(new_mapper)
