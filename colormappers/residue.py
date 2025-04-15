import numpy as np
import vtk

# Use vtk directly instead of vtkmodules
vtkLookupTable = vtk.vtkLookupTable
vtkActor = vtk.vtkActor
vtkPolyDataMapper = vtk.vtkPolyDataMapper

from .base import BaseColorMapper

class ResidueColorMapper(BaseColorMapper):
    """Color mapper that colors by residue type"""
    
    def apply_to_ball_and_stick(self, reader, renderer):
        """Apply residue-based coloring to Ball and Stick visualization"""
        # Create lookup table for residue coloring
        lut = vtkLookupTable()
        lut.SetNumberOfTableValues(20)  # Common amino acids
        
        # Define colors for common amino acids
        # Each color is RGBA (Red, Green, Blue, Alpha)
        lut.SetTableValue(0, 1.0, 0.0, 0.0, 1.0)  # ALA - Red
        lut.SetTableValue(1, 0.0, 1.0, 0.0, 1.0)  # ARG - Green
        lut.SetTableValue(2, 0.0, 0.0, 1.0, 1.0)  # ASN - Blue
        lut.SetTableValue(3, 1.0, 1.0, 0.0, 1.0)  # ASP - Yellow
        lut.SetTableValue(4, 1.0, 0.0, 1.0, 1.0)  # CYS - Magenta
        lut.SetTableValue(5, 0.0, 1.0, 1.0, 1.0)  # GLN - Cyan
        lut.SetTableValue(6, 0.5, 0.0, 0.0, 1.0)  # GLU - Dark Red
        lut.SetTableValue(7, 0.0, 0.5, 0.0, 1.0)  # GLY - Dark Green
        lut.SetTableValue(8, 0.0, 0.0, 0.5, 1.0)  # HIS - Dark Blue
        lut.SetTableValue(9, 0.5, 0.5, 0.0, 1.0)  # ILE - Olive
        lut.SetTableValue(10, 0.5, 0.0, 0.5, 1.0)  # LEU - Purple
        lut.SetTableValue(11, 0.0, 0.5, 0.5, 1.0)  # LYS - Teal
        lut.SetTableValue(12, 0.7, 0.5, 0.5, 1.0)  # MET - Salmon
        lut.SetTableValue(13, 0.5, 0.7, 0.5, 1.0)  # PHE - Light Green
        lut.SetTableValue(14, 0.5, 0.5, 0.7, 1.0)  # PRO - Light Blue
        lut.SetTableValue(15, 0.8, 0.7, 0.6, 1.0)  # SER - Tan
        lut.SetTableValue(16, 0.6, 0.8, 0.7, 1.0)  # THR - Mint
        lut.SetTableValue(17, 0.7, 0.6, 0.8, 1.0)  # TRP - Lavender
        lut.SetTableValue(18, 0.8, 0.8, 0.8, 1.0)  # TYR - Light Grey
        lut.SetTableValue(19, 0.3, 0.3, 0.3, 1.0)  # VAL - Dark Grey
        
        lut.Build()
        
        # Atoms mapper with residue coloring
        atoms_mapper = vtkPolyDataMapper()
        atoms_mapper.SetInputConnection(reader.GetOutputPort(0))
        atoms_mapper.SetScalarModeToUsePointFieldData()
        atoms_mapper.SelectColorArray("residue")
        atoms_mapper.SetLookupTable(lut)
        atoms_mapper.SetScalarRange(0, 19)
        
        # Bonds mapper
        bonds_mapper = vtkPolyDataMapper()
        bonds_mapper.SetInputConnection(reader.GetOutputPort(1))
        
        # Create actors
        atoms_actor = vtkActor()
        atoms_actor.SetMapper(atoms_mapper)
        
        bonds_actor = vtkActor()
        bonds_actor.SetMapper(bonds_mapper)
        bonds_actor.GetProperty().SetColor(0.8, 0.8, 0.8)  # Grey bonds
        
        # Add actors to renderer
        renderer.AddActor(atoms_actor)
        renderer.AddActor(bonds_actor)
        
        return [atoms_actor, bonds_actor]
    
    def apply_to_protein_ribbon(self, ribbon, renderer):
        """Apply residue-based coloring to Protein Ribbon visualization"""
        # Create lookup table for residue coloring
        lut = vtkLookupTable()
        lut.SetNumberOfTableValues(20)  # Common amino acids
        
        # Define colors for common amino acids (same as above)
        lut.SetTableValue(0, 1.0, 0.0, 0.0, 1.0)  # ALA - Red
        lut.SetTableValue(1, 0.0, 1.0, 0.0, 1.0)  # ARG - Green
        lut.SetTableValue(2, 0.0, 0.0, 1.0, 1.0)  # ASN - Blue
        lut.SetTableValue(3, 1.0, 1.0, 0.0, 1.0)  # ASP - Yellow
        lut.SetTableValue(4, 1.0, 0.0, 1.0, 1.0)  # CYS - Magenta
        lut.SetTableValue(5, 0.0, 1.0, 1.0, 1.0)  # GLN - Cyan
        lut.SetTableValue(6, 0.5, 0.0, 0.0, 1.0)  # GLU - Dark Red
        lut.SetTableValue(7, 0.0, 0.5, 0.0, 1.0)  # GLY - Dark Green
        lut.SetTableValue(8, 0.0, 0.0, 0.5, 1.0)  # HIS - Dark Blue
        lut.SetTableValue(9, 0.5, 0.5, 0.0, 1.0)  # ILE - Olive
        lut.SetTableValue(10, 0.5, 0.0, 0.5, 1.0)  # LEU - Purple
        lut.SetTableValue(11, 0.0, 0.5, 0.5, 1.0)  # LYS - Teal
        lut.SetTableValue(12, 0.7, 0.5, 0.5, 1.0)  # MET - Salmon
        lut.SetTableValue(13, 0.5, 0.7, 0.5, 1.0)  # PHE - Light Green
        lut.SetTableValue(14, 0.5, 0.5, 0.7, 1.0)  # PRO - Light Blue
        lut.SetTableValue(15, 0.8, 0.7, 0.6, 1.0)  # SER - Tan
        lut.SetTableValue(16, 0.6, 0.8, 0.7, 1.0)  # THR - Mint
        lut.SetTableValue(17, 0.7, 0.6, 0.8, 1.0)  # TRP - Lavender
        lut.SetTableValue(18, 0.8, 0.8, 0.8, 1.0)  # TYR - Light Grey
        lut.SetTableValue(19, 0.3, 0.3, 0.3, 1.0)  # VAL - Dark Grey
        
        lut.Build()
        
        # Create mapper with residue coloring
        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(ribbon.GetOutputPort())
        mapper.SetLookupTable(lut)
        mapper.SetScalarRange(0, 19)
        
        # Create actor
        actor = vtkActor()
        actor.SetMapper(mapper)
        
        # Add actor to renderer
        renderer.AddActor(actor)
        
        return [actor]
