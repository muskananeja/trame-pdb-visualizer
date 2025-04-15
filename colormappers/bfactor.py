import numpy as np
import vtk

# Use vtk directly instead of vtkmodules
vtkLookupTable = vtk.vtkLookupTable
vtkActor = vtk.vtkActor
vtkPolyDataMapper = vtk.vtkPolyDataMapper

from .base import BaseColorMapper

class BFactorColorMapper(BaseColorMapper):
    """Color mapper that colors by B-factor (temperature factor)"""
    
    def apply_to_ball_and_stick(self, reader, renderer):
        """Apply B-factor coloring to Ball and Stick visualization"""
        # Get the output and determine B-factor range
        atoms = reader.GetOutput(0)
        b_factors = atoms.GetPointData().GetArray("b_factor")
        
        if b_factors:
            b_factor_range = b_factors.GetRange()
        else:
            # Default range if B-factors not available
            b_factor_range = (0, 100)
        
        # Create lookup table for B-factor coloring
        lut = vtkLookupTable()
        lut.SetHueRange(0.0, 0.667)  # Blue to red
        lut.SetTableRange(b_factor_range)
        lut.Build()
        
        # Atoms mapper with B-factor coloring
        atoms_mapper = vtkPolyDataMapper()
        atoms_mapper.SetInputConnection(reader.GetOutputPort(0))
        atoms_mapper.SetScalarModeToUsePointFieldData()
        atoms_mapper.SelectColorArray("b_factor")
        atoms_mapper.SetLookupTable(lut)
        atoms_mapper.SetScalarRange(b_factor_range)
        
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
        """Apply B-factor coloring to Protein Ribbon visualization"""
        # Create lookup table for B-factor coloring
        lut = vtkLookupTable()
        lut.SetHueRange(0.0, 0.667)  # Blue to red
        lut.SetTableRange(0, 100)  # Typical B-factor range
        lut.Build()
        
        # Create mapper with B-factor coloring
        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(ribbon.GetOutputPort())
        mapper.SetScalarModeToUsePointFieldData()
        mapper.SelectColorArray("b_factor")
        mapper.SetLookupTable(lut)
        mapper.SetScalarRange(0, 100)
        
        # Create actor
        actor = vtkActor()
        actor.SetMapper(mapper)
        
        # Add actor to renderer
        renderer.AddActor(actor)
        
        return [actor]
