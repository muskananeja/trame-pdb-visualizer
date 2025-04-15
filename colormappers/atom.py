import numpy as np
import vtk

# Use vtk directly instead of vtkmodules
vtkColorSeries = vtk.vtkColorSeries
vtkLookupTable = vtk.vtkLookupTable
vtkActor = vtk.vtkActor
vtkPolyDataMapper = vtk.vtkPolyDataMapper
vtkPeriodicTable = vtk.vtkPeriodicTable

from .base import BaseColorMapper

class AtomColorMapper(BaseColorMapper):
    """Color mapper that colors by atom type"""
    
    def __init__(self):
        # Create a periodic table for atom colors
        self.periodic_table = vtkPeriodicTable()
    
    def apply_to_ball_and_stick(self, reader, renderer):
        """Apply atom-based coloring to Ball and Stick visualization"""
        # Create mappers and actors
        atoms_mapper = vtkPolyDataMapper()
        atoms_mapper.SetInputConnection(reader.GetOutputPort(0))
        
        bonds_mapper = vtkPolyDataMapper()
        bonds_mapper.SetInputConnection(reader.GetOutputPort(1))
        
        # Create atom actor with per-atom coloring
        atoms_actor = vtkActor()
        atoms_actor.SetMapper(atoms_mapper)
        atoms_actor.GetProperty().SetRepresentationToSurface()
        atoms_actor.GetProperty().SetInterpolationToGouraud()
        
        # Create bond actor
        bonds_actor = vtkActor()
        bonds_actor.SetMapper(bonds_mapper)
        bonds_actor.GetProperty().SetColor(0.8, 0.8, 0.8)  # Grey bonds
        
        # Add actors to renderer
        renderer.AddActor(atoms_actor)
        renderer.AddActor(bonds_actor)
        
        return [atoms_actor, bonds_actor]
    
    def apply_to_protein_ribbon(self, ribbon, renderer):
        """Apply atom-based coloring to Protein Ribbon visualization"""
        # Create lookup table for coloring
        lut = vtkLookupTable()
        color_series = vtkColorSeries()
        color_series.SetColorSchemeByName("Brewer Qualitative Paired")
        color_series.BuildLookupTable(lut)
        
        # Create mapper with the lookup table
        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(ribbon.GetOutputPort())
        mapper.SetLookupTable(lut)
        mapper.SetScalarRange(0, 19)  # Typical range for amino acids
        
        # Create actor
        actor = vtkActor()
        actor.SetMapper(mapper)
        
        # Add actor to renderer
        renderer.AddActor(actor)
        
        return [actor]
