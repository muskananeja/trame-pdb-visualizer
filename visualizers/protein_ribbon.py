from trame.widgets import vtk as vtk_widgets
import vtk

# Use vtk directly instead of vtkmodules
vtkDataObject = vtk.vtkDataObject
vtkPDBReader = vtk.vtkPDBReader
vtkActor = vtk.vtkActor
vtkRenderer = vtk.vtkRenderer
vtkRenderWindow = vtk.vtkRenderWindow
vtkRenderWindowInteractor = vtk.vtkRenderWindowInteractor
vtkPolyDataMapper = vtk.vtkPolyDataMapper
vtkCellPicker = vtk.vtkCellPicker
vtkCommand = vtk.vtkCommand
vtkTransform = vtk.vtkTransform
vtkSphereSource = vtk.vtkSphereSource
vtkCylinderSource = vtk.vtkCylinderSource
vtkGlyphSource2D = vtk.vtkGlyphSource2D

# Protein ribbon specific
vtkProteinRibbonFilter = vtk.vtkProteinRibbonFilter
vtkMoleculeMapper = vtk.vtkMoleculeMapper

import numpy as np

from .base import BaseVisualizer
from utils import extract_data_arrays

class ProteinRibbonVisualizer(BaseVisualizer):
    """Protein Ribbon visualization style for PDB files"""
    
    def create_visualization(self, pdb_file, color_mapper=None, state=None, ctrl=None):
        """Create a Protein Ribbon visualization of the PDB file"""
        # Create renderer and window
        renderer = vtkRenderer()
        renderWindow = vtkRenderWindow()
        renderWindow.AddRenderer(renderer)
        
        # PDB reader
        reader = vtkPDBReader()
        reader.SetFileName(pdb_file)
        reader.Update()
        
        # Get data arrays
        data_arrays = extract_data_arrays(reader)
        
        # Create protein ribbon filter
        ribbon = vtkProteinRibbonFilter()
        ribbon.SetInputConnection(reader.GetOutputPort())
        ribbon.Update()
        
        # Apply color mapping if provided
        if color_mapper:
            actor = color_mapper.apply_to_protein_ribbon(ribbon, renderer)
        else:
            # Default visualization
            mapper = vtkPolyDataMapper()
            mapper.SetInputConnection(ribbon.GetOutputPort())
            
            actor = vtkActor()
            actor.SetMapper(mapper)
            
            # Add actor to renderer
            renderer.AddActor(actor)
        
        # Setup picker for interaction
        picker = vtkCellPicker()
        picker.SetTolerance(0.005)
        
        # Setup callback for hover
        def handle_mouse_move(obj, event):
            global_point = obj.GetEventPosition()
            result = picker.Pick(global_point[0], global_point[1], 0, renderer)
            
            if result != 0:
                picked_position = picker.GetPickPosition()
                cell_id = picker.GetCellId()
                
                if cell_id >= 0:
                    # For ribbon visualization, we can't directly map to atoms
                    # So we just display the position
                    info = f"Position: ({picked_position[0]:.2f}, {picked_position[1]:.2f}, {picked_position[2]:.2f})"
                    state.hover_info = info
                else:
                    state.hover_info = ""
            else:
                state.hover_info = ""
        
        # Set up interaction (if we have an interactor)
        interactor = renderWindow.GetInteractor()
        if interactor is not None:
            interactor.AddObserver(vtkCommand.MouseMoveEvent, handle_mouse_move)
        
        # Reset camera
        renderer.ResetCamera()
        
        # Create simple HTML view reference for our wrapper
        view_html = f'<div id="vtk-view-protein-ribbon"></div>'
        
        # Update state with the view
        state.view_protein_ribbon = view_html
        
        return view_html
