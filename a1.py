from trame.app import get_server
from trame.ui.vuetify import SinglePageLayout
from trame.widgets import vtk, vuetify

import os
import tempfile

from vtkmodules.vtkIOChemistry import vtkPDBReader
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkPolyDataMapper,
    vtkRenderer,
    vtkRenderWindow,
    vtkRenderWindowInteractor
)
from vtkmodules.vtkCommonDataModel import vtkMolecule
from vtkmodules.vtkDomainsChemistry import vtkMoleculeMapper
from vtkmodules.vtkDomainsChemistry import vtkProteinRibbonFilter
from vtkmodules.vtkInteractionStyle import vtkInteractorStyleTrackballCamera

# ---------------------------------------------------------------------
# VTK Pipeline for PDB Rendering
# ---------------------------------------------------------------------
CURRENT_DIRECTORY = os.path.abspath(os.path.dirname(__file__))

# Create a renderer and render window
renderer = vtkRenderer()
renderWindow = vtkRenderWindow()
renderWindow.AddRenderer(renderer)

renderWindowInteractor = vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)

# Set up the camera interactor style
interactor_style = vtkInteractorStyleTrackballCamera()
renderWindowInteractor.SetInteractorStyle(interactor_style)

# Load an initial PDB File (could be a default placeholder)
pdb_reader = vtkPDBReader()
default_pdb = os.path.join(CURRENT_DIRECTORY, "../data/102m.pdb")
pdb_reader.SetFileName(default_pdb)
pdb_reader.Update()


# Get Output Data
molecule = pdb_reader.GetOutput()
num_atoms = molecule.GetNumberOfPoints()
print(f"Number of Atoms: {num_atoms}")

# Print available data arrays in PDB file
point_data = molecule.GetPointData()
print("Available Data Arrays:")
for i in range(point_data.GetNumberOfArrays()):
    print(f"- {point_data.GetArrayName(i)}")

# Convert PDB to Molecule Data
molecule = vtkMolecule()
molecule.DeepCopy(pdb_reader.GetOutput())

# Set up the molecule mapper and actor (atoms + bonds)
molecule_mapper = vtkMoleculeMapper()
molecule_mapper.SetInputData(molecule)

molecule_actor = vtkActor()
molecule_actor.SetMapper(molecule_mapper)

# Set up the ribbon representation for proteins
ribbon_filter = vtkProteinRibbonFilter()
ribbon_filter.SetInputConnection(pdb_reader.GetOutputPort())

ribbon_mapper = vtkPolyDataMapper()
ribbon_mapper.SetInputConnection(ribbon_filter.GetOutputPort())

ribbon_actor = vtkActor()
ribbon_actor.SetMapper(ribbon_mapper)

# Add both actors to the renderer
renderer.AddActor(molecule_actor)  # Atoms & Bonds Representation
renderer.AddActor(ribbon_actor)    # Ribbon Representation

renderer.ResetCamera()

# ---------------------------------------------------------------------
# Trame Server and File Upload Callback
# ---------------------------------------------------------------------
server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

@state.change("pdb_file")
def update_pdb_file(pdb_file, **kwargs):
    if not pdb_file:
        return

    # Get the uploaded file bytes (join if a list)
    pdb_bytes = pdb_file.get("content")
    if isinstance(pdb_bytes, list):
        pdb_bytes = b"".join(pdb_bytes)

    # Write the uploaded bytes to a temporary file (so vtkPDBReader can read it)
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as tmp:
        tmp.write(pdb_bytes)
        tmp.flush()
        pdb_filename = tmp.name

    # Update the VTK pipeline with the new file
    pdb_reader.SetFileName(pdb_filename)
    pdb_reader.Update()

    # Update the molecule data and mapper
    molecule.DeepCopy(pdb_reader.GetOutput())
    molecule_mapper.SetInputData(molecule)
    molecule_mapper.Update()

    # Update the ribbon filter and mapper
    ribbon_filter.SetInputConnection(pdb_reader.GetOutputPort())
    ribbon_filter.Update()
    ribbon_mapper.Update()

    # Reset the camera and update the view
    renderer.ResetCamera()
    ctrl.view_update()

# ---------------------------------------------------------------------
# Trame UI
# ---------------------------------------------------------------------
with SinglePageLayout(server) as layout:
    layout.title.set_text("PDB Viewer")

    with layout.toolbar:
        # File input for PDB upload
        vuetify.VFileInput(
            v_model=("pdb_file", None),
            label="Upload PDB File",
            accept=".pdb",
            dense=True,
            hide_details=True,
        )

    with layout.content:
        with vuetify.VContainer(
            fluid=True,
            classes="pa-0 fill-height",
        ):
            # VTK local view using the render window
            view = vtk.VtkLocalView(renderWindow)
            ctrl.view_update = view.update
            ctrl.view_reset_camera = view.reset_camera

# ---------------------------------------------------------------------
# Main Execution
# ---------------------------------------------------------------------
if __name__ == "__main__":
    server.start()