from trame.widgets import vtk as vtk_widgets  # ← Comment this if unused
import vtk
import numpy as np

from .base import BaseVisualizer
from utils import extract_data_arrays

class BallAndStickVisualizer(BaseVisualizer):
    """Ball and Stick visualization style for PDB files."""

    def create_visualization(self, pdb_file, color_mapper=None, state=None, ctrl=None):
        renderer = vtk.vtkRenderer()
        render_window = vtk.vtkRenderWindow()
        render_window.AddRenderer(renderer)

        reader = vtk.vtkPDBReader()
        reader.SetFileName(pdb_file)
        reader.Update()

        data_arrays = extract_data_arrays(reader)
        atoms = reader.GetOutput().GetPoints()
        num_atoms = atoms.GetNumberOfPoints()

        # --- Atom rendering ---
        atom_source = vtk.vtkSphereSource()
        atom_source.SetRadius(0.3)  # Atom size
        atom_source.SetPhiResolution(20)
        atom_source.SetThetaResolution(20)

        atom_glyph = vtk.vtkGlyph3D()
        atom_glyph.SetSourceConnection(atom_source.GetOutputPort())
        atom_glyph.SetInputConnection(reader.GetOutputPort())
        atom_glyph.SetColorModeToColorByScalar()
        atom_glyph.ScalingOff()
        atom_glyph.Update()

        atom_mapper = vtk.vtkPolyDataMapper()
        atom_mapper.SetInputConnection(atom_glyph.GetOutputPort())
        atom_actor = vtk.vtkActor()
        atom_actor.SetMapper(atom_mapper)
        renderer.AddActor(atom_actor)

        # --- Bond rendering ---
        bond_append = vtk.vtkAppendPolyData()

        # Commented out old bond fetch:
        # bonds = reader.GetBonds()
        # bond_count = bonds.GetNumberOfBonds()
        # for i in range(bond_count): ...

        # Manually compute bonds using distance threshold
        bond_threshold = 1.8  # Max distance for bond in angstroms

        for i in range(num_atoms):
            pos1 = np.array(atoms.GetPoint(i))
            for j in range(i + 1, num_atoms):
                pos2 = np.array(atoms.GetPoint(j))
                distance = np.linalg.norm(pos1 - pos2)

                if 0.4 < distance < bond_threshold:  # Avoid overlapping atoms
                    center = (pos1 + pos2) / 2
                    direction = pos2 - pos1
                    length = np.linalg.norm(direction)

                    cylinder = vtk.vtkCylinderSource()
                    cylinder.SetRadius(0.07)
                    cylinder.SetResolution(20)
                    cylinder.SetHeight(length)
                    cylinder.CappingOn()
                    cylinder.Update()

                    transform = vtk.vtkTransform()
                    transform.Translate(center)

                    z_axis = np.array([0, 0, 1])
                    axis = np.cross(z_axis, direction)
                    angle = np.degrees(np.arccos(np.dot(z_axis, direction) / length))

                    if np.linalg.norm(axis) > 1e-6:
                        transform.RotateWXYZ(angle, *axis)

                    transform_filter = vtk.vtkTransformPolyDataFilter()
                    transform_filter.SetTransform(transform)
                    transform_filter.SetInputConnection(cylinder.GetOutputPort())
                    transform_filter.Update()

                    bond_append.AddInputData(transform_filter.GetOutput())

        bond_append.Update()

        bond_mapper = vtk.vtkPolyDataMapper()
        bond_mapper.SetInputConnection(bond_append.GetOutputPort())

        bond_actor = vtk.vtkActor()
        bond_actor.SetMapper(bond_mapper)
        bond_actor.GetProperty().SetColor(0.8, 0.8, 0.8)
        renderer.AddActor(bond_actor)

        # --- Picking & hover ---
        picker = vtk.vtkCellPicker()
        picker.SetTolerance(0.005)

        def handle_mouse_move(obj, event):
            x, y = obj.GetEventPosition()
            if picker.Pick(x, y, 0, renderer):
                cell_id = picker.GetCellId()
                if cell_id >= 0:
                    atoms = reader.GetAtoms()
                    residues = reader.GetResidues()
                    chains = reader.GetChains()
                    elements = reader.GetAtomType()

                    if atoms and cell_id < atoms.GetNumberOfTuples():
                        atom_name = atoms.GetValue(cell_id)
                        residue = residues.GetValue(cell_id)
                        chain = chains.GetValue(cell_id)
                        element = elements.GetValue(cell_id) if elements else "Unknown"

                        state.hover_info = (
                            f"Atom: {atom_name}, Residue: {residue}, "
                            f"Chain: {chain}, Element: {element}"
                        )
                    else:
                        state.hover_info = f"Cell ID: {cell_id}"
                else:
                    state.hover_info = ""
            else:
                state.hover_info = ""

        interactor = render_window.GetInteractor()
        if interactor:
            interactor.AddObserver(vtk.vtkCommand.MouseMoveEvent, handle_mouse_move)

        renderer.ResetCamera()

        view_html = '<div id="vtk-view-ball-and-stick"></div>'
        state.view_ball_and_stick = view_html

        return view_html
