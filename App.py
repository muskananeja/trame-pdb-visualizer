import os
import vtk
from trame.app import get_app
from trame.widgets import html, vtk
from trame.ui.vtk import VtkRemoteViewer
from Bio import PDB
from trame.widgets import upload

# Step 1: Load the PDB file and parse it with Biopython
def parse_pdb(pdb_file):
    # Initialize the PDB parser
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('Protein', pdb_file)
    
    # Extract atoms and bonds
    atoms = []
    bonds = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atoms.append(atom)
                    # Create bonds based on atom connectivity
                    if atom.get_bonded_atoms():
                        for bonded_atom in atom.get_bonded_atoms():
                            bonds.append((atom, bonded_atom))
    return atoms, bonds

# Step 2: Create VTK representations for atoms and bonds
def create_vtk_representation(atoms, bonds):
    # Create a VTK renderer
    renderer = vtk.vtkRenderer()
    atoms_polydata = vtk.vtkPolyData()

    # Create VTK points and a list of atoms
    points = vtk.vtkPoints()
    for atom in atoms:
        points.InsertNextPoint(atom.coord)
    atoms_polydata.SetPoints(points)

    # Create VTK spheres for atoms
    atom_spheres = vtk.vtkGlyph3D()
    atom_spheres.SetSourceConnection(vtk.vtkSphereSource().GetOutputPort())
    atom_spheres.SetInputData(atoms_polydata)
    atom_spheres.Update()

    # Create bonds as lines between atoms
    bond_lines = vtk.vtkCellArray()
    for bond in bonds:
        atom1, atom2 = bond
        id1 = points.InsertNextPoint(atom1.coord)
        id2 = points.InsertNextPoint(atom2.coord)
        bond_lines.InsertNextCell(2)
        bond_lines.InsertCellPoint(id1)
        bond_lines.InsertCellPoint(id2)

    # Create VTK polydata for bonds
    bond_polydata = vtk.vtkPolyData()
    bond_polydata.SetPoints(points)
    bond_polydata.SetLines(bond_lines)

    # Visualize bonds as lines and atoms as spheres
    bond_mapper = vtk.vtkPolyDataMapper()
    bond_mapper.SetInputData(bond_polydata)
    bond_actor = vtk.vtkActor()
    bond_actor.SetMapper(bond_mapper)

    atom_mapper = vtk.vtkPolyDataMapper()
    atom_mapper.SetInputData(atom_spheres.GetOutput())
    atom_actor = vtk.vtkActor()
    atom_actor.SetMapper(atom_mapper)

    # Add actors to the renderer
    renderer.AddActor(bond_actor)
    renderer.AddActor(atom_actor)
    
    return renderer

# Step 3: Setup the Trame App with file upload
def main():
    # Set up the Trame app
    app = get_app()

    # Upload handler
    def handle_file_upload(file_data, filename):
        if filename.endswith('.pdb'):
            # Step 4: Parse the uploaded PDB file
            atoms, bonds = parse_pdb(file_data)

            # Create the VTK renderer
            renderer = create_vtk_representation(atoms, bonds)

            # Set up VTK viewer using Trame
            viewer = VtkRemoteViewer()
            viewer.renderer = renderer
            viewer.camera.view_up = [0, 0, 1]

            # Update the viewer
            viewer.update()

            # Set the renderer and viewer
            viewer.renderer = renderer
            with app.layout:
                html.H1("3D Molecular Visualization")
                viewer

        else:
            print("Please upload a valid PDB file.")

    # File upload widget
    upload_widget = upload.FileUpload(
        on_file=handle_file_upload,
        label="Upload PDB File"
    )

    # Layout with file upload widget
    with app.layout:
        html.H1("Upload a PDB File for Visualization")
        upload_widget

    # Run the app
    app.run()

if __name__ == '__main__':
    main()
