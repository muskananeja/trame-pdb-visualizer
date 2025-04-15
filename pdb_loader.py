import vtk
from Bio import PDB

def load_pdb_to_vtk(file_path):
    """
    Loads a PDB file and converts it into a VTK Molecule object for rendering.
    """
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("molecule", file_path)

    molecule = vtk.vtkMolecule()

    atom_index_map = {}  # Maps atom ID to VTK atom index

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_vtk_id = molecule.AppendAtom(atom.element, *atom.coord)
                    atom_index_map[atom.serial_number] = atom_vtk_id

    # Add bonds (simple distance-based check)
    for model in structure:
        for chain in model:
            for residue in chain:
                atom_list = list(residue)
                for i, atom1 in enumerate(atom_list):
                    for j, atom2 in enumerate(atom_list[i + 1 :]):
                        if atom1 - atom2 < 1.6:  # Bond threshold in Ångstroms
                            molecule.AppendBond(
                                atom_index_map[atom1.serial_number], atom_index_map[atom2.serial_number], 1
                            )

    return molecule
