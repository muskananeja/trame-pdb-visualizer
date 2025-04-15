import vtk

# Use vtk directly instead of vtkmodules
vtkDataObject = vtk.vtkDataObject

def extract_data_arrays(reader):
    """
    Extract data arrays from a VTK reader (vtkPDBReader)
    
    Parameters
    ----------
    reader : vtkPDBReader
        The PDB reader to extract data from
        
    Returns
    -------
    list
        A list of dictionaries with metadata about the arrays
    """
    result = []
    
    # Function to process a dataset
    def process_dataset(dataset, association):
        if not dataset:
            return
            
        # Point data
        if association == vtkDataObject.POINT:
            data = dataset.GetPointData()
        # Cell data    
        else:
            data = dataset.GetCellData()
            
        # Iterate through arrays
        for i in range(data.GetNumberOfArrays()):
            array = data.GetArray(i)
            if not array:
                continue
                
            name = array.GetName()
            
            # Get range if array has tuples
            if array.GetNumberOfTuples() > 0:
                scalar_range = array.GetRange()
            else:
                scalar_range = None
                
            result.append({
                "name": name,
                "index": i,
                "range": scalar_range,
                "association": "point" if association == vtkDataObject.POINT else "cell"
            })
    
    # Process point and cell data for output 0 (atoms)
    output0 = reader.GetOutput(0)
    if output0:
        process_dataset(output0, vtkDataObject.POINT)
        process_dataset(output0, vtkDataObject.CELL)
    
    # Process point and cell data for output 1 (bonds)
    output1 = reader.GetOutput(1)
    if output1:
        process_dataset(output1, vtkDataObject.POINT)
        process_dataset(output1, vtkDataObject.CELL)
    
    return result
