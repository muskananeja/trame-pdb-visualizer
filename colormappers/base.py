from abc import ABC, abstractmethod

class BaseColorMapper(ABC):
    """Base class for all color mappers"""
    
    @abstractmethod
    def apply_to_ball_and_stick(self, reader, renderer):
        """
        Apply color mapping to Ball and Stick visualization
        
        Parameters
        ----------
        reader : vtkPDBReader
            The PDB reader
        renderer : vtkRenderer
            The renderer to add actors to
            
        Returns
        -------
        list
            List of actors added to the renderer
        """
        pass
    
    @abstractmethod
    def apply_to_protein_ribbon(self, ribbon, renderer):
        """
        Apply color mapping to Protein Ribbon visualization
        
        Parameters
        ----------
        ribbon : vtkProteinRibbonFilter
            The protein ribbon filter
        renderer : vtkRenderer
            The renderer to add actors to
            
        Returns
        -------
        list
            List of actors added to the renderer
        """
        pass
