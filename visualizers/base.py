from abc import ABC, abstractmethod

class BaseVisualizer(ABC):
    """Base class for all visualizers"""
    
    @abstractmethod
    def create_visualization(self, pdb_file, color_mapper=None, state=None, ctrl=None):
        """
        Create a visualization of the PDB file
        
        Parameters
        ----------
        pdb_file : str
            Path to the PDB file
        color_mapper : BaseColorMapper
            Color mapper to use for the visualization
        state : trame.state
            Trame state object
        ctrl : trame.controller
            Trame controller object
            
        Returns
        -------
        str
            HTML representation of the view
        """
        pass
