from trame.widgets import vuetify, html

def create_viewer(state, ctrl):
    """
    Create the main viewer component
    
    Parameters
    ----------
    state : trame.state
        Trame state object
    ctrl : trame.controller
        Trame controller object
    """
    with vuetify.VContainer(fluid=True, classes="fill-height pa-0"):
        with html.Div(style="position: relative; width: 100%; height: 100%;"):
            # Simple welcome message (always shown for now)
            with html.Div(
                style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%); text-align: center;"
            ):
                html.H3("Welcome to Molecule Visualizer")
                html.P("Upload a PDB file or load the example to get started")
                vuetify.VIcon("mdi-molecule", size="100", color="primary")
