from trame.widgets import vuetify, html

def create_drawer(drawer, visualization_styles, color_mappings, state, ctrl):
    """
    Create the drawer UI with controls
    
    Parameters
    ----------
    drawer : layout.drawer or None
        The drawer layout component (can be None for direct embedding)
    visualization_styles : list
        List of available visualization styles
    color_mappings : list
        List of available color mappings
    state : trame.state
        Trame state object
    ctrl : trame.controller
        Trame controller object
    """
    with vuetify.VContainer(fluid=True, classes="pa-0"):
        # Title
        vuetify.VSubheader("Molecule Visualization Controls")
        
        # Separator
        vuetify.VDivider()
        
        # File Upload
        with vuetify.VCard(flat=True, classes="pa-2"):
            with vuetify.VCardTitle(classes="text-subtitle-1 py-1"):
                html.Div("PDB File Upload")
            
            with vuetify.VCardText(classes="py-1"):
                with vuetify.VRow(classes="pa-0", dense=True):
                    with html.Div(style="width: 100%;"):
                        vuetify.VFileInput(
                            label="Upload PDB File",
                            accept=".pdb",
                            show_size=True,
                            truncate_length=15,
                            v_model=("file_content", None),
                            on_change=ctrl.on_file_upload,
                        )
            
            with vuetify.VCardActions():
                vuetify.VSpacer()
                vuetify.VBtn(
                    "Load Example",
                    small=True,
                    color="primary",
                    click=ctrl.load_example,
                )
        
        # Separator
        vuetify.VDivider()
        
        # Visualization Style Dropdown
        with vuetify.VCard(flat=True, classes="pa-2"):
            with vuetify.VCardTitle(classes="text-subtitle-1 py-1"):
                html.Div("Visualization Style")
            
            with vuetify.VCardText(classes="py-1"):
                with vuetify.VRow(classes="pa-0", dense=True):
                    with html.Div(style="width: 100%;"):
                        # Convert styles to display-friendly format
                        items = [
                            {"text": style.replace("_", " ").title(), "value": style}
                            for style in visualization_styles
                        ]
                        
                        vuetify.VSelect(
                            v_model=("visualization_style", visualization_styles[0]),
                            items=items,
                            label="Select Style",
                            outlined=True,
                            dense=True,
                            hide_details=True,
                            on_change=lambda style: ctrl.update_visualization_style(style=style),
                        )
        
        # Separator
        vuetify.VDivider()
        
        # Color Mapping Dropdown
        with vuetify.VCard(flat=True, classes="pa-2"):
            with vuetify.VCardTitle(classes="text-subtitle-1 py-1"):
                html.Div("Color Mapping")
            
            with vuetify.VCardText(classes="py-1"):
                with vuetify.VRow(classes="pa-0", dense=True):
                    with html.Div(style="width: 100%;"):
                        # Convert color mappings to display-friendly format
                        items = [
                            {"text": mapping.replace("_", " ").title(), "value": mapping}
                            for mapping in color_mappings
                        ]
                        
                        vuetify.VSelect(
                            v_model=("color_mapping", color_mappings[0]),
                            items=items,
                            label="Select Color Mapping",
                            outlined=True,
                            dense=True,
                            hide_details=True,
                            on_change=lambda mapping: ctrl.update_color_mapping(mapping=mapping),
                        )
        
        # Separator
        vuetify.VDivider()
        
        # Info Panel
        with vuetify.VCard(flat=True, classes="pa-2"):
            with vuetify.VCardTitle(classes="text-subtitle-1 py-1"):
                html.Div("About")
            
            with vuetify.VCardText(classes="py-1"):
                html.P("This is a molecular visualization tool for PDB files.")
                html.P("Hover over atoms to see detailed information.")
                html.P("Use different visualization styles and color mappings to explore the molecule.")
