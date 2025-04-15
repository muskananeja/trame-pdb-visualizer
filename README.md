# Molecular Visualizer - 3D Interactive Structure Viewer

**Muskan Aneja's Capstone Project**

## Project Overview

Molecular Visualizer is an interactive web-based platform for visualizing protein structures from Protein Data Bank (PDB) files. This tool allows researchers, students, and molecular biology enthusiasts to explore molecular structures in both 2D and 3D representations with various visualization styles and color mapping options.

## Key Features

- **3D Interactive Visualization**: Rotate, zoom, and pan to explore molecular structures from any angle
- **Multiple Visualization Styles**: 
  - Ball and Stick model 
  - Protein Ribbon model
- **Diverse Color Mapping Options**:
  - Element-based coloring (CPK coloring)
  - B-factor (temperature) coloring
  - Residue-based coloring
- **Interactive Atom Information**: Hover over atoms to view detailed information (element, residue, chain, B-factor, coordinates)
- **2D Projections**: Top view (X-Y plane) and side view (X-Z plane) for complementary visualization
- **File Upload**: Upload and visualize your own PDB files
- **Responsive Design**: Works on desktop and mobile devices

## Technologies Used

- **Backend**: Flask (Python)
- **Frontend**: 
  - Three.js for 3D rendering
  - Vue.js and Vuetify for UI components
  - HTML5 Canvas for 2D visualization
- **Data Processing**: Custom PDB file parsing algorithms
- **Containerization**: Docker for easy deployment

## How to Run the Application

### Method 1: Using Python directly

#### Windows:
1. Make sure Python 3.6+ is installed
2. Double-click on `run.bat` or open a command prompt and run:
   ```
   run.bat
   ```
3. Open http://localhost:5000 in your browser

#### Mac/Linux:
1. Make sure Python 3.6+ is installed
2. Open a terminal and run:
   ```
   chmod +x run.sh
   ./run.sh
   ```
3. Open http://localhost:5000 in your browser

### Method 2: Using Docker

1. Make sure Docker is installed on your system
2. Build the Docker image:
   ```
   docker build -t molecular-visualizer .
   ```
3. Run the Docker container:
   ```
   docker run -p 5000:5000 molecular-visualizer
   ```
4. Open http://localhost:5000 in your browser

### Method 3: Using Python manually

1. Install the required dependencies:
   ```
   pip install -e .
   ```
2. Run the application:
   ```
   python app.py
   ```
3. Open http://localhost:5000 in your browser

## How to Use

1. Once the application is running, open http://localhost:5000 in your browser
2. Click "Load Example Structure" or upload your own PDB file
3. Use mouse controls to interact with the 3D view:
   - Left-click + drag: Rotate the molecule
   - Right-click + drag: Pan the view
   - Scroll wheel: Zoom in/out
4. Try different visualization styles and color mappings using the dropdown options
5. Hover over atoms to see detailed information
6. Click the fullscreen button for a more immersive experience

## Future Enhancements

- Support for more visualization styles (Space-filling, Wireframe)
- Molecular surface visualization
- Distance and angle measurement tools
- Ligand highlighting
- Animation capabilities for dynamics visualization
- Export options for images and 3D models

## License

Open source for educational and research purposes.

## Project Structure

```
/ (root)
├── static/                  # Static assets
│   ├── css/                 # CSS stylesheets
│   ├── js/                  # JavaScript files
│   │   ├── 3d_visualizer.js # Three.js based 3D visualization
│   │   └── simple_visualizer.js # 2D canvas visualization
│   ├── lib/                 # Third-party libraries
│   └── examples/            # Example PDB files
├── templates/               # HTML templates
│   └── index.html           # Main application page
├── app.py                   # Flask application
├── setup.py                 # Package setup and dependencies
├── run.bat                  # Windows run script
├── run.sh                   # Mac/Linux run script
└── Dockerfile               # Docker configuration file
```