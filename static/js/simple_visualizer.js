/**
 * A simpler 2D PDB Visualizer for the molecular structures
 * 
 * Features:
 * - Renders atoms with proper colors
 * - Shows bonds between nearby atoms
 * - Provides top view (X-Y) and side view (X-Z) projections
 * - Handles window resizing for responsive display
 */

// Global variables for responsive visualization
let currentTopAtoms = [];
let currentSideAtoms = [];
let currentBondsTop = [];
let currentBondsSide = [];
let currentTopOptions = {};
let currentSideOptions = {};
let canResize = false;

// Responsive visualization - redraw when window is resized
window.addEventListener('resize', function() {
    if (canResize && currentTopAtoms.length > 0) {
        drawMoleculeTopView('topViewCanvas', currentTopAtoms, currentBondsTop, currentTopOptions);
        drawMoleculeSideView('sideViewCanvas', currentSideAtoms, currentBondsSide, currentSideOptions);
    }
});

// Function to create a static top view of the molecule
function createTopView(pdbData, canvasId, options = {}) {
    // Parse the PDB data to extract coordinates
    const atoms = parsePdbAtoms(pdbData);
    console.log(`Parsed ${atoms.length} atoms for static rendering`);
    
    // Store for responsive redrawing
    currentTopAtoms = atoms;
    
    // Apply options
    const bondThreshold = options.bondThreshold || 3.0; // Default 3Å threshold for bonds
    const atomSizeMultiplier = options.atomSizeMultiplier || 1.0; // Default size multiplier
    const colorMapping = options.colorMapping || 'Atom'; // Default color mapping
    
    // Calculate bonds between atoms with the specified threshold
    const bonds = calculateBonds(atoms, bondThreshold);
    currentBondsTop = bonds;
    
    // Store the options for responsive redrawing
    currentTopOptions = {
        bondThreshold: bondThreshold,
        atomSizeMultiplier: atomSizeMultiplier,
        colorMapping: colorMapping
    };
    
    // Draw the molecule
    drawMoleculeTopView(canvasId, atoms, bonds, currentTopOptions);
    
    // Enable responsive redrawing
    canResize = true;
}

// Function to draw the molecule in top view (separating calculation from drawing)
function drawMoleculeTopView(canvasId, atoms, bonds, options = {}) {
    // Get the canvas and drawing context
    const canvas = document.getElementById(canvasId);
    if (!canvas) {
        console.error('Canvas not found:', canvasId);
        return;
    }
    
    const ctx = canvas.getContext('2d');
    if (!ctx) {
        console.error('Could not get canvas context');
        return;
    }
    
    // Apply options
    const atomSizeMultiplier = options.atomSizeMultiplier || 1.0;
    const colorMapping = options.colorMapping || 'Atom';
    const visualStyle = options.visualStyle || 'Ball and Stick';
    
    // Clear canvas and set background
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = '#f0f0f0';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    
    // Draw title with visualization style and color mapping
    ctx.fillStyle = '#000000';
    ctx.font = '14px Arial';
    ctx.textAlign = 'center';
    ctx.fillText('Top View (X-Y) - ' + visualStyle + ' / ' + colorMapping, canvas.width / 2, 20);
    
    // Calculate the bounds of the molecule
    let minX = Infinity, maxX = -Infinity;
    let minY = Infinity, maxY = -Infinity;
    
    for (const atom of atoms) {
        minX = Math.min(minX, atom.x);
        maxX = Math.max(maxX, atom.x);
        minY = Math.min(minY, atom.y);
        maxY = Math.max(maxY, atom.y);
    }
    
    // Calculate scale to fit canvas (with some padding)
    const padding = 40;
    const scaleX = (canvas.width - padding * 2) / (maxX - minX);
    const scaleY = (canvas.height - padding * 2) / (maxY - minY);
    const scale = Math.min(scaleX, scaleY) * 0.8;
    
    // Store transformed coordinates for bonds
    const transformedCoords = {};
    
    // First: Draw bonds - thicker for "Ball and Stick", thinner for "Protein Ribbon"
    ctx.lineWidth = visualStyle === 'Ball and Stick' ? 1.5 : 0.8;
    ctx.strokeStyle = '#cccccc';
    
    for (const bond of bonds) {
        const atom1 = atoms[bond.atom1];
        const atom2 = atoms[bond.atom2];
        
        const x1 = (atom1.x - minX) * scale + padding;
        const y1 = (atom1.y - minY) * scale + padding;
        const x2 = (atom2.x - minX) * scale + padding;
        const y2 = (atom2.y - minY) * scale + padding;
        
        // Store coords for atom drawing
        transformedCoords[bond.atom1] = {x: x1, y: y1};
        transformedCoords[bond.atom2] = {x: x2, y: y2};
        
        // Draw the bond
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();
    }
    
    // Second: Draw atoms on top of bonds
    for (let i = 0; i < atoms.length; i++) {
        const atom = atoms[i];
        
        // Use precomputed coordinates if available (from bond drawing)
        let x, y;
        if (transformedCoords[i]) {
            x = transformedCoords[i].x;
            y = transformedCoords[i].y;
        } else {
            x = (atom.x - minX) * scale + padding;
            y = (atom.y - minY) * scale + padding;
        }
        
        // Set base circle size based on atom type and apply style multiplier
        let radius = 4;  // Base size
        if (atom.element === 'C') radius = 4;
        else if (atom.element === 'N') radius = 4.5;
        else if (atom.element === 'O') radius = 4.5;
        else if (atom.element === 'S') radius = 5;
        else if (atom.element === 'P') radius = 5;
        else radius = 4;
        
        // Apply size multiplier from visualization style
        radius *= atomSizeMultiplier;
        
        // Apply color based on selected color mapping
        let color;
        switch (colorMapping) {
            case 'Atom':
                // Color by atom element
                color = getAtomColor(atom.element);
                break;
            case 'B-Factor':
                // Simple "mock" B-factor coloring - in reality would use actual B-factor values
                const zNormalized = (atom.z - Math.min(...atoms.map(a => a.z))) / 
                                   (Math.max(...atoms.map(a => a.z)) - Math.min(...atoms.map(a => a.z)));
                color = getBFactorColor(zNormalized);
                break;
            case 'Residue':
                // Simple mock residue coloring - alternating colors by index
                color = getResidueColor(i % 5);
                break;
            default:
                color = getAtomColor(atom.element);
        }
        ctx.fillStyle = color;
        
        // Draw the atom
        ctx.beginPath();
        ctx.arc(x, y, radius, 0, Math.PI * 2);
        ctx.fill();
        
        // Add a slight outline for better visibility
        ctx.strokeStyle = 'rgba(0,0,0,0.3)';
        ctx.lineWidth = 0.5;
        ctx.stroke();
    }
}

// Function to create a static side view of the molecule
function createSideView(pdbData, canvasId, options = {}) {
    // Parse the PDB data to extract coordinates
    const atoms = parsePdbAtoms(pdbData);
    console.log(`Parsed ${atoms.length} atoms for side view rendering`);
    
    // Store for responsive redrawing
    currentSideAtoms = atoms;
    
    // Apply options
    const bondThreshold = options.bondThreshold || 3.0; // Default 3Å threshold for bonds
    const atomSizeMultiplier = options.atomSizeMultiplier || 1.0; // Default size multiplier
    const colorMapping = options.colorMapping || 'Atom'; // Default color mapping
    const visualStyle = options.visualStyle || 'Ball and Stick'; // Default style
    
    // Calculate bonds between atoms with the specified threshold
    const bonds = calculateBonds(atoms, bondThreshold);
    currentBondsSide = bonds;
    
    // Store the options for responsive redrawing
    currentSideOptions = {
        bondThreshold: bondThreshold,
        atomSizeMultiplier: atomSizeMultiplier,
        colorMapping: colorMapping,
        visualStyle: visualStyle
    };
    
    // Draw the molecule
    drawMoleculeSideView(canvasId, atoms, bonds, currentSideOptions);
}

// Function to draw the molecule in side view (separating calculation from drawing)
function drawMoleculeSideView(canvasId, atoms, bonds, options = {}) {
    // Get the canvas and drawing context
    const canvas = document.getElementById(canvasId);
    if (!canvas) {
        console.error('Canvas not found:', canvasId);
        return;
    }
    
    const ctx = canvas.getContext('2d');
    if (!ctx) {
        console.error('Could not get canvas context');
        return;
    }
    
    // Apply options
    const atomSizeMultiplier = options.atomSizeMultiplier || 1.0;
    const colorMapping = options.colorMapping || 'Atom';
    const visualStyle = options.visualStyle || 'Ball and Stick';
    
    // Clear canvas and set background
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.fillStyle = '#f0f0f0';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    
    // Draw title with visualization style and color mapping
    ctx.fillStyle = '#000000';
    ctx.font = '14px Arial';
    ctx.textAlign = 'center';
    ctx.fillText('Side View (X-Z) - ' + visualStyle + ' / ' + colorMapping, canvas.width / 2, 20);
    
    // Calculate the bounds of the molecule
    let minX = Infinity, maxX = -Infinity;
    let minZ = Infinity, maxZ = -Infinity;
    
    for (const atom of atoms) {
        minX = Math.min(minX, atom.x);
        maxX = Math.max(maxX, atom.x);
        minZ = Math.min(minZ, atom.z);
        maxZ = Math.max(maxZ, atom.z);
    }
    
    // Calculate scale to fit canvas (with some padding)
    const padding = 40;
    const scaleX = (canvas.width - padding * 2) / (maxX - minX);
    const scaleZ = (canvas.height - padding * 2) / (maxZ - minZ);
    const scale = Math.min(scaleX, scaleZ) * 0.8;
    
    // Store transformed coordinates for bonds
    const transformedCoords = {};
    
    // First: Draw bonds - thicker for "Ball and Stick", thinner for "Protein Ribbon"
    ctx.lineWidth = visualStyle === 'Ball and Stick' ? 1.5 : 0.8;
    ctx.strokeStyle = '#cccccc';
    
    for (const bond of bonds) {
        const atom1 = atoms[bond.atom1];
        const atom2 = atoms[bond.atom2];
        
        const x1 = (atom1.x - minX) * scale + padding;
        const y1 = (atom1.z - minZ) * scale + padding;
        const x2 = (atom2.x - minX) * scale + padding;
        const y2 = (atom2.z - minZ) * scale + padding;
        
        // Store coords for atom drawing
        transformedCoords[bond.atom1] = {x: x1, y: y1};
        transformedCoords[bond.atom2] = {x: x2, y: y2};
        
        // Draw the bond
        ctx.beginPath();
        ctx.moveTo(x1, y1);
        ctx.lineTo(x2, y2);
        ctx.stroke();
    }
    
    // Second: Draw atoms on top of bonds
    for (let i = 0; i < atoms.length; i++) {
        const atom = atoms[i];
        
        // Use precomputed coordinates if available (from bond drawing)
        let x, y;
        if (transformedCoords[i]) {
            x = transformedCoords[i].x;
            y = transformedCoords[i].y;
        } else {
            x = (atom.x - minX) * scale + padding;
            y = (atom.z - minZ) * scale + padding;
        }
        
        // Set base circle size based on atom type
        let radius = 4;  // Base size
        if (atom.element === 'C') radius = 4;
        else if (atom.element === 'N') radius = 4.5;
        else if (atom.element === 'O') radius = 4.5;
        else if (atom.element === 'S') radius = 5;
        else if (atom.element === 'P') radius = 5;
        else radius = 4;
        
        // Apply size multiplier from visualization style
        radius *= atomSizeMultiplier;
        
        // Add small adjustment for depth effect (based on Y position)
        const depthFactor = 0.5 + 0.5 * (atom.y - Math.min(...atoms.map(a => a.y))) / 
                              (Math.max(...atoms.map(a => a.y)) - Math.min(...atoms.map(a => a.y)));
        radius *= depthFactor;
        
        // Apply color based on selected color mapping
        let color;
        switch (colorMapping) {
            case 'Atom':
                // Color by atom element
                color = getAtomColor(atom.element);
                break;
            case 'B-Factor':
                // Simple "mock" B-factor coloring - in reality would use actual B-factor values
                const yNormalized = (atom.y - Math.min(...atoms.map(a => a.y))) / 
                                   (Math.max(...atoms.map(a => a.y)) - Math.min(...atoms.map(a => a.y)));
                color = getBFactorColor(yNormalized);
                break;
            case 'Residue':
                // Simple mock residue coloring - alternating colors by index
                color = getResidueColor(i % 5);
                break;
            default:
                color = getAtomColor(atom.element);
        }
        ctx.fillStyle = color;
        
        // Draw the atom
        ctx.beginPath();
        ctx.arc(x, y, radius, 0, Math.PI * 2);
        ctx.fill();
        
        // Add a slight outline for better visibility
        ctx.strokeStyle = 'rgba(0,0,0,0.3)';
        ctx.lineWidth = 0.5;
        ctx.stroke();
    }
}

// Function to calculate bonds between atoms based on distance
function calculateBonds(atoms, threshold) {
    const bonds = [];
    
    // Calculate bonds (simple distance-based approach)
    for (let i = 0; i < atoms.length; i++) {
        for (let j = i + 1; j < atoms.length; j++) {
            const atom1 = atoms[i];
            const atom2 = atoms[j];
            
            // Calculate 3D distance
            const dx = atom1.x - atom2.x;
            const dy = atom1.y - atom2.y;
            const dz = atom1.z - atom2.z;
            const distance = Math.sqrt(dx*dx + dy*dy + dz*dz);
            
            // If atoms are close enough, they're bonded
            if (distance < threshold) {
                bonds.push({
                    atom1: i,
                    atom2: j,
                    distance: distance
                });
            }
        }
    }
    
    return bonds;
}

// Helper function to parse PDB data
function parsePdbAtoms(pdbData) {
    const atoms = [];
    const lines = pdbData.split('\n');
    
    for (const line of lines) {
        if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
            try {
                const x = parseFloat(line.substring(30, 38));
                const y = parseFloat(line.substring(38, 46));
                const z = parseFloat(line.substring(46, 54));
                const atomName = line.substring(12, 16).trim();
                const element = line.substring(76, 78).trim() || atomName.charAt(0);
                
                atoms.push({
                    x: x,
                    y: y,
                    z: z,
                    element: element
                });
            } catch (e) {
                console.error('Error parsing PDB line:', e);
            }
        }
    }
    
    return atoms;
}

// Helper function to get color for an atom
function getAtomColor(element) {
    const atomColors = {
        'C': '#909090', // Carbon - gray
        'N': '#3050F8', // Nitrogen - blue
        'O': '#FF0D0D', // Oxygen - red
        'S': '#FFFF30', // Sulfur - yellow
        'P': '#FF8000', // Phosphorus - orange
        'H': '#FFFFFF', // Hydrogen - white
        'F': '#90E050', // Fluorine - light green
        'CL': '#1FF01F', // Chlorine - green
        'BR': '#A62929', // Bromine - brown
        'I': '#940094',  // Iodine - purple
        'CA': '#00FF00', // Calcium - bright green
        'MG': '#00FFFF', // Magnesium - cyan
        'NA': '#0000FF', // Sodium - blue
        'K': '#8000FF',  // Potassium - purple
        'ZN': '#808080', // Zinc - gray
        'FE': '#FFA500'  // Iron - orange
    };
    
    return atomColors[element] || '#808080'; // Default to gray
}

// Helper function to get color for B-factor
function getBFactorColor(value) {
    // Color gradient from blue (cold) to red (hot)
    const r = Math.floor(255 * value);
    const b = Math.floor(255 * (1 - value));
    return `rgb(${r}, 100, ${b})`;
}

// Helper function to get color for residue types
function getResidueColor(index) {
    // Different colors for different residue types
    const residueColors = [
        '#33FF33', // Green
        '#3333FF', // Blue
        '#FF3333', // Red
        '#33FFFF', // Cyan
        '#FF33FF'  // Magenta
    ];
    
    return residueColors[index % residueColors.length];
}