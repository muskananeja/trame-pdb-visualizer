/**
 * 3D PDB Visualizer using Three.js
 * This script provides 3D visualization of PDB files with interactivity:
 * - Rotation, zoom, and pan controls
 * - Atom labels and tooltips
 * - Different visualization styles and color mappings
 */

// Global variables for 3D visualizer
let container, camera, scene, renderer, controls;
let moleculeGroup, labelGroup;
let raycaster, mouse;
let hoveredAtom = null;
let selectedAtom = null;
let atomLabels = [];

// Current visualization settings
let currentSettings = {
    visualStyle: 'Ball and Stick',
    colorMapping: 'Atom',
    bondThreshold: 3.0,
    atomSizeMultiplier: 1.0
    // Removed showLabels option since we now use hover information
};

// Atom colors by element (CPK coloring)
const atomColors = {
    'C': 0x909090, // Carbon - gray
    'N': 0x3050F8, // Nitrogen - blue
    'O': 0xFF0D0D, // Oxygen - red
    'S': 0xFFFF30, // Sulfur - yellow
    'P': 0xFF8000, // Phosphorus - orange
    'H': 0xFFFFFF, // Hydrogen - white
    'F': 0x90E050, // Fluorine - light green
    'CL': 0x1FF01F, // Chlorine - green
    'BR': 0xA62929, // Bromine - brown
    'I': 0x940094,  // Iodine - purple
    'CA': 0x00FF00, // Calcium - bright green
    'MG': 0x00FFFF, // Magnesium - cyan
    'NA': 0x0000FF, // Sodium - blue
    'K': 0x8000FF,  // Potassium - purple
    'ZN': 0x808080, // Zinc - gray
    'FE': 0xFFA500  // Iron - orange
};

// Element names (for display)
const elementNames = {
    'C': 'Carbon',
    'N': 'Nitrogen',
    'O': 'Oxygen',
    'S': 'Sulfur',
    'P': 'Phosphorus',
    'H': 'Hydrogen',
    'F': 'Fluorine',
    'CL': 'Chlorine',
    'BR': 'Bromine',
    'I': 'Iodine',
    'CA': 'Calcium',
    'MG': 'Magnesium',
    'NA': 'Sodium',
    'K': 'Potassium',
    'ZN': 'Zinc',
    'FE': 'Iron'
};

// Atom radii in Ångström scale (van der Waals)
const atomRadii = {
    'H': 1.2,
    'C': 1.7,
    'N': 1.55,
    'O': 1.52,
    'S': 1.8,
    'P': 1.8,
    'F': 1.47,
    'CL': 1.75,
    'BR': 1.85,
    'I': 1.98,
    'CA': 2.31,
    'MG': 1.73,
    'NA': 2.27,
    'K': 2.75,
    'ZN': 1.39,
    'FE': 1.56
};

// Residue colors
const residueColors = [
    0x33FF33, // Green
    0x3333FF, // Blue
    0xFF3333, // Red
    0x33FFFF, // Cyan
    0xFF33FF  // Magenta
];

/**
 * Initialize the 3D visualization environment
 * @param {string} containerId - The ID of the container element
 */
function init3DVisualizer(containerId) {
    // Get the container
    container = document.getElementById(containerId);
    if (!container) {
        console.error('Container not found:', containerId);
        return;
    }
    
    // Create the scene
    scene = new THREE.Scene();
    scene.background = new THREE.Color(0xf0f0f0);
    
    // Add lights
    const ambientLight = new THREE.AmbientLight(0x404040);
    scene.add(ambientLight);
    
    const directionalLight1 = new THREE.DirectionalLight(0xffffff, 0.8);
    directionalLight1.position.set(1, 1, 1);
    scene.add(directionalLight1);
    
    const directionalLight2 = new THREE.DirectionalLight(0xffffff, 0.5);
    directionalLight2.position.set(-1, -1, -1);
    scene.add(directionalLight2);
    
    // Create camera
    const width = container.clientWidth;
    const height = container.clientHeight;
    camera = new THREE.PerspectiveCamera(45, width / height, 0.1, 2000);
    camera.position.set(0, 0, 100);
    
    // Create renderer
    renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(width, height);
    container.appendChild(renderer.domElement);
    
    // Add OrbitControls
    controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;
    controls.dampingFactor = 0.25;
    
    // Groups to organize the scene
    moleculeGroup = new THREE.Group();
    labelGroup = new THREE.Group();
    scene.add(moleculeGroup);
    scene.add(labelGroup);
    
    // Setup for raycasting (for interactivity)
    raycaster = new THREE.Raycaster();
    mouse = new THREE.Vector2();
    
    // Add event listeners
    window.addEventListener('resize', onWindowResize);
    renderer.domElement.addEventListener('mousemove', onMouseMove);
    renderer.domElement.addEventListener('click', onMouseClick);
    
    // Start animation loop
    animate();
}

/**
 * Load and visualize a PDB molecule
 * @param {string} pdbContent - The content of a PDB file
 * @param {Object} options - Visualization options
 */
function visualizeMolecule3D(pdbContent, options = {}) {
    // Clear previous molecule
    while (moleculeGroup.children.length > 0) {
        const object = moleculeGroup.children[0];
        object.geometry.dispose();
        object.material.dispose();
        moleculeGroup.remove(object);
    }
    
    // Clear previous labels
    while (labelGroup.children.length > 0) {
        const label = labelGroup.children[0];
        labelGroup.remove(label);
    }
    atomLabels = [];
    
    // Update settings with provided options
    Object.assign(currentSettings, options);
    
    // Parse PDB content
    const atoms = parsePdbAtoms(pdbContent);
    console.log(`Parsed ${atoms.length} atoms for 3D rendering`);
    
    // Calculate molecule center to center the camera
    const center = calculateMoleculeCenter(atoms);
    
    // Calculate bonds
    const bonds = calculateBonds(atoms, currentSettings.bondThreshold);
    console.log(`Calculated ${bonds.length} bonds`);
    
    // Create 3D objects based on visual style
    if (currentSettings.visualStyle === 'Ball and Stick') {
        createBallAndStickModel(atoms, bonds, center);
    } else if (currentSettings.visualStyle === 'Protein Ribbon') {
        createProteinRibbonModel(atoms, bonds, center);
    } else {
        createSpaceFillModel(atoms, center);
    }
    
    // Reposition camera to view the whole molecule
    positionCameraForMolecule(atoms);
}

/**
 * Parse PDB file content into atom objects
 * @param {string} pdbContent - Content of the PDB file
 * @returns {Array} Array of atom objects with coordinates and properties
 */
function parsePdbAtoms(pdbContent) {
    const atoms = [];
    // Handle both Unix and Windows line endings and properly handle escaped newlines
    const lines = pdbContent.replace(/\\n/g, '\n').split(/\r?\n/);
    let residueIndex = 0;
    let currentResidue = null;
    
    for (const line of lines) {
        if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
            try {
                const x = parseFloat(line.substring(30, 38));
                const y = parseFloat(line.substring(38, 46));
                const z = parseFloat(line.substring(46, 54));
                const atomName = line.substring(12, 16).trim();
                const residueName = line.substring(17, 20).trim();
                const chainId = line.substring(21, 22).trim();
                const residueNumber = parseInt(line.substring(22, 26).trim());
                // Get element from columns 76-78 if available, otherwise use first character of atom name
                let element = line.substring(76, 78).trim();
                if (!element) {
                    // Extract first letter as the element (usually correct for common atoms)
                    element = atomName.replace(/[^A-Za-z]/g, '')[0];
                    // Convert to uppercase for consistency
                    element = element ? element.toUpperCase() : 'C'; // Default to Carbon if unknown
                }
                const bFactor = parseFloat(line.substring(60, 66).trim() || '0');
                
                // Assign residue index (for coloring)
                if (currentResidue !== residueNumber) {
                    currentResidue = residueNumber;
                    residueIndex++;
                }
                
                atoms.push({
                    x: x,
                    y: y,
                    z: z,
                    name: atomName,
                    element: element,
                    residue: {
                        name: residueName,
                        number: residueNumber,
                        index: residueIndex % 5,  // For coloring
                        chain: chainId
                    },
                    bFactor: bFactor
                });
            } catch (e) {
                console.error('Error parsing PDB line:', e, line);
            }
        }
    }
    
    console.log(`Parsed ${atoms.length} atoms for 3D rendering`);
    return atoms;
}

/**
 * Calculate bonds between atoms based on distance
 * @param {Array} atoms - Array of atom objects
 * @param {number} threshold - Distance threshold for bonding in Ångströms
 * @returns {Array} Array of bond objects
 */
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

/**
 * Calculate the center point of a molecule
 * @param {Array} atoms - Array of atom objects
 * @returns {Object} Center coordinates {x, y, z}
 */
function calculateMoleculeCenter(atoms) {
    let minX = Infinity, maxX = -Infinity;
    let minY = Infinity, maxY = -Infinity;
    let minZ = Infinity, maxZ = -Infinity;
    
    // Find min/max coordinates
    for (const atom of atoms) {
        minX = Math.min(minX, atom.x);
        maxX = Math.max(maxX, atom.x);
        minY = Math.min(minY, atom.y);
        maxY = Math.max(maxY, atom.y);
        minZ = Math.min(minZ, atom.z);
        maxZ = Math.max(maxZ, atom.z);
    }
    
    // Calculate center
    return {
        x: (minX + maxX) / 2,
        y: (minY + maxY) / 2,
        z: (minZ + maxZ) / 2
    };
}

/**
 * Create a ball-and-stick model for the molecule
 * @param {Array} atoms - Array of atom objects
 * @param {Array} bonds - Array of bond objects
 * @param {Object} center - Center coordinates {x, y, z}
 */
function createBallAndStickModel(atoms, bonds, center) {
    // Create spheres for atoms
    for (let i = 0; i < atoms.length; i++) {
        const atom = atoms[i];
        const radius = (atomRadii[atom.element] || 1.5) * 0.4 * currentSettings.atomSizeMultiplier;
        
        // Set color based on selected mapping
        let color;
        switch (currentSettings.colorMapping) {
            case 'Atom':
                color = atomColors[atom.element] || 0x808080;
                break;
            case 'B-Factor':
                // Normalize B-factor for this dataset
                const bFactorNorm = normalizeBFactor(atom.bFactor, atoms);
                color = getBFactorColor(bFactorNorm);
                break;
            case 'Residue':
                color = residueColors[atom.residue.index];
                break;
            default:
                color = atomColors[atom.element] || 0x808080;
        }
        
        const geometry = new THREE.SphereGeometry(radius, 16, 16);
        const material = new THREE.MeshPhongMaterial({ color: color });
        const sphere = new THREE.Mesh(geometry, material);
        
        // Position relative to center
        sphere.position.set(
            atom.x - center.x,
            atom.y - center.y,
            atom.z - center.z
        );
        
        // Add metadata to the sphere
        sphere.userData = {
            type: 'atom',
            index: i,
            atom: atom
        };
        
        moleculeGroup.add(sphere);
        
        // Add invisible text label
        addAtomLabel(atom, sphere.position, i);
    }
    
    // Create cylinders for bonds
    for (const bond of bonds) {
        const atom1 = atoms[bond.atom1];
        const atom2 = atoms[bond.atom2];
        
        const start = new THREE.Vector3(
            atom1.x - center.x,
            atom1.y - center.y,
            atom1.z - center.z
        );
        
        const end = new THREE.Vector3(
            atom2.x - center.x,
            atom2.y - center.y,
            atom2.z - center.z
        );
        
        const direction = new THREE.Vector3().subVectors(end, start);
        const length = direction.length();
        
        // Create a thin cylinder
        const bondGeometry = new THREE.CylinderGeometry(0.1, 0.1, length, 8);
        bondGeometry.translate(0, length / 2, 0);
        
        const bondMaterial = new THREE.MeshPhongMaterial({ color: 0xcccccc });
        const bondCylinder = new THREE.Mesh(bondGeometry, bondMaterial);
        
        // Position and rotate the cylinder to connect the atoms
        bondCylinder.position.copy(start);
        bondCylinder.lookAt(end);
        
        moleculeGroup.add(bondCylinder);
    }
}

/**
 * Create a space-fill model for the molecule
 * @param {Array} atoms - Array of atom objects
 * @param {Object} center - Center coordinates {x, y, z}
 */
function createSpaceFillModel(atoms, center) {
    // Create spheres for atoms with full van der Waals radii
    for (let i = 0; i < atoms.length; i++) {
        const atom = atoms[i];
        const radius = (atomRadii[atom.element] || 1.5) * currentSettings.atomSizeMultiplier;
        
        // Set color based on selected mapping
        let color;
        switch (currentSettings.colorMapping) {
            case 'Atom':
                color = atomColors[atom.element] || 0x808080;
                break;
            case 'B-Factor':
                const bFactorNorm = normalizeBFactor(atom.bFactor, atoms);
                color = getBFactorColor(bFactorNorm);
                break;
            case 'Residue':
                color = residueColors[atom.residue.index];
                break;
            default:
                color = atomColors[atom.element] || 0x808080;
        }
        
        const geometry = new THREE.SphereGeometry(radius, 32, 32);
        const material = new THREE.MeshPhongMaterial({ color: color });
        const sphere = new THREE.Mesh(geometry, material);
        
        // Position relative to center
        sphere.position.set(
            atom.x - center.x,
            atom.y - center.y,
            atom.z - center.z
        );
        
        // Add metadata to the sphere
        sphere.userData = {
            type: 'atom',
            index: i,
            atom: atom
        };
        
        moleculeGroup.add(sphere);
        
        // Add invisible text label
        addAtomLabel(atom, sphere.position, i);
    }
}

/**
 * Create a simplified ribbon model for protein visualization
 * @param {Array} atoms - Array of atom objects
 * @param {Array} bonds - Array of bond objects
 * @param {Object} center - Center coordinates {x, y, z}
 */
function createProteinRibbonModel(atoms, bonds, center) {
    // Extract alpha carbon atoms to trace the protein backbone
    const backboneAtoms = atoms.filter(atom => atom.name === 'CA');
    
    if (backboneAtoms.length < 2) {
        // Fall back to ball-and-stick if not enough backbone atoms
        console.warn('Not enough backbone atoms for ribbon model, using ball-and-stick');
        createBallAndStickModel(atoms, bonds, center);
        return;
    }
    
    // Create a smooth curve through alpha carbon positions
    const points = [];
    for (const atom of backboneAtoms) {
        points.push(new THREE.Vector3(
            atom.x - center.x,
            atom.y - center.y,
            atom.z - center.z
        ));
    }
    
    // Create a catmull-rom spline curve
    const curve = new THREE.CatmullRomCurve3(points);
    curve.curveType = 'catmullrom';
    curve.tension = 0.5;
    
    // Create a tube geometry along the curve
    const tubeGeometry = new THREE.TubeGeometry(curve, 
        Math.min(100, points.length * 4), // Segments
        0.4, // Radius
        8, // Radial segments
        false // Closed
    );
    
    // Create material based on color mapping
    let material;
    if (currentSettings.colorMapping === 'Residue') {
        // Use vertex colors for residue-based coloring
        const colors = [];
        const positions = tubeGeometry.attributes.position.array;
        const pointCount = positions.length / 3;
        
        for (let i = 0; i < pointCount; i++) {
            // Map position along tube to corresponding residue
            const t = i / pointCount;
            const index = Math.floor(t * (backboneAtoms.length - 1));
            const residueIndex = backboneAtoms[index].residue.index;
            const color = new THREE.Color(residueColors[residueIndex]);
            
            colors.push(color.r, color.g, color.b);
        }
        
        tubeGeometry.setAttribute('color', 
            new THREE.Float32BufferAttribute(colors, 3));
            
        material = new THREE.MeshPhongMaterial({ 
            vertexColors: true,
            shininess: 30
        });
    } else if (currentSettings.colorMapping === 'B-Factor') {
        // Use gradient material for B-factor coloring
        const bFactors = backboneAtoms.map(atom => atom.bFactor);
        const texture = createBFactorGradientTexture(bFactors);
        
        material = new THREE.MeshPhongMaterial({ 
            map: texture,
            shininess: 30
        });
    } else {
        // Default to a simple color
        material = new THREE.MeshPhongMaterial({ 
            color: 0x3050F8,
            shininess: 30
        });
    }
    
    const tube = new THREE.Mesh(tubeGeometry, material);
    moleculeGroup.add(tube);
    
    // Add key atoms (like CA) as small spheres
    for (let i = 0; i < backboneAtoms.length; i++) {
        const atom = backboneAtoms[i];
        const sphere = new THREE.Mesh(
            new THREE.SphereGeometry(0.3, 8, 8),
            new THREE.MeshPhongMaterial({ color: 0x909090 })
        );
        
        sphere.position.set(
            atom.x - center.x,
            atom.y - center.y,
            atom.z - center.z
        );
        
        // Add metadata to the sphere
        sphere.userData = {
            type: 'atom',
            index: atoms.indexOf(atom),
            atom: atom
        };
        
        moleculeGroup.add(sphere);
        
        // Add invisible text label
        addAtomLabel(atom, sphere.position, atoms.indexOf(atom));
    }
}

/**
 * Add an invisible text label for an atom
 * @param {Object} atom - Atom object
 * @param {THREE.Vector3} position - Position of the atom
 * @param {number} index - Index of the atom
 */
function addAtomLabel(atom, position, index) {
    // Create a div element for the label
    const labelDiv = document.createElement('div');
    labelDiv.className = 'atom-label';
    labelDiv.textContent = `${atom.element} - ${atom.name}`;
    labelDiv.style.display = 'none';
    
    // Add to document
    document.body.appendChild(labelDiv);
    
    // Store label information
    atomLabels.push({
        element: labelDiv,
        position: position.clone(),
        index: index,
        atom: atom
    });
}

/**
 * Show atom label at a specific position
 * @param {number} index - Index of the atom
 * @param {number} clientX - X position in screen coordinates
 * @param {number} clientY - Y position in screen coordinates
 */
function showAtomLabel(index, clientX, clientY) {
    // Hide any visible labels
    hideAllLabels();
    
    // Find the label for this atom
    const label = atomLabels.find(l => l.index === index);
    if (!label) return;
    
    // Update label content with detailed info
    const atom = label.atom;
    label.element.innerHTML = `
        <div class="atom-label-content">
            <span class="atom-element">${atom.element}</span>
            <span class="atom-name">${atom.name}</span>
            <div class="atom-details">
                Residue: ${atom.residue.name} ${atom.residue.number}
                <br>Chain: ${atom.residue.chain}
                <br>B-Factor: ${atom.bFactor.toFixed(2)}
                <br>Position: ${atom.x.toFixed(1)}, ${atom.y.toFixed(1)}, ${atom.z.toFixed(1)}
            </div>
        </div>
    `;
    
    // Position the label near the mouse pointer
    label.element.style.display = 'block';
    label.element.style.position = 'absolute';
    label.element.style.left = (clientX + 15) + 'px';
    label.element.style.top = (clientY - 15) + 'px';
}

/**
 * Hide all atom labels
 */
function hideAllLabels() {
    for (const label of atomLabels) {
        label.element.style.display = 'none';
    }
}

/**
 * Position camera to view the whole molecule
 * @param {Array} atoms - Array of atom objects
 */
function positionCameraForMolecule(atoms) {
    let minX = Infinity, maxX = -Infinity;
    let minY = Infinity, maxY = -Infinity;
    let minZ = Infinity, maxZ = -Infinity;
    
    // Find min/max coordinates
    for (const atom of atoms) {
        minX = Math.min(minX, atom.x);
        maxX = Math.max(maxX, atom.x);
        minY = Math.min(minY, atom.y);
        maxY = Math.max(maxY, atom.y);
        minZ = Math.min(minZ, atom.z);
        maxZ = Math.max(maxZ, atom.z);
    }
    
    // Calculate the bounding sphere radius
    const center = new THREE.Vector3(
        (minX + maxX) / 2,
        (minY + maxY) / 2,
        (minZ + maxZ) / 2
    );
    
    const radius = Math.max(
        maxX - minX,
        maxY - minY,
        maxZ - minZ
    ) / 2;
    
    // Set camera position based on bounding sphere
    const fov = camera.fov * (Math.PI / 180);
    const distanceToFit = (radius * 1.5) / Math.sin(fov / 2);
    
    camera.position.set(0, 0, distanceToFit);
    controls.target.set(0, 0, 0);
    controls.update();
}

/**
 * Handle window resize events
 */
function onWindowResize() {
    const width = container.clientWidth;
    const height = container.clientHeight;
    
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
    
    renderer.setSize(width, height);
}

/**
 * Handle mouse move events for hover effects
 * @param {Event} event - Mouse event
 */
function onMouseMove(event) {
    // Calculate mouse position in normalized device coordinates (-1 to +1)
    const rect = renderer.domElement.getBoundingClientRect();
    mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
    mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
    
    // Check for atom hover
    raycaster.setFromCamera(mouse, camera);
    const intersects = raycaster.intersectObjects(moleculeGroup.children);
    
    if (intersects.length > 0) {
        const object = intersects[0].object;
        if (object.userData && object.userData.type === 'atom') {
            if (hoveredAtom !== object.userData.index) {
                // New atom hovered
                hoveredAtom = object.userData.index;
                showAtomLabel(hoveredAtom, event.clientX, event.clientY);
                
                // Highlight the atom
                object.material.emissive.set(0x555555);
                document.body.style.cursor = 'pointer';
            }
        }
    } else if (hoveredAtom !== null) {
        // No longer hovering over an atom
        const previousObj = moleculeGroup.children.find(
            obj => obj.userData && obj.userData.type === 'atom' && 
                  obj.userData.index === hoveredAtom
        );
        
        if (previousObj) {
            previousObj.material.emissive.set(0x000000);
        }
        
        hoveredAtom = null;
        hideAllLabels();
        document.body.style.cursor = 'default';
    }
}

/**
 * Handle mouse click events for selecting atoms
 * @param {Event} event - Mouse event
 */
function onMouseClick(event) {
    // Use the same raycasting logic as mousemove
    const rect = renderer.domElement.getBoundingClientRect();
    mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
    mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
    
    raycaster.setFromCamera(mouse, camera);
    const intersects = raycaster.intersectObjects(moleculeGroup.children);
    
    if (intersects.length > 0) {
        const object = intersects[0].object;
        if (object.userData && object.userData.type === 'atom') {
            const atomIndex = object.userData.index;
            
            // If this atom is already selected, deselect it
            if (selectedAtom === atomIndex) {
                selectedAtom = null;
                object.material.color.set(
                    getAtomColor(object.userData.atom, currentSettings.colorMapping)
                );
            } else {
                // Deselect previous atom if any
                if (selectedAtom !== null) {
                    const previousObj = moleculeGroup.children.find(
                        obj => obj.userData && obj.userData.type === 'atom' && 
                              obj.userData.index === selectedAtom
                    );
                    
                    if (previousObj) {
                        previousObj.material.color.set(
                            getAtomColor(previousObj.userData.atom, currentSettings.colorMapping)
                        );
                    }
                }
                
                // Select new atom
                selectedAtom = atomIndex;
                object.material.color.set(0xFFFF00); // Highlight in yellow
                
                // Display information in a more permanent way
                const atom = object.userData.atom;
                const atomInfo = `Element: ${atom.element} (${elementNames[atom.element] || 'Unknown'})<br>
                               Name: ${atom.name}<br>
                               Residue: ${atom.residue.name} ${atom.residue.number}<br>
                               Chain: ${atom.residue.chain}<br>
                               B-Factor: ${atom.bFactor.toFixed(2)}<br>
                               Position: ${atom.x.toFixed(2)}, ${atom.y.toFixed(2)}, ${atom.z.toFixed(2)}`;
                
                // Update the atom information display
                const infoElement = document.getElementById('atom-info-display');
                if (infoElement) {
                    infoElement.innerHTML = atomInfo;
                    infoElement.style.display = 'block';
                }
            }
        }
    }
}

/**
 * Animation loop
 */
function animate() {
    requestAnimationFrame(animate);
    
    controls.update();
    renderer.render(scene, camera);
    
    // Update any labels that should follow atoms
    updateAtomLabelsPosition();
}

/**
 * Update the positions of visible atom labels to follow their atoms
 */
function updateAtomLabelsPosition() {
    for (const label of atomLabels) {
        if (label.element.style.display === 'block') {
            // Project the 3D position to screen coordinates
            const position = label.position.clone();
            position.project(camera);
            
            // Convert to screen coordinates
            const x = (position.x * 0.5 + 0.5) * container.clientWidth;
            const y = (1 - (position.y * 0.5 + 0.5)) * container.clientHeight;
            
            // Update label position
            label.element.style.left = (x + 15) + 'px';
            label.element.style.top = (y - 15) + 'px';
        }
    }
}

/**
 * Get atom color based on the selected color mapping
 * @param {Object} atom - Atom object
 * @param {string} colorMapping - The color mapping method
 * @returns {number} - Color in hexadecimal format
 */
function getAtomColor(atom, colorMapping) {
    switch (colorMapping) {
        case 'Atom':
            return atomColors[atom.element] || 0x808080;
        case 'B-Factor':
            const bFactorNorm = normalizeBFactor(atom.bFactor, [atom]);
            return getBFactorColor(bFactorNorm);
        case 'Residue':
            return residueColors[atom.residue.index];
        default:
            return atomColors[atom.element] || 0x808080;
    }
}

/**
 * Normalize B-factor values for coloring
 * @param {number} bFactor - B-factor value
 * @param {Array} atoms - Array of all atoms for min/max calculation
 * @returns {number} - Normalized value between 0 and 1
 */
function normalizeBFactor(bFactor, atoms) {
    // Find min and max B-factors
    let minB = Infinity, maxB = -Infinity;
    for (const atom of atoms) {
        minB = Math.min(minB, atom.bFactor);
        maxB = Math.max(maxB, atom.bFactor);
    }
    
    // Normalize to range [0,1]
    if (maxB === minB) return 0.5; // Default to middle if all values are the same
    return (bFactor - minB) / (maxB - minB);
}

/**
 * Get color for B-factor visualization
 * @param {number} value - Normalized B-factor value (0-1)
 * @returns {number} - Color in hexadecimal format
 */
function getBFactorColor(value) {
    // Cold (blue) to hot (red) gradient
    const r = Math.floor(255 * value);
    const g = Math.floor(255 * (1 - Math.abs(2 * value - 1)));
    const b = Math.floor(255 * (1 - value));
    
    return (r << 16) | (g << 8) | b;
}

/**
 * Create a gradient texture for B-factor visualization
 * @param {Array} bFactors - Array of B-factor values
 * @returns {THREE.Texture} - The generated texture
 */
function createBFactorGradientTexture(bFactors) {
    const canvas = document.createElement('canvas');
    canvas.width = 256;
    canvas.height = 1;
    
    const ctx = canvas.getContext('2d');
    const gradient = ctx.createLinearGradient(0, 0, 256, 0);
    
    gradient.addColorStop(0, 'blue');
    gradient.addColorStop(0.5, 'green');
    gradient.addColorStop(1, 'red');
    
    ctx.fillStyle = gradient;
    ctx.fillRect(0, 0, 256, 1);
    
    const texture = new THREE.CanvasTexture(canvas);
    texture.needsUpdate = true;
    
    return texture;
}

/**
 * Toggle atom labels display - used to be for the Show Labels toggle, 
 * now adapted to work with hover functionality
 * @param {boolean} show - If true, allows visibility on hover
 */
function toggleLabels(show) {
    // We now handle all label visibility through the hover functionality only
    // This function is kept for backward compatibility
    
    // Hide all labels except for hovered atoms
    for (const label of atomLabels) {
        if (label.index !== hoveredAtom) {
            label.element.style.display = 'none';
        }
    }
}

/**
 * Update visualization style and color mapping
 * @param {Object} options - Visualization options
 */
function updateVisualization(options) {
    // Get the current PDB content from the Vue app
    const appElement = document.getElementById('app');
    if (!appElement || !appElement.__vue__) {
        console.error('Vue app not found');
        return;
    }
    
    const pdbContent = appElement.__vue__.pdbContent;
    if (!pdbContent) {
        console.error('No PDB content available');
        return;
    }
    
    // Visualize with new options
    visualizeMolecule3D(pdbContent, options);
}