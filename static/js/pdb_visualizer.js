/**
 * Simple 2D PDB Visualizer
 * 
 * This script creates a simple 2D representation of a protein from PDB data
 * It extracts the 3D coordinates and projects them onto a 2D canvas with a
 * simple orthographic projection
 */

class PDB2DVisualizer {
    constructor(canvasId) {
        console.log('Initializing PDB2DVisualizer for', canvasId);
        this.canvas = document.getElementById(canvasId);
        console.log('Canvas element:', this.canvas);
        if (!this.canvas) {
            console.error('Canvas element not found:', canvasId);
            return;
        }
        this.ctx = this.canvas.getContext('2d');
        console.log('Canvas context:', this.ctx);
        this.atomColors = {
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
            'FE': '#FFA500', // Iron - orange
            'default': '#808080' // Default - gray
        };
        
        this.reset();
    }
    
    reset() {
        this.atoms = [];
        this.xMin = Number.POSITIVE_INFINITY;
        this.xMax = Number.NEGATIVE_INFINITY;
        this.yMin = Number.POSITIVE_INFINITY;
        this.yMax = Number.NEGATIVE_INFINITY;
        this.zMin = Number.POSITIVE_INFINITY;
        this.zMax = Number.NEGATIVE_INFINITY;
    }
    
    parsePDB(pdbContent) {
        console.log('Parsing PDB data...');
        this.reset();
        
        const lines = pdbContent.split('\n');
        
        for (const line of lines) {
            if (line.startsWith('ATOM') || line.startsWith('HETATM')) {
                try {
                    const x = parseFloat(line.substring(30, 38));
                    const y = parseFloat(line.substring(38, 46));
                    const z = parseFloat(line.substring(46, 54));
                    const atomName = line.substring(12, 16).trim();
                    const atomType = atomName.charAt(0); // First letter is element type
                    
                    this.atoms.push({
                        x: x,
                        y: y,
                        z: z,
                        type: atomType,
                        radius: this.getAtomRadius(atomType)
                    });
                    
                    // Update bounds
                    this.xMin = Math.min(this.xMin, x);
                    this.xMax = Math.max(this.xMax, x);
                    this.yMin = Math.min(this.yMin, y);
                    this.yMax = Math.max(this.yMax, y);
                    this.zMin = Math.min(this.zMin, z);
                    this.zMax = Math.max(this.zMax, z);
                } catch (e) {
                    console.error('Error parsing PDB line:', e);
                }
            }
        }
        
        console.log(`Parsed ${this.atoms.length} atoms`);
        console.log(`X range: ${this.xMin} to ${this.xMax}`);
        console.log(`Y range: ${this.yMin} to ${this.yMax}`);
        console.log(`Z range: ${this.zMin} to ${this.zMax}`);
    }
    
    getAtomRadius(atomType) {
        // Simple atomic radii in Angstroms (scaled for visualization)
        const radii = {
            'C': 0.7,
            'N': 0.65,
            'O': 0.6,
            'S': 1.0,
            'P': 1.0,
            'H': 0.25,
            'F': 0.5,
            'CL': 1.0,
            'BR': 1.15,
            'I': 1.4,
            'default': 0.7
        };
        
        return radii[atomType] || radii['default'];
    }
    
    getAtomColor(atomType) {
        return this.atomColors[atomType] || this.atomColors['default'];
    }
    
    renderTopView() {
        this.render('top');
    }
    
    renderSideView() {
        this.render('side');
    }
    
    render(viewType = 'top') {
        console.log(`Rendering ${viewType} view with ${this.atoms.length} atoms`);
        
        if (this.atoms.length === 0) {
            console.error('No atoms to render');
            return;
        }
        
        console.log('Canvas dimensions:', this.canvas.width, this.canvas.height);
        
        // Clear canvas
        this.ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
        
        // Set background
        this.ctx.fillStyle = '#f5f5f5';
        this.ctx.fillRect(0, 0, this.canvas.width, this.canvas.height);
        
        // Calculate scale to fit canvas
        const padding = 40;
        let scaleX, scaleY;
        let centerX, centerY;
        
        // Choose projection based on view type
        let projection;
        if (viewType === 'top') {
            // Top view: x, y coordinates
            scaleX = (this.canvas.width - padding * 2) / (this.xMax - this.xMin);
            scaleY = (this.canvas.height - padding * 2) / (this.yMax - this.yMin);
            centerX = this.xMin;
            centerY = this.yMin;
            projection = (atom) => ({
                x: (atom.x - centerX) * scaleX + padding,
                y: (atom.y - centerY) * scaleY + padding,
                z: atom.z
            });
        } else if (viewType === 'side') {
            // Side view: x, z coordinates
            scaleX = (this.canvas.width - padding * 2) / (this.xMax - this.xMin);
            scaleY = (this.canvas.height - padding * 2) / (this.zMax - this.zMin);
            centerX = this.xMin;
            centerY = this.zMin;
            projection = (atom) => ({
                x: (atom.x - centerX) * scaleX + padding,
                y: (atom.z - centerY) * scaleY + padding,
                z: atom.y
            });
        } else {
            console.error('Invalid view type');
            return;
        }
        
        // Use the smaller scale to maintain aspect ratio
        const scale = Math.min(scaleX, scaleY) * 0.9;
        
        // Sort atoms by z-coordinate for depth effect
        const sortedAtoms = [...this.atoms].sort((a, b) => a.z - b.z);
        
        // Draw atoms
        for (const atom of sortedAtoms) {
            const projected = projection(atom);
            const x = (atom.x - centerX) * scale + this.canvas.width / 2;
            const y = (projected.y - centerY) * scale + this.canvas.height / 2;
            
            // Adjust radius based on z-coordinate for depth effect
            let depthFactor = 1;
            if (viewType === 'top') {
                depthFactor = 0.5 + 0.5 * (atom.z - this.zMin) / (this.zMax - this.zMin);
            } else {
                depthFactor = 0.5 + 0.5 * (atom.y - this.yMin) / (this.yMax - this.yMin);
            }
            
            const radius = atom.radius * scale * depthFactor * 3; // Scale up for visibility
            
            // Draw atom
            this.ctx.beginPath();
            this.ctx.arc(x, y, radius, 0, Math.PI * 2);
            
            // Get atom color
            this.ctx.fillStyle = this.getAtomColor(atom.type);
            
            // Add gradient for 3D effect
            const gradient = this.ctx.createRadialGradient(
                x - radius * 0.3, y - radius * 0.3, radius * 0.1,
                x, y, radius
            );
            gradient.addColorStop(0, '#ffffff');
            gradient.addColorStop(1, this.getAtomColor(atom.type));
            this.ctx.fillStyle = gradient;
            
            this.ctx.fill();
            
            // Add outline
            this.ctx.strokeStyle = 'rgba(0,0,0,0.2)';
            this.ctx.lineWidth = 1;
            this.ctx.stroke();
        }
        
        // Add view label
        this.ctx.fillStyle = 'rgba(0,0,0,0.5)';
        this.ctx.font = '14px Arial';
        this.ctx.textAlign = 'left';
        this.ctx.fillText(`${viewType.toUpperCase()} VIEW`, 10, 20);
    }
}