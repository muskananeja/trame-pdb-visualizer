import os
import io
from flask import Flask, send_from_directory, request, jsonify, render_template, send_file

# Initialize Flask app
app = Flask(__name__)

def parse_pdb_info(pdb_content):
    """Parse basic information from PDB file"""
    lines = pdb_content.splitlines()
    
    # Count atoms, residues, and chains
    atoms = 0
    residues = set()
    chains = set()
    
    for line in lines:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atoms += 1
            try:
                residue_id = line[22:27].strip()
                chain_id = line[21]
                residues.add(residue_id + chain_id)
                chains.add(chain_id)
            except:
                pass
    
    return {
        "atoms": atoms,
        "residues": len(residues),
        "chains": len(chains)
    }

@app.route('/')
def index():
    # Serve the index.html file directly as static content
    return send_file('templates/index.html')

@app.route('/load_example')
def load_example():
    """Load example PDB file"""
    example_path = os.path.join(os.path.dirname(__file__), "static/examples/1cbs.pdb")
    
    if os.path.exists(example_path):
        with open(example_path, 'r') as f:
            pdb_content = f.read()
        
        # Parse PDB info
        info = parse_pdb_info(pdb_content)
        info["filename"] = "1cbs.pdb"
        
        print(f"Loaded PDB with {info['atoms']} atoms, {info['residues']} residues, {info['chains']} chains")
        
        return jsonify({
            "info": info,
            "content": pdb_content
        })
    else:
        return jsonify({"error": f"Example file not found: {example_path}"}), 404

@app.route('/upload_pdb', methods=['POST'])
def upload_pdb():
    """Handle PDB file upload"""
    # Check if request has a file
    if 'file' not in request.files:
        return jsonify({"error": "No file part in the request"}), 400
    
    file = request.files['file']
    
    # If user doesn't select a file
    if file.filename == '':
        return jsonify({"error": "No file selected"}), 400
    
    # Check if file is a PDB
    if not file.filename.lower().endswith('.pdb'):
        return jsonify({"error": "Not a PDB file (must end with .pdb)"}), 400
    
    # Read the file content
    pdb_content = file.read().decode('utf-8')
    
    # Validate the file has atom entries
    if not 'ATOM' in pdb_content and not 'HETATM' in pdb_content:
        return jsonify({"error": "Invalid PDB file format (no ATOM or HETATM entries)"}), 400
    
    # Parse PDB info
    info = parse_pdb_info(pdb_content)
    info["filename"] = file.filename
    
    print(f"Uploaded PDB with {info['atoms']} atoms, {info['residues']} residues, {info['chains']} chains")
    
    return jsonify({
        "info": info,
        "content": pdb_content
    })

@app.route('/static/<path:path>')
def send_static(path):
    return send_from_directory('static', path)

if __name__ == '__main__':
    # Create static directory if it doesn't exist
    os.makedirs('static/examples', exist_ok=True)
    
    # Check if example file exists, download if not
    example_path = 'static/examples/1cbs.pdb'
    if not os.path.exists(example_path):
        import urllib.request
        url = 'https://files.rcsb.org/download/1CBS.pdb'
        try:
            urllib.request.urlretrieve(url, example_path)
            print(f"Downloaded example PDB file to {example_path}")
        except Exception as e:
            print(f"Error downloading example file: {e}")
    
    # Run the Flask app
    app.run(host='0.0.0.0', port=5000, debug=True)
