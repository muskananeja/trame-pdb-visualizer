#!/bin/bash

# This script runs the Molecular Visualizer on Mac/Linux systems

# Check if Python is installed
if ! command -v python3 &> /dev/null; then
    echo "Python 3 is not installed. Please install Python 3 and try again."
    exit 1
fi

# Check if pip is installed
if ! command -v pip3 &> /dev/null; then
    echo "pip is not installed. Please install pip and try again."
    exit 1
fi

# Create a virtual environment if it doesn't exist
if [ ! -d "venv" ]; then
    echo "Creating virtual environment..."
    python3 -m venv venv
fi

# Activate the virtual environment
echo "Activating virtual environment..."
source venv/bin/activate

# Install dependencies
echo "Installing dependencies..."
pip install -e .

# Create necessary directories
mkdir -p static/examples

# Download example PDB file if it doesn't exist
if [ ! -f "static/examples/1cbs.pdb" ]; then
    echo "Downloading example PDB file..."
    curl -o static/examples/1cbs.pdb https://files.rcsb.org/download/1CBS.pdb
fi

# Run the application
echo "Starting the Molecular Visualizer..."
echo "Open your browser and go to: http://localhost:5000"
python app.py

# Deactivate the virtual environment when done
deactivate