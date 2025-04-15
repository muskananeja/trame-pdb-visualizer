@echo off
REM This script runs the Molecular Visualizer on Windows systems

REM Check if Python is installed
python --version >nul 2>&1
if %errorlevel% neq 0 (
    echo Python is not installed or not in PATH. Please install Python and try again.
    exit /b
)

REM Create a virtual environment if it doesn't exist
if not exist venv (
    echo Creating virtual environment...
    python -m venv venv
)

REM Activate the virtual environment
echo Activating virtual environment...
call venv\Scripts\activate

REM Install dependencies
echo Installing dependencies...
pip install -e .

REM Create necessary directories
if not exist static\examples mkdir static\examples

REM Download example PDB file if it doesn't exist
if not exist static\examples\1cbs.pdb (
    echo Downloading example PDB file...
    powershell -Command "Invoke-WebRequest -Uri https://files.rcsb.org/download/1CBS.pdb -OutFile static\examples\1cbs.pdb"
)

REM Run the application
echo Starting the Molecular Visualizer...
echo Open your browser and go to: http://localhost:5000
python app.py

REM Deactivate the virtual environment when done
call venv\Scripts\deactivate