FROM python:3.9-slim

WORKDIR /app

# Copy requirements and setup files
COPY setup.py .
COPY app.py .

# Copy static files and templates
COPY static/ ./static/
COPY templates/ ./templates/

# Create examples directory and install dependencies
RUN mkdir -p static/examples && \
    pip install --no-cache-dir -e .

# Download example PDB file
RUN apt-get update && \
    apt-get install -y curl && \
    curl -o static/examples/1cbs.pdb https://files.rcsb.org/download/1CBS.pdb && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Expose the port the app runs on
EXPOSE 5000

# Run the application
CMD ["python", "app.py"]