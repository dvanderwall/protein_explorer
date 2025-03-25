# Use Python 3.9 slim as base image
FROM python:3.9-slim

# Set working directory
WORKDIR /app

RUN echo "Starting Docker build for Protein Explorer"

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copy requirements first to leverage Docker cache
RUN echo "Copying requirements.txt"
COPY requirements.txt .

# Install pandas explicitly first to ensure it's available
RUN pip install --no-cache-dir pandas numpy scipy scikit-learn networkx gunicorn flask flask-cors feather-format pyarrow requests tqdm matplotlib jinja2 seaborn flask-caching biopython
RUN echo "Installing Python dependencies from requirements.txt"
RUN pip install --no-cache-dir -r requirements.txt

# Copy application code
RUN echo "Copying application code"
COPY web_app/ ./web_app/
COPY protein_explorer/ ./protein_explorer/
COPY setup.py .
COPY README.md .
RUN echo "Application code copied successfully"

# Create cache directory
RUN echo "Creating cache directory"
RUN mkdir -p /root/.protein_explorer/cache
RUN echo "Cache directory created at /root/.protein_explorer/cache"

# Copy data files
RUN echo "Copying data files"
# Note: These files should exist in your build context
COPY Combined_Kinome_10A_Master_Filtered_2.feather .
RUN echo "Copied Combined_Kinome_10A_Master_Filtered_2.feather"
COPY Sequence_Similarity_Edges.parquet .
RUN echo "Copied Sequence_Similarity_Edges.parquet"
COPY Structure_Kinase_Scores.feather .
RUN echo "Copied Structure_Kinase_Scores.feather"
COPY Sequence_Kinase_Scores.feather .
RUN echo "Copied Sequence_Kinase_Scores.feather"
COPY protein_explorer/PhosphositeSuppData.feather .
RUN echo "Copied PhosphositeSuppData.feather"
RUN echo "All data files copied successfully"

# Set environment variables
RUN echo "Setting environment variables"
ENV FLASK_APP=web_app/app.py
ENV FLASK_ENV=production
ENV PYTHONPATH=/app
RUN echo "Environment setup complete"

# Expose port
EXPOSE 8080
RUN echo "Exposing port 8080"

# Final directory listing for verification
RUN echo "Final directory contents:" && ls -la
RUN echo "Python path:" && python -c "import sys; print(sys.path)"
RUN echo "Checking for pandas:" && python -c "import pandas; print(f'Pandas version: {pandas.__version__}')"

# Command to run the application
CMD ["gunicorn", "--bind", "0.0.0.0:8080", "web_app.app:app"]