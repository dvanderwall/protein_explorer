from setuptools import setup, find_packages

setup(
    name="protein-explorer",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "requests>=2.25.0",
        "numpy>=1.20.0",
        "scipy>=1.6.0",
        "networkx>=2.5.0",
        "plotly>=4.14.0",
        "scikit-learn>=0.24.0",
        "flask>=2.0.0",
        "biopython>=1.78",
    ],
    author="Your Name",
    author_email="your.email@example.com",
    description="A package for exploring protein structures and interactions",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/class-account/protein-explorer",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)