from setuptools import setup, find_packages

setup(
    name="molecular-visualizer",
    version="1.0.0",
    description="3D Interactive Molecular Structure Viewer",
    author="Muskan Aneja",
    author_email="muskan.aneja@example.com",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "flask>=2.0.1",
        "numpy>=1.21.0",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Web Environment",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    python_requires=">=3.6",
)