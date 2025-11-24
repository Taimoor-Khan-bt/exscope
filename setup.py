from setuptools import setup, find_packages

setup(
    name="exscope",
    version="1.2.1",
    description="A genomics visualization and analysis tool integrating established bioinformatics tools.",
    author="Taimoor Khan",
    author_email="taimoorkhan@scichores.com",
    packages=find_packages(),
    install_requires=[
        # Core dependencies - use flexible versions for Python 3.12 compatibility
        "pysam>=0.22.0",
        "matplotlib>=3.5.0",
        "pandas>=1.5.0",
        "numpy>=1.24.0",
        "gffutils>=0.12",
        "tqdm>=4.65.0",
        
        # Visualization dependencies
        "plotly>=5.14.0",
        "kaleido>=0.2.0",
        "pygenometracks>=3.7",
        "pybedtools>=0.9.0"
    ],
    extras_require={
        'dev': [
            'pytest>=7.0.0',
            'pytest-cov>=4.0.0',
            'black>=23.0.0',
            'isort>=5.12.0',
            'flake8>=6.0.0',
        ],
    },
    entry_points={
        'console_scripts': [
            'exscope=exscope.cli:main',
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Operating System :: OS Independent",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.8',
    project_urls={
        'Source': 'https://github.com/Taimoor-Khan-bt/exscope',
        'Bug Reports': 'https://github.com/Taimoor-Khan-bt/exscope/issues',
    }
)
