from setuptools import setup, find_packages

setup(
    name="exscope",
    version="0.2.0",
    description="A tool to extract read counts from BAM files and visualize them.",
    author="Taimoor Khan",
    author_email="taimoorkhan@scichores.com",
    packages=find_packages(),
    install_requires=[
        "pysam",
        "matplotlib",
        "pandas",
        "numpy",
        "gffutils",
        "tqdm",
        "scikit-learn",
        "seaborn"
    ],
    entry_points={
        'console_scripts': [
            'exscope=exscope.main:main',
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
