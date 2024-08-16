from setuptools import setup, find_packages

setup(
    name='EXSCOPE',
    version='0.1.0',
    description='A tool for visualizing read counts at specified genomic regions.',
    author='Taimoor Khan',
    author_email='taimoorkhan@scichores.com',
    url='https://github.com/taimoorkhan-bt/exscope',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'exscope = exscope.exscope:main',
        ],
    },
    install_requires=[
        'pysam',
        'matplotlib',
        'pandas',
        'numpy',
        'scipy'
    ],
    python_requires='>=3.6',
)
