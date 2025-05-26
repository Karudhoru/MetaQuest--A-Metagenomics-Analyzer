from setuptools import setup, find_packages

setup(
    name='MetaQuest',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'plotly',
        'seaborn',
        'biopython',
        'matplotlib'
    ],
    entry_points={
        'console_scripts': [
            'metagenomics-cli = app:main'
        ]
    }
)
