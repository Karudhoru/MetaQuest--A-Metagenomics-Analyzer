from setuptools import setup, find_packages

setup(
    name='MetaQuest',
    version='0.1.2',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'plotly',
        'seaborn',
        'Bio',
        'matplotlib',
        'scikit-learn',
        'scipy',
        'numpy',
        'xgboost',
        'lightgbm',
        'catboost'
    ],
    entry_points={
        'console_scripts': [
            'metagenomics-cli = app:main'
        ]
    }
)
