from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name='MetaQuest',
    version='3.5.0',
    author='MetaQuest Development Team',
    author_email='devpatelcambay@gmail.com',
    description='Comprehensive metagenomics analysis pipeline with ML-based pathogen detection',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/yourusername/MetaQuest',
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    python_requires='>=3.8',
    install_requires=requirements,
    entry_points={
        'console_scripts': [
            'metaquest=metaquest.cli:main',
        ],
    },
    include_package_data=True,
    package_data={
        'metaquest.ml': ['model_artifacts/*'],
    },
)