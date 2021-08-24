import setuptools
from setuptools.command.install import install
from io import open
import os

install_requires = ["numpy", "scipy", "pysam", "pandas", "biopython", "bcbio-gff"]
setuptools.setup(
    name='CMSeq',
    version='1.0.4',
    author='Moreno Zolfo',
    author_email='moreno.zolfo@unitn.it',
    url='http://github.com/SegataLab/cmseq/',
    license = 'LICENSE.txt',
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': [
            'breadth_depth.py = cmseq.breadth_depth:bd_from_file',
            'consensus.py = cmseq.consensus:consensus_from_file',
            'consensus_aDNA.py = cmseq.consensus_aDNA:consensus_from_file',
            'polymut.py = cmseq.polymut:polymut_from_file',
            'poly.py = cmseq.poly:poly_from_file'
        ]
    },
    long_description_content_type='text/markdown',
    long_description=open('README.md').read(),
    description='Set of utilities on sequences and BAM files',
    install_requires=install_requires
)
