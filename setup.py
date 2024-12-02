from setuptools import setup, find_packages
import pathlib
here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.rst").read_text(encoding="utf-8")

def read_requirements(fname):
    with open(fname, 'r', encoding='utf-8') as file:
        return [line.rstrip() for line in file]

setup(
    name='Augusta',
    version='1.0.6',
    packages=find_packages(),
    url='https://github.com/JanaMus/Augusta',
    license='MIT',
    author='Jana Musilova, Zdenek Vafek, Karel Sedlar',
    author_email='musilovajana@vut.cz',
    description='Python package for inference of the gene regulatory network and the boolean network using RNA-Seq data.',
    long_description=long_description,
    keywords = ['Computational biology', 'Bioinformatics', 'RNA-Seq', 'mutual information', 'database', 'Boolean network', 'Gene Regulatory network', 'SBML'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 4 - Beta"],
    python_requires='>=3.7, <3.9',
    install_requires=read_requirements('requirements.txt'),
)
