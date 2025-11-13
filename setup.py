from setuptools import find_packages, setup

setup(
    name="aphylogeo",
    url="https://github.com/tahiri-lab/aPhyloGeo",
    author="Tahiri Lab",
    description="A phylogenetic and geographic analysis tool",
    version="2.0.0",
    include_package_data=True,
    install_requires=[
        "bio==1.5.9",
        "multiprocess==0.70.15",
        "pandas==2.1.1",
        "psutil==5.9.5",
        "PyYAML==6.0.1",
        "pytest==7.4.2",
        "numpy==1.26.0",
        "cython==3.0.2",
        "robinson-foulds==1.2",
        "ete3==3.1.3",
        "dendropy==4.6.1",
        "textdistance==4.6.0",
        "typer==0.13.1"
    ],
    python_requires=">=3.9.0",
    packages=find_packages(),
)
