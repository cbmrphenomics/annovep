from setuptools import find_packages, setup

setup(
    name="annovep",
    version="2",
    description="Annotation pipeline for VCF files using VEP",
    url="https://github.com/cbmrphenomics/annovep",
    author="Mikkel Schubert",
    author_email="Mikkel.Schubert@sund.ku.dk",
    license="MIT",
    packages=find_packages(),
    install_requires=[
        "coloredlogs>=7.3",
        "isal>=0.11.1",
        "liftover>=1.1.16",
        "pydantic>=2.4.2",
        "pysam>=0.16.0.1",
        "ruamel.yaml>=0.18.5",
        "typing_extensions>=4.7.1",
    ],
    entry_points={"console_scripts": ["annovep=annovep.__main__:main"]},
    include_package_data=True,
)
