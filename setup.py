from setuptools import setup
import versioneer

requirements = [
    "pandas>=1.3.3",
    "tqdm==4.50.0",
    "jsonapi_client",
    "requests",
    "urllib3",
]

setup(
    setup_requires=[
        # Setuptools 18.0 properly handles Cython extensions.
        "setuptools>=18.0",
        "Cython>=0.29.24",
        "pandas>=1.3.3",
        "tqdm==4.50.0",
        "jsonapi_client",
        "requests",
        "urllib3",
    ],
    name="get-biomes",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="A simple tool to get biomes related data from MGnify and ENA",
    license="GNUv3",
    author="Antonio Fernandez-Guerra",
    author_email="antonio@metagenomics.eu",
    url="https://github.com/genomewalker/get-biomes",
    packages=["get_biomes"],
    entry_points={"console_scripts": ["getBiomes=get_biomes.__main__:main"]},
    install_requires=requirements,
    keywords="get-biomes",
    classifiers=[
        "Programming Language :: Python :: 3.9",
    ],
)
