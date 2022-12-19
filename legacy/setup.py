import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gasification-NicolasAnese",
    version="1.0",
    author="Rodolfo Rodrigues",
    author_email="rodolfo.rodrigues@ufsm.br",
    description="Gasification equilibrium model",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/NicolasAnese/gasification",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    python_requires=">=3.6",
)