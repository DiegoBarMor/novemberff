from setuptools import setup, find_packages

setup(
    name="novemberff",
    version="0.3.0",
    description="Custom implementation of the AMBER forcefield, focused on simplicity and verbose outputs.",
    keywords="Amber Forcefield MolecularDynamics MD Simulation",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="DiegoBarMor",
    author_email="diegobarmor42@gmail.com",
    url="https://github.com/diegobarmor/novemberff",
    license="MIT",
    packages=find_packages(),
    install_requires=[],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
)
