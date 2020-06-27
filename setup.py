import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="topo_lab_tools-gaberber", # Replace with your own username
    version="0.0.1",
    author="Guanzhong Wang",
    author_email="cnhajzwgz@gmail.com",
    description="Simple codes for data analysis in the topo lab",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/gaberber/topo_lab_tools",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)