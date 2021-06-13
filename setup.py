import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="MACA-Python",
    version="1.0.1",
    author="Yang Xu",
    author_email="yxu71@vols.utk.edu",
    description="MACA: Marker-based automatic cell-type annotation for single cell expression data",
    #long_description=long_description,
    #long_description_content_type="text/markdown",
    url="https://github.com/ImXman/MACA",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
