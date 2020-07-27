import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="KACA",
    version="0.0.1",
    author="Yang Xu",
    author_email="yxu71@vols.utk.edu",
    description="KACA: Knowledge-based automatic cell-type annotation for single cell expression data",
    #long_description=long_description,
    #long_description_content_type="text/markdown",
    url="https://github.com/ImXman/KACA",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU GPLv3",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)