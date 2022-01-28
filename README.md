# MACA
MACA: Marker-based automatic cell-type annotation for single cell expression data

MACA was published on Bioinformatics but a free preprint version can be found at https://www.biorxiv.org/content/10.1101/2021.10.25.465734v1.

## 1. Installment
MACA works anndata format and is compatible with pipeline analysis through scanpy

     pip install scanpy, anndata, scikit-learn ##install prerequisite Scanpy 1.6.0 and Anndata 0.7.5
     
     pip install MACA-Python ##it's the only version of MACA. We update installment once we update new MACA version

## 2. Tutorial for basic use of MACA
See /Tutorial/Basic use of MACA/

     MACA_tutorial.ipynb
     
## 3. Tutorial for integrated annotation
See /Tutorial/Integrated annotation via MACA/

     MACA_integrated_annotation_humanheart.ipynb

     MACA_integrated_annotation_humanpancreas.ipynb
    
     MACA_integrated_annotation_humanPBMC.ipynb
     
## 4. Standardization of cell type annotation across COVID19 datasets via MACA

See /Tutorial/Integrated annotation via MACA/

     MACA_integrated_annotation_COVID19.ipynb

![alt text](https://github.com/ImXman/MACA/blob/master/Tutorial/Integrated%20annotation%20via%20MACA/Figure%201.jpg?raw=true)

## 5. Cell-type annotation for 10X Visium data

See /Tutorial/

     MACA_transfer_annotation_spatialbrain10xVisium.ipynb

![alt text](https://github.com/ImXman/MACA/blob/master/Tutorial/Figure2.jpg?raw=true)

# Update 03/12/2021

MACA was modified for parallel computing. For combined ~647K single nuclei human heart data (Tucker et al, Circulation 2020 and Litviňuková et al, Nature 2020), annotation through MACA takes 24 mins with NMI as 0.739 and ARI as 0.818 against authors' annotations.

# Update 11/14/2021

We established a new github repo named MASI (https://github.com/ImXman/MASI), which combines reference data and MACA for fast label transferring.
     
# Statement

GPU-supported research has speeded up integrative discoveries across single-cell studies. However, access to a good graphic card for model training is not taken granted, especially in undeveloped and developing countries. Even renting a gpu instance on the cloud is costy for researchers.

We devote to make integrative single-cell analysis accessible for most people, and MACA is a cheap solution to label transferring for large single-cell data. MACA annotates 1 million cells for 40 minutes, on a personal laptop with i7-8550U CPU, 16GB memory, and no GPU support.
