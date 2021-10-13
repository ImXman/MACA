# MACA
MACA: Marker-based automatic cell-type annotation for single cell expression data

While our manuscript is under review, if you are using our tool for your research, we kindly ask you to cite MACA github. Thanks!

## 1. Installment
MACA works anndata format and is compatible with pipeline analysis through scanpy

     pip install scanpy, anndata, scikit-learn, cosg ##install prerequisite Scanpy 1.6.0 and Anndata 0.7.5
     
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

## 5. Label transferring via MACA (Ongoing)

See /Tutorial/Label transferring via MACA/

     MACA_transfer_annotation.ipynb

![alt text](https://github.com/ImXman/MACA/blob/master/Tutorial/Label%20transferring%20via%20MACA/Figure%202.jpg?raw=true)

The 1st row shows MACA's annotations, and the 2nd row shows author's annotations. In the 3rd and 4th rows, cells are colored by study source and used as either reference or target data.
    

# Update 03/12/2021

MACA was modified for parallel computing. For combined ~647K single nuclei human heart data (Tucker et al, Circulation 2020 and Litviňuková et al, Nature 2020), annotation through MACA takes 24 mins with NMI as 0.739 and ARI as 0.818 against authors' annotations.

# Update 10/07/2021

Marker identification function was added to identify robust differential markers across cell-types, and it is used for label transferring based on reference data
     
# Statement

GPU-supported research has speeded up integrative discoveries across single-cell studies. However, access to a good graphic card for model training is not taken granted, especially in undeveloped and developing countries. Even renting a gpu instance on the cloud is costy for researchers.

We devote to make integrative single-cell analysis accessible for most people, and MACA is a cheap solution to label transferring for large single-cell data. MACA annotates 1 million cells for 40 minutes, on a personal laptop with i7-8550U CPU, 16GB memory, and no GPU support.

 
