# MACA
MACA: Marker-based automatic cell-type annotation for single cell expression data

## 1. Installment
     
     pip install scanpy ##install prerequisite 
     
     pip install MACA-Python

## 2. Tutorial for basic use of MACA
See /Tutorial/Basic use of MACA/

     |--------|----------------|---MACA_tutorial.ipynb
     
## 3. Tutorial for integrated annotation
See /Tutorial/Integrated annotation via MACA/

     |--------|------------------------------|---MACA_integrated_annotation_humanheart.ipynb

     |--------|------------------------------|---MACA_integrated_annotation_humanpancreas.ipynb
    
     |--------|------------------------------|---MACA_integrated_annotation_humanPBMC.ipynb
     
## 4. Standardization of cell type annotation across COVID19 datasets via MACA

See /Tutorial/Integrated annotation via MACA/

     |--------|------------------------------|---MACA_integrated_annotation_COVID19.ipynb

![alt text](https://github.com/ImXman/MACA/blob/master/Tutorial/Integrated%20annotation%20via%20MACA/Figure%201.jpg?raw=true)
    

# Update 03/12/2021

     MACA was modified for parallel computing. 
     
     For combined ~647K single nuclei human heart data (Tucker et al., Circulation 2020 and Litviňuková et al., Nature 2020), annotation through MACA takes 24 mins with NMI as 0.739 and ARI as 0.818 against authors' annotations. 
     
     See /Tutorial/Integrated annotation via MACA/MACA_integrated_annotation_humanheart.ipynb
