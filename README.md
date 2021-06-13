# MACA
MACA: Marker-based automatic cell-type annotation for single cell expression data

## Installment
pip install MACA

## Tutorial for basic use of MACA
See /Tutorial/Basic use of MACA/MACA_tutorial.ipynb

## Tutorial for integrated annotation
See /Tutorial/Integrated annotation via MACA/

    |--------|------------------------------|---MACA_integrated_annotation_humanheart.ipynb

    |--------|------------------------------|---MACA_integrated_annotation_humanpancreas.ipynb
    
    |--------|------------------------------|---MACA_integrated_annotation_humanPBMC.ipynb
    

# Update 03/12/2021

MACA was modified for parallel computing.

For combined ~647K single nuclei human heart data 

(Tucker et al., Circulation 2020 and Litviňuková et al., Nature 2020), 

annotation through MACA takes 24 mins with NMI as 0.739 and ARI as 0.818 against authors' annotations. 

See /Tutorial/Integrated annotation via MACA/MACA_integrated_annotation_humanheart.ipynb
