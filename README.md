## Presc_pre_post_covid
Repo to freshly analyse prescription data

Code: 
The Code/ directory streamlines all the work done till date on extraction of the geographical trends in prescriptions from the NHS data. 
The scripts sequentially help in matching conditions/drug categories to a set of drugs, 
which are then matched to a set of BNF codes, which in-turn are used to map prescriptions to LSOAs in England. 
The code is supposed to be run on a pre-serialized set of NHS prescriptions files (which you can download if you don't want to serialize again from here: https://www.dropbox.com/sh/09367w8uf4hbla5/AAAJIHQPq4Lw_4R--AyY7xIua?dl=0)
If you do wish to serialize form scratch, download the raw prescriptions form the NHS website and run the serialize.py on all of them. 

## Notebooks: 
The notebooks are there to do some exploratory analysis of the prescriptions. 

## Mappings: 
This directory contains all the mappings required for generating spatio-temporal trends from raw prescroptions. 
