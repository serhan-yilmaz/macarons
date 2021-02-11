# Macarons
#### Macarons: a fast and simple algorithm to select complementary SNPs
## Requirements
Running Macarons requires some packages such as `numpy` and `scipy`, which are documented in the `requirements.txt` file. To install them, run
```
pip install -r requirements.txt
```
Note that the current version of Macarons has been tested on Python3.8, but should work with any Python3 version.

## Running Macarons on Python
A demo is given in `demo.py`. To run it, simply use
```
python demo.py
```
The demo uses the 4W phenotype data from the Arabidopsis Thaliana (AT) obtained from [Atwell et. al. (2010)](https://www.ncbi.nlm.nih.gov/pubmed/20336072). 


## References
Yilmaz, S., Fakhouri, M., Koyuturk, M., Cicek, A. E. and Tastan, O. (2020). [Uncovering complementary sets of variants for predicting quantitative phenotypes](https://doi.org/10.1101/2020.12.11.419952). bioRxiv

Atwell, S. et al. (2010) [Genome-wide association study of 107 phenotypes
in Arabidopsis thaliana inbred lines](https://www.ncbi.nlm.nih.gov/pubmed/20336072). Nature, 465(7298), 627–631.

