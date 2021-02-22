# PRER: A Patient Representation with Pairwise Relative Expression of Proteins on Biological Networks

---
## Abstract
Alterations in protein and gene expression levels are often used as features to predictive models such as clinical outcome prediction. A common strategy to combine signals on individual proteins is to integrate alterations with biological knowledge. In this work, we propose a novel patient representation where we integrate the expression levels of proteins with the biological networks. Patient representation with PRER (Pairwise Relative Expressions with Random walks) operates in the neighborhood of a protein and aims to capture the dysregulation patterns in protein abundance for proteins that are known to interact. This neighborhood of the source protein is derived using a biased random-walk strategy on the network. Specifically, PRER computes a feature vector for a patient by comparing the protein expression level of the source protein with other proteins’ levels in its neighborhood. We test PRER’s performance through a survival prediction task in 10 different cancers using random forest survival models. PRER representation yields a statistically significant predictive performance in 8 out of 10 cancer types when compared to a representation based on individual protein expression. We also identify the set of proteins that are important not because of alteration of its expression values but due to the alteration in their pairwise relative expression values. The set of identified relations provides a valuable collection of biomarkers with high prognostic value. PRER representation can be used for other complex diseases and prediction tasks that use molecular expression profiles as input.

---

## Authors
Halil Ibrahim Kuru, Mustafa Buyukozkan, Oznur Tastan

---

## Instructions Manual

### Requirements
Required R packages
- randomForestSRC
- survival 
- glmnet
- foreach
- doParallel

### Random walk data
Download preprocessed random walk data from <a href="https://drive.google.com/file/d/1KRjSuWVj9INrauBOSMHY8vupNM-QmF92/view?usp=sharing">**link**</a>, save it to `randomWalk/`

You can also generate your own random walk data by using the following command.
```shell
python randomWalk/main.py --input randomWalk/ppi_edgelist.csv --output randomWalk/randomWalks.txt --walk-length 100 --num-walks 18 --p 0.25 --q 0.25 --workers 8
```
### Reading and processing raw data
```shell
Run main_reader.R file
```

### Training and testing PRER model
Use the processed data that we read in main_reader.R file. Then:
```shell
Run main_rsf.R file
```

---

## License


- **[CC BY-NC-SA 2.0](https://creativecommons.org/licenses/by-nc-sa/2.0/)**
- **[MIT License](https://github.com/hikuru/matchmaker/blob/master/LICENSE)**
