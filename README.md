README FILE: Bayesian Variable Selection Under High-dimensional Settings with Grouped Covariates
Authors: Pranay Agarwal, Subhajit Dutta, Minerva Mukhopadhyay

1.	Abstract: This folder contains the codes corresponding to the following methods proposed in the paper “Bayesian Variable Selection Under High-dimensional Settings with Grouped Covariates”:
   
(a)	M-GiVSA: pMTM algorithm applied to the Modified-g prior setup
(b)	G-GiVSA: pMTM algorithm applied to the Group informed-g prior setup
(c)	GSISd: The GSIS group inclusion probabilities for a chosen d

The codes corresponding to (a) and (b) are provided in the GiVSA_codes sub-folder, and the codes corresponding to (c) are provided in the GSIS_codes subfolder. 

3.	Application of M-GiVSA and G-GiVSA:
   
2A.	To apply M-GiVSA or G-GiVSA methods, one needs to execute the following steps:

(a)	Set the working directory to the GiVSA_codes folder.

(b)	Install the following R-packages:
a.	corpcor
b.	Rcpp
c.	lme4
d.	parallel
e.	Matrix

(c)	Save the data in a RData file named input.RData, which contains the following objects:
a.	 X: The design matrix in nXp format, where n is the no. of observations, and p is the number of covariates.
b.	 y: The n vector of responses
c.	true_group_structure: A p-length vector indicating the groups indices as numbers. 
     The input.RData file should be placed in GiVSA_codes folder.

(d)	Execute the steps described in Application.R file.

2B.	The GiVSA_codes folder contains the following five R files,

i.	GiVSA.R: This file contains the Givsa function, which executes the GiVSA algorithm applied to the Modified g-prior as well as the Group-informed g-prior setup. The inputs and outputs of the function Givsa are given in GiVSA.R.

ii.	Output_GiVSA.R: This file contains the function Output_GiVSA, which evaluates the results obtained from the function Givsa. The inputs and outputs of the function Output_GiVSA are described in Output_GiVSA.R.

iii.	Application.R: This file contains the codes for applying M-GiVSA or G-GiVSA methods on a given data set. 

iv.	Readme_GiVSA_codes.R: This file contains a description of all the codes provided in the GiVSA_codes folder. 
                 the following C++ file,

v.	Hash.CPP: This is a helper function of Givsa and is required for faster implementation of the codes.
       and the following RData file.

vi.	input.RData: This file is required for the Application.R file and it contains an example data set on which the codes are applied.

2C.	The Readme_GiVSA_codes.R contains a description of all the functions and files related to M-GiVSA and G-GiVSA methods and is present in the GiVSA_codes folder. Further, the execution of M-GiVSA and G-GiVSA with an example data set (contained in input.RData) is shown in the Application.R file.

3.	Application of GSISd: 

3A.	To apply GSISd methods, one needs to execute the following steps:

(a)	Set the working directory to the GSIS_codes folder.

(b)	Install the following R-package:	corpcor

(c)	Save the data in a RData file named input.RData, which contains the following objects:

a.	 X: The design matrix in nXp format, where n is the no. of observations, and p is the number of covariates.

b.	 y: The n vector of responses

c.	true_group_structure: A p-length vector indicating the groups indices as numbers. 
     The input.RData file should be placed in GSIS_codes folder.

(d)	Execute the steps described in Application.R file.

3B.	The GSIS_codes folder contains the following five R files,

i.	GSIS.R: This file contains the GSISd function, which calculates the GSIS group inclusion probabilities for a given data, group structure and choice of d. The inputs and outputs of the function GSIS are given in GSIS.R.

ii.	Application.R: This file contains the codes for applying GSIS methods on a given data set. 

iii.	Readme_GSIS_codes.R: This file contains a description of all the codes provided in the GSIS_codes folder. 
       and the following RData file.

iv.	input.RData: This file is required for the Application.R file and it contains an example data set on which the codes are applied.

3C.	The Readme_GSIS.R contains a description of all the functions and files related to GSIS method and is present in the GSIS_codes folder. Further, the execution of GSIS with an example data set (contained in input.RData) is shown in the Application.R file.
