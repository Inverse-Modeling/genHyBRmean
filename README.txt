  ***************************************************************************
  * All the software  contained in this library  is protected by copyright. *
  * Permission  to use, copy, modify, and  distribute this software for any *
  * purpose without fee is hereby granted, provided that this entire notice *
  * is included  in all copies  of any software which is or includes a copy *
  * or modification  of this software  and in all copies  of the supporting *
  * documentation for such software.                                        *
  ***************************************************************************
  * THIS SOFTWARE IS BEING PROVIDED "AS IS", WITHOUT ANY EXPRESS OR IMPLIED *
  * WARRANTY. IN NO EVENT, NEITHER  THE AUTHORS, NOR THE PUBLISHER, NOR ANY *
  * MEMBER  OF THE EDITORIAL BOARD OF  THE JOURNAL  "NUMERICAL ALGORITHMS", *
  * NOR ITS EDITOR-IN-CHIEF, BE  LIABLE FOR ANY ERROR  IN THE SOFTWARE, ANY *
  * MISUSE  OF IT  OR ANY DAMAGE ARISING OUT OF ITS USE. THE ENTIRE RISK OF *
  * USING THE SOFTWARE LIES WITH THE PARTY DOING SO.                        *
  ***************************************************************************
  * ANY USE  OF THE SOFTWARE  CONSTITUTES  ACCEPTANCE  OF THE TERMS  OF THE *
  * ABOVE STATEMENT AND OF THE ACCOMPANYING FILE LICENSE.txt.               *
  ***************************************************************************

   AUTHORS:
        
        Taewon Cho, Department of Mathematics, 
        Virginia Tech

        Julianne Chung, Department of Mathematics, 
        Virginia Tech
        
        Scot M. Miller, Department of Environmental Health and Engineering, 
        Johns Hopkins University
        
        Arvind K. Saibaba, Department of Mathematics, 
        North Carolina State University
   
   REFERENCE:

       "Hybrid Projection Method for Large-scale Inverse Problems with Mean 
        Estimation in Hierarchical Gaussian Priors". 
            2021.

   SOFTWARE LANGUAGE:

       MATLAB 9.6 (R2019a)


=====================================================================
SOFTWARE and DATA SETS
=====================================================================
The MainDrivers require the following package:

    "genHyBR"
    by Julianne Chung and Arvind Saibaba
    https://github.com/juliannechung/genHyBR

and require the data sets:
    
    "Geostatistical inverse modeling with large atmospheric data: 
    data files for a case study from OCO-2"
    by Scot M. Miller, Arvind K. Saibaba, Michael E. Trudeau, 
    Marikate E. Mountain, and Arlyn E. Andrews
    https://doi.org/10.5281/zenodo.3241466

- Include genHyBR package in the current path.
- Put all 'H_i' matrices and 'distmat.mat' in the folders "6wk/H/" and "1yr/H/"
- Put all 's_i' vectors and corresponding Z_5, Z_10, Z_50, and Z_100
  in the folder "6wk/s/" and "1yr/s/"

MainDrivers for each numerical experiments of inverse problems

  Driver_Hybrid.m          Run genHyBRs or genHyBRmean to 6wk & 1yr case study 
                           with Tikhonov regularization for various noise levels.
                           Solve either for only unknowns(CO2 fluxes)
                           or for both unknowns(CO2 fluxes) and mean 
                           coefficients simultaneously.
  
  Driver_Direct.m          Run direct method to 6wk case study.
                           Solve for unknowns(CO2 fluxes).
                           1yr case is not recommended with direct method.
