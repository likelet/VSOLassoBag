# Simulation Experiments for LassoBag Framework

This is a README file for **Simulation Experiments for LassoBag Framework**. For your further understanding of our Simulation results, we document the __detailed data format__ we used in simulation experiments as well as the way to amend your own version of simulation data according to different choice of distribution. More information can be referred to . 

We store the simulation data as __RDS__ file during the simulation. However, we will also support using __cvs__ or __txt__ format in the future. 


## Simulation Data Format 
1. __X__ (Independent variables, stored in `X.rds`):

   - `$X_Matrix`: 
   
      A matrix with `ncol` and `nrow` that can be defined by input. Simplely change `p` (for number of column, aka **Predictors**) and `n` (for number of rows, aka **Sample numbers**)
   - `$X_info`: 
   
      A list each of which is corresponded to a dimension in __X_Matrix__ with sub-attribute `$type` as the distribution **type** that generated the column data and `$para` as corresponding **parameters** used.
   

2. **coY_level** (Dependent Variable and corresponding coeffients, stored in `coY_level_n.rds`):

    - `$matrix_coeffs`:
    
       A **matrix** that contain the coefficients that used to generate the __Y values__. It was generated in a __cross-validation__ fashion in order to test the performance of the framework in handling high dimensional data with **multipul overlap** across dimensions.  
       
       **Note**: there are some basic level coefficient vectors from which multipul level of combination was used to generate the rest of the coefficient vectors. The way they combine is documented as the column name in the matrix.
     
    - `$matrix_Y`:
    
       A **matrix** that comtain the __Y value__ (*Dependent variables*) that generated from **X_Matrix** and **cofficient vectors**. Each column represents a Y vector that is used for learning and training. Therefore the number of columns of this matrix should be equal to the number of **Predictors** and the number of rows should be the same as the **Sample numbers**.   
       
       **Note**: The column name also contained the information with which you can trace back to the exact coefficient vector that generates it.
       
## Try yourself
