
EvalCertification
-------------

Copyright (C) 2021 [Parker
Edwards](https://sites.nd.edu/parker-edwards)

External requirements
---------------------

1. [Maple](https://www.maplesoft.com/)


Description
-----------

Maple library for evaluating functions at roots of polynomials. 
See the article http://arxiv.org/abs/2102.00115 for
theoretical details and Section 5 of that article for discussion
of the included examples. It has been tested on Maple 2020.


The included workbooks for detailed examples highlight
different input options. The main library is *EvalCertificaiton.mpl*.

Version 1.0.0
-------------

Basic usage for *EvalCertification.mpl*
---------------------------
To read the library into your project, assuming the file is in your
current Maple directory.
``` perl
    read("EvalCertification.mpl")
```
The main export is the procedure `EstimateRootsAndCertifyEvaluations`
```perl
    S := EstimateRootsAndCertifyEvaluations(
    polynomial_to_solve,
    [function_to_evaluate_1,function_to_evaluate_2,...,function_to_evaluate_n],
    HolderFunction,
    EstimationPrecision)
```

-  **polynomial_to_solve** is a polynomial with rational (Maple data type fraction) coefficients.
-  **[function_to_evaluate_1,...,function_to_evaluate_n]** is a list of       locally Hölder continuous functions.
-  **HolderFunction** is a procedure which computes Hölder information for functions
   in the evaluation list. For example, the functions are polynomials, then    **HolderFunction** could be 
   the built in procedure **HolderInformationForPolynomial**  provided by **EvalCertification**.
-  **EstimationPrecision** is a positive rational number.

The object **S** is a Record with fields

-  **S:-root-values** a list of the estimated roots of **polynomial_to_solve**.
-  **S:-evaluations_functions_i** for i=1,2,..., n is a list of certified 
   evaluations of **function_to_evaluate_i** at the roots of **polynomial_to_solve**. 
   The evaluations are in the same order as the roots in **S:-root_values**.

License
-------

EvalCertification is licensed under an MIT license.
