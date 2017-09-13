<!-- README.md is generated from README.Rmd. Please edit that file -->
PSTR
====

The PSTR package implements the Panel Smooth Transition Regression (PSTR) modelling.

The modelling procedure consists of three stages: Specification, Estimation and Evaluation. The package offers tools helping the users to conduct model specification tests, to do PSTR model estimation, and to do model evaluation.

The cluster-dependency and heteroskedasticity-consistent tests are implemented in the package.

The wild bootstrap and cluster wild bootstrap tests are also implemented.

Parallel computation (as an option) is implemented in some functions, especially the bootstrap tests. Therefore, the package suits tasks running many cores on super-computation servers.

Example
-------

After installing and attaching the package, you can take a look at all the available functions and data in the package

``` r
ls( grep("PSTR", search()) ) 
#> [1] "EstPSTR"     "EvalTest"    "Hansen99"    "LinTest"     "NewPSTR"    
#> [6] "version"     "WCB_HETest"  "WCB_LinTest" "WCB_TVTest"
```

In the package, a data set called "Hansen99" is offered to give prompt example. For details of the data set, you can run

``` r
?Hansen99 
```

You can create a new object of the class PSTR by doing

``` r
pstr = NewPSTR(Hansen99, dep='inva', indep=4:20, indep_k=c('vala','debta','cfa','sales'),
               tvars=c('vala','debta','cfa','sales'), iT=14)
print(pstr)
#> #########################################################################
#> ## package name: PSTR
#> ## Version 1.0.0, Sep. 2017
#> #########################################################################
#> ***********************************************************************
#> Summary of the model:
#> -----------------------------------------------------------------------
#>   time horizon sample size = 14,  number of individuals = 560
#> -----------------------------------------------------------------------
#> Dependent variable:  inva
#> -----------------------------------------------------------------------
#> Explanatory variables in the linear part:
#>   dt_75 dt_76 dt_77 dt_78 dt_79 dt_80 dt_81 dt_82 dt_83 dt_84 dt_85 dt_86 dt_87 vala debta cfa sales
#> -----------------------------------------------------------------------
#> Explanatory variables in the non-linear part:
#>   vala debta cfa sales
#> -----------------------------------------------------------------------
#> Potential transition variable(s) to be tested:
#>   vala debta cfa sales
#> ***********************************************************************
#> #########################################################################
```

It says that the data set "Hansen99" is used, the dependent variable is "inva", the variables in the data from column 4 to 20 are the explanatory variable in the linear part (though you can write down the names of them instead), the explanatory variables in the nonlinear part are the four ones in "indep\_k", and the potential transition variables are those four.

Note that you can print the object of the class PSTR. By default, it gives you a summary of the PSTR model. They are mainly about which one is the dependent variable, which ones are explanatory variables and etc..

The following code does linearity tests

``` r
pstr = LinTest(use=pstr,im=3) 
print(pstr, "tests")
#> #########################################################################
#> ## package name: PSTR
#> ## Version 1.0.0, Sep. 2017
#> #########################################################################
#> ***********************************************************************
#> Results of the linearity (homogeneity) tests:
#> -----------------------------------------------------------------------
#> LM tests based on transition variable 'vala'
#>   m  LM_X PV  LM_F PV HAC_X        PV HAC_F        PV
#>   1 125.3  0 28.99  0 30.03 4.819e-06 6.952 1.396e-05
#>   2 217.4  0 25.15  0 55.01 4.404e-09 6.363 2.939e-08
#>   3 290.8  0 22.42  0 76.52 1.895e-11 5.897 2.610e-10
#> -----------------------------------------------------------------------
#> LM tests based on transition variable 'debta'
#>   m  LM_X        PV   LM_F        PV HAC_X       PV HAC_F       PV
#>   1 37.81 1.227e-07  8.752 4.809e-07 13.72 0.008238 3.176 0.012860
#>   2 86.87 1.998e-15 10.050 4.929e-14 21.87 0.005159 2.530 0.009535
#>   3 89.71 5.618e-14  6.914 1.341e-12 24.22 0.018990 1.867 0.033470
#> -----------------------------------------------------------------------
#> LM tests based on transition variable 'cfa'
#>   m  LM_X PV  LM_F PV HAC_X        PV HAC_F        PV
#>   1 128.9  0 29.83  0 20.41 4.140e-04 4.725 8.307e-04
#>   2 142.3  0 16.46  0 35.77 1.935e-05 4.138 6.059e-05
#>   3 206.0  0 15.88  0 47.18 4.347e-06 3.636 1.836e-05
#> -----------------------------------------------------------------------
#> LM tests based on transition variable 'sales'
#>   m   LM_X PV  LM_F PV HAC_X        PV HAC_F        PV
#>   1  94.83  0 21.95  0 15.10 0.0045070 3.494 0.0074070
#>   2 116.50  0 13.47  0 29.91 0.0002193 3.460 0.0005485
#>   3 136.30  0 10.50  0 31.59 0.0016000 2.435 0.0037060
#> ***********************************************************************
#> Sequence of homogeneity tests for selecting number of switches 'm':
#> -----------------------------------------------------------------------
#> LM tests based on transition variable 'vala'
#>   m   LM_X        PV  LM_F        PV HAC_X        PV HAC_F        PV
#>   1 125.30 0.000e+00 28.99 0.000e+00 30.03 4.819e-06 6.952 1.396e-05
#>   2  93.68 0.000e+00 21.67 0.000e+00 22.12 1.896e-04 5.118 4.079e-04
#>   3  75.50 1.554e-15 17.46 2.887e-14 24.62 6.010e-05 5.692 1.431e-04
#> -----------------------------------------------------------------------
#> LM tests based on transition variable 'debta'
#>   m   LM_X        PV    LM_F        PV  HAC_X       PV  HAC_F      PV
#>   1 37.810 1.227e-07  8.7520 4.809e-07 13.720 0.008238 3.1760 0.01286
#>   2 49.300 5.046e-10 11.4100 3.147e-09 10.960 0.027010 2.5360 0.03818
#>   3  2.867 5.804e-01  0.6628 6.178e-01  1.316 0.858700 0.3043 0.87530
#> -----------------------------------------------------------------------
#> LM tests based on transition variable 'cfa'
#>   m   LM_X        PV   LM_F        PV  HAC_X       PV  HAC_F        PV
#>   1 128.90 0.000e+00 29.830 0.000e+00 20.410 0.000414 4.7250 0.0008307
#>   2  13.67 8.418e-03  3.163 1.316e-02  1.894 0.755200 0.4382 0.7811000
#>   3  64.83 2.791e-13 14.990 3.316e-12  6.682 0.153700 1.5450 0.1862000
#> -----------------------------------------------------------------------
#> LM tests based on transition variable 'sales'
#>   m  LM_X        PV   LM_F        PV  HAC_X       PV  HAC_F       PV
#>   1 94.83 0.0000000 21.950 0.0000000 15.100 0.004507 3.4940 0.007407
#>   2 21.89 0.0002104  5.065 0.0004488  5.262 0.261500 1.2170 0.301100
#>   3 20.11 0.0004747  4.650 0.0009505  2.907 0.573500 0.6722 0.611200
#> ***********************************************************************
#> #########################################################################
```

You can see that the function Lintest takes the PSTR object "pstr" and overwrites it when return. This is the way I recommend as the functions handling the PSTR object in the package update the object by adding new atrributes or members. However, the same function will change the values of the attributes which it updates.

In order to avoid the loss of some important information in the object due to overwriting, consider the following example

``` r
pstr1 = LinTest(use=pstr,im=1) 
pstr2 = LinTest(use=pstr,im=2) 
```

You can of course create new PSTR objects to take the return values.

You can do the wild bootstrap and wild cluster bootstrap by running the following code. (Warning! Don't run it except that you have at least 50 cores!)

``` r
iB = 5000 # the number of repetitions in the bootstrap
pstr = WCB_LinTest(use=pstr,im=3,iB=iB,parallel=T,cpus=50)
```

It takes a long long time to run the bootstrap. This function is developed for those who work on some super-computation server with many cores and a large memory.

When you determine which transition variable to use for the estimation, in this case, the first one in tvars "inva", you can estimate the PSTR model

``` r
pstr = EstPSTR(use=pstr,im=1,iq=1,par=c(1,1), vLower=4, vUpper=4)
print(pstr,"estimates")
#> #########################################################################
#> ## package name: PSTR
#> ## Version 1.0.0, Sep. 2017
#> #########################################################################
#> ***********************************************************************
#> Results of the PSTR estimation:
#> -----------------------------------------------------------------------
#> Transition variable 'vala' is used in the estimation.
#> -----------------------------------------------------------------------
#> Parameter estimates in the linear part (first extreme regime) are
#>        dt_75_0   dt_76_0   dt_77_0   dt_78_0 dt_79_0  dt_80_0  dt_81_0
#> Est  -0.003951 -0.007785 -0.005575 0.0005957 0.00280 0.006312 0.000949
#> s.e.  0.002403  0.002552  0.002656 0.0027880 0.00271 0.002920 0.002934
#>        dt_82_0   dt_83_0  dt_84_0  dt_85_0   dt_86_0   dt_87_0   vala_0
#> Est  -0.008131 -0.014400 0.000482 0.004950 0.0009812 -0.005942 0.039070
#> s.e.  0.002589  0.002668 0.003047 0.003216 0.0030900  0.003099 0.006902
#>       debta_0   cfa_0    sales_0
#> Est  -0.03852 0.06654 -0.0009666
#> s.e.  0.01496 0.01506  0.0044910
#> -----------------------------------------------------------------------
#> Parameter estimates in the non-linear part are
#>         vala_1 debta_1    cfa_1  sales_1
#> Est  -0.031630 0.04593 -0.02962 0.014680
#> s.e.  0.006834 0.01748  0.02047 0.004748
#> -----------------------------------------------------------------------
#> Parameter estimates in the second extreme regime are
#>      vala_{0+1} debta_{0+1} cfa_{0+1} sales_{0+1}
#> Est    0.007434    0.007403   0.03692    0.013710
#> s.e.   0.001380    0.011630   0.01406    0.003183
#> -----------------------------------------------------------------------
#> Non-linear parameter estimates are
#>      gamma   c_1
#> Est  2.718 1.000
#> s.e. 0.723 0.358
#> -----------------------------------------------------------------------
#> Estimated standard deviation of the residuals is 0.04328
#> ***********************************************************************
#> #########################################################################
```

Thus, the evaluation tests can be done based on the estimated model

``` r
## evaluatio tests
pstr1 = EvalTest(use=pstr,im=1,vq=pstr$mQ[,1])
pstr2 = EvalTest(use=pstr,im=1,vq=pstr$mQ[,2])
pstr3 = EvalTest(use=pstr,im=1,vq=pstr$mQ[,3])
pstr4 = EvalTest(use=pstr,im=1,vq=pstr$mQ[,4])
```

Note that in the "EvalTest", only one transition variable is taken each time for the no remaining nonlinearity test. This is different from the "LinTest" function which can take several transition variables. This is the reason why I save the results into new PSTR objects instead of overwriting.

The user can also do the wild bootstrap and wild cluster bootstrap in the following way, provided that he or she has the super-computation resources.

``` r
iB = 5000
cpus = 50

## wild bootstrap time-varyint evaluation test 
pstr = WCB_TVTest(use=pstr,im=1,iB=iB,parallel=T,cpus=cpus)

## wild bootstrap heterogeneity evaluation test
pstr1 = WCB_HETest(use=pstr1,im=1,vq=pstr$mQ[,1],iB=iB,parallel=T,cpus=cpus)
pstr2 = WCB_HETest(use=pstr2,im=1,vq=pstr$mQ[,2],iB=iB,parallel=T,cpus=cpus)
pstr3 = WCB_HETest(use=pstr3,im=1,vq=pstr$mQ[,3],iB=iB,parallel=T,cpus=cpus)
pstr4 = WCB_HETest(use=pstr4,im=1,vq=pstr$mQ[,4],iB=iB,parallel=T,cpus=cpus)
```
