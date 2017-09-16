<!-- README.md is generated from README.Rmd. Please edit that file -->
PSTR
====

The PSTR package implements the Panel Smooth Transition Regression (PSTR) modelling.

The modelling procedure consists of three stages: Specification, Estimation and Evaluation. The package offers tools helping the users to conduct model specification tests, to do PSTR model estimation, and to do model evaluation.

The cluster-dependency and heteroskedasticity-consistent tests are implemented in the package.

The wild bootstrap and cluster wild bootstrap tests are also implemented.

Parallel computation (as an option) is implemented in some functions, especially the bootstrap tests. Therefore, the package suits tasks running many cores on super-computation servers.

How to install
--------------

You can either install the stable version from CRAN

``` r
install.packages("PSTR")
```

or install the develop version from GitHub

``` r
devtools::install_github{"yukai-yang/PSTR"}
```

provided that the package "devtools" has been installed beforehand.

Example
-------

After installing the package, you need to load (attach better say) it by running the code

``` r
library(PSTR)
```

You can first check the information and the current version number by running

``` r
version()
#> #########################################################################
#> ## package name: PSTR
#> ## author: Yukai Yang
#> ## Department of Statistics
#> ## Uppsala University
#> ## yukai.yang@statistik.uu.se
#> ## Version 1.0.2 Sep. 2017
#> #########################################################################
```

Then you can take a look at all the available functions and data in the package

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
               tvars=c('vala'), im=1, iT=14)
print(pstr)
#> #########################################################################
#> ## package name: PSTR
#> ## Version 1.0.2 Sep. 2017
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
#>   vala
#> ***********************************************************************
#> #########################################################################
```

It says that the data set "Hansen99" is used, the dependent variable is "inva", the variables in the data from column 4 to 20 are the explanatory variables in the linear part (though you can write down the names of them), the explanatory variables in the nonlinear part are the four ones in "indep\_k", and the potential transition variable is "vala" (Tobin's Q).

Now you can see that the "NewPSTR" is basically defining the settings of the model.

Note that you can print the object of the class PSTR. By default, it gives you a summary of the PSTR model. They are mainly about which one is the dependent variable, which ones are explanatory variables and etc..

The following code does linearity tests

``` r
pstr = LinTest(use=pstr) 
print(pstr, "tests")
#> #########################################################################
#> ## package name: PSTR
#> ## Version 1.0.2 Sep. 2017
#> #########################################################################
#> ***********************************************************************
#> Results of the linearity (homogeneity) tests:
#> -----------------------------------------------------------------------
#> LM tests based on transition variable 'vala'
#>   m  LM_X PV  LM_F PV HAC_X        PV HAC_F        PV
#>   1 125.3  0 28.99  0 30.03 4.819e-06 6.952 1.396e-05
#> ***********************************************************************
#> Sequence of homogeneity tests for selecting number of switches 'm':
#> -----------------------------------------------------------------------
#> LM tests based on transition variable 'vala'
#>   m  LM_X PV  LM_F PV HAC_X        PV HAC_F        PV
#>   1 125.3  0 28.99  0 30.03 4.819e-06 6.952 1.396e-05
#> ***********************************************************************
#> #########################################################################
```

You can see that the function "LinTest" takes the PSTR object "pstr" and overwrites it when return. This is the way I recommend as the functions handling the PSTR object in the package update the object by adding new atrributes or members. However, the same function will change the values of the attributes it adds. You can of course create new PSTR objects to take the return values in order to save the results from different settings of the model.

You can do the wild bootstrap and wild cluster bootstrap by running the following code. (Warning! Don't run it except that you have at least 50 cores!)

``` r
iB = 5000 # the number of repetitions in the bootstrap
library(snowfall)
pstr = WCB_LinTest(use=pstr,iB=iB,parallel=T,cpus=50)
```

It takes a long long time to run the bootstrap. This function is developed for those who work on some super-computation server with many cores and a large memory. Note that you will have to attach the "snowfall" package manually.

But of course, you can try the function on your personal computer by reducing the number of repetitions and the cores.

``` r
pstr = WCB_LinTest(use=pstr,iB=4,parallel=T,cpus=2)
```

When you determine which transition variable to use for the estimation, in this case "inva", you can estimate the PSTR model

``` r
pstr = EstPSTR(use=pstr,im=1,iq=1,par=c(1.6,.5), vLower=4, vUpper=4)
print(pstr,"estimates")
```

By default, the "optim" method "L-BFGS-B" is used, but you can change the method for estimation by doing

``` r
pstr = EstPSTR(use=pstr,im=1,iq=1,par=c(1.6,.5), method="CG")
print(pstr,"estimates")
#> #########################################################################
#> ## package name: PSTR
#> ## Version 1.0.2 Sep. 2017
#> #########################################################################
#> ***********************************************************************
#> Results of the PSTR estimation:
#> -----------------------------------------------------------------------
#> Transition variable 'vala' is used in the estimation.
#> -----------------------------------------------------------------------
#> Parameter estimates in the linear part (first extreme regime) are
#>        dt_75_0   dt_76_0   dt_77_0  dt_78_0  dt_79_0  dt_80_0  dt_81_0
#> Est  -0.004332 -0.007436 -0.005040 0.001092 0.003012 0.006406 0.001119
#> s.e.  0.002502  0.002586  0.002738 0.002799 0.002697 0.002862 0.002944
#>        dt_82_0   dt_83_0  dt_84_0  dt_85_0  dt_86_0   dt_87_0   vala_0
#> Est  -0.007943 -0.014010 0.000836 0.004903 0.000628 -0.006537 -0.01722
#> s.e.  0.002677  0.002757 0.003085 0.003234 0.003104  0.003192  0.06771
#>       debta_0   cfa_0   sales_0
#> Est  -0.06389 0.08598 -0.000464
#> s.e.  0.02035 0.03117  0.004060
#> -----------------------------------------------------------------------
#> Parameter estimates in the non-linear part are
#>       vala_1 debta_1    cfa_1 sales_1
#> Est  0.02403 0.06768 -0.04671 0.01033
#> s.e. 0.06772 0.02021  0.03373 0.00494
#> -----------------------------------------------------------------------
#> Parameter estimates in the second extreme regime are
#>      vala_{0+1} debta_{0+1} cfa_{0+1} sales_{0+1}
#> Est    0.006810    0.003784   0.03927    0.009865
#> s.e.   0.001193    0.011720   0.01235    0.003116
#> -----------------------------------------------------------------------
#> Non-linear parameter estimates are
#>      gamma    c_1
#> Est  4.953 0.4949
#> s.e. 1.211 0.2538
#> -----------------------------------------------------------------------
#> Estimated standard deviation of the residuals is 0.04323
#> ***********************************************************************
#> #########################################################################
```

Thus, the evaluation tests can be done based on the estimated model

``` r
## evaluatio tests
pstr1 = EvalTest(use=pstr,vq=pstr$mQ[,1])
```

Note that in the "EvalTest", only one transition variable is taken each time for the no remaining nonlinearity test. This is different from the "LinTest" function which can take several transition variables. This is the reason why I save the results into new PSTR objects "pstr1" instead of overwriting. By doing so, I can save more test results from different transition variables in new objects.

The user can also do the wild bootstrap and wild cluster bootstrap in the following way, provided that he or she has the super-computation resources.

``` r
iB = 5000
cpus = 50

## wild bootstrap time-varyint evaluation test 
pstr = WCB_TVTest(use=pstr,iB=iB,parallel=T,cpus=cpus)

## wild bootstrap heterogeneity evaluation test
pstr1 = WCB_HETest(use=pstr1,vq=pstr$mQ[,1],iB=iB,parallel=T,cpus=cpus)
```
