The gfweakivtest.ado file extends the routine weakivtest.ado of

Pflueger, C.E. and S. Wang (2015), "A Robust Test for Weak Instruments in Stata,” The Stata Journal, 216–225,

which produces the effective F-statistic and Nagar bias related critical values developed in 

Montiel Olea, J. L. and C. Pflueger (2013), “A Robust Test for Weak Instruments,” Journal of Business & Economic Statistics, 31, 358–369.

The routine gfweakivtest.ado incorporates further results for the robust F-statistic as a test for weak instruments related to the Nagar bias of the GMMf estimator, see

Windmeijer F. (2025), "The Robust F-Statistic as a Test for Weak Instruments,” Journal of Econometrics.

It outputs estimation results for 2SLS, LIML and GMMf, the non-robust F-statistic, the Effective F-statistic and the robust F-statistic as well as the critical values 
for the effective and robust F-statistics.

The gfweakivtestols.ado file outputs the critical values related to the worst-case OLS bias, instead of the benchmark bias of Montiel-Olea and Pflueger (2013).


For example:

. use http://fmwww.bc.edu/ec-p/data/wooldridge/mroz.dta

. qui ivreg2 lwage exper expersq (educ=age kidslt6 kidsge6), rob
 
. gfweakivtest

Montiel-Pflueger robust weak instrument test

--------------------------------------------
btsls:                       0.0964

sebtsls:                     0.0865

bliml:                       0.0958

sebliml:                     0.0913

kappa:                       1.0016

Non-Robust F statistic:       4.342

Effective F statistic:        4.552

Confidence level alpha:          5%

--------------------------------------------

--------------------------------------------
Critical Values             TSLS      LIML

--------------------------------------------
% of Worst Case Bias

tau=5%                    15.711    15.406

tau=10%                    9.957     9.789

tau=20%                    6.749     6.654

tau=30%                    5.560     5.491

--------------------------------------------
--------------------------------------------
bgmmf:                       0.0948

sebgmmf:                     0.0868

Robust F statistic:           5.021

Confidence level alpha:          5%

--------------------------------------------

--------------------------------------------
Critical Values             GMMf

--------------------------------------------
% of Worst Case Bias

tau=5%                    13.651

tau=10%                    8.745

tau=20%                    6.021

tau=30%                    5.018

--------------------------------------------

. gfweakivtestols

Montiel-Pflueger robust weak instrument test

--------------------------------------------
btsls:                       0.0964

sebtsls:                     0.0865

bliml:                       0.0958

sebliml:                     0.0913

kappa:                       1.0016

Non-Robust F statistic:       4.342

Effective F statistic:        4.552

Confidence level alpha:          5%

--------------------------------------------

--------------------------------------------
Critical Values             TSLS      LIML

--------------------------------------------
% of Worst Case Bias

tau=5%                    15.900    15.406

tau=10%                   10.062     9.789

tau=20%                    6.808     6.654

tau=30%                    5.602     5.491

--------------------------------------------
--------------------------------------------
bgmmf:                       0.0948

sebgmmf:                     0.0868

Robust F statistic:           5.021

Confidence level alpha:          5%

--------------------------------------------

--------------------------------------------
Critical Values             GMMf

--------------------------------------------
% of Worst Case Bias

tau=5%                    13.901

tau=10%                    8.882

tau=20%                    6.098

tau=30%                    5.073

--------------------------------------------

