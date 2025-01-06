//	Authors: 
//	    Carolin E. Pflueger, University of British Columbia, Vancouver BC V6T 1Z2, Canada
//      carolin.pflueger@sauder.ubc.ca
//      Su Wang, London School of Economics, London WC2A 2AE, United Kingdom
//      s.wang50@lse.ac.uk.

//  The weakivtest ado file creates a STATA command, -weakivtest-, implementing
//	the test for weak instruments of Montiel Olea, J. L. and C. E. Pflueger.(2013) 
//	"A robust test for weak instruments". Journal of Business and Economic Statistics 31:358-369.
//	For more details of weakivtest routine, please refer to Pflueger, C. E. and Su Wang(2013)
//	"A robust test for weak instruments in Stata." http://papers.ssrn.com/sol3/papers.cfm?abstract_id=2323012.
* this version: 10/28/2020 - include noconstant option
 
program gfweakivtestols,rclass
	version 10

	*Allowed options: significant level, optimization parameter eps and parameter for the inverse non-central F distribution n2.
	syntax [,level(string) eps(string) n2(string)]
	
	*If "avar" command is not installed in user's computer, display the error message.
	capture which avar
	if _rc==111 {
		di as err `"User contributed command avar is needed to run weakivtest. Install avar by typing "ssc install avar"."'
		exit
	}
	
	*weakivtest is a postestimation command to ivreg2 or ivregress; otherwise display error message.
	if `"`e(cmd)'"' != "ivreg2" & `"`e(cmd)'"' != "ivregress" {
		di as err `"Weakivtest is a postestimation command after running ivreg2 or ivregress."'
		exit
	}
	
	*CP 10/28/2020: allow for noconstant
	/*
	*weakivtest is a postestimation command to ivreg2 or ivregress with constant; otherwise display error message.
	if ( "`e(constant)'"=="noconstant") | ( "`e(cons)'"=="0") {
		di as err `"Noconstant option is not supported by weakivtest."'
		exit
	}
	*/
	*Su 10/30/2020: check for noconst_flag
	local noconst_flag=( "`e(constant)'"=="noconstant") | ( "`e(cons)'"=="0")
	*weakivtest only allows one single endogenous regressor; otherwise display error message.
	if  ( wordcount("`e(instd)'")!=1) {
		di as err `"Weakivtest requires one endogenous regressor."'
		exit
	}
	
	*Save inputs and set default values in syntax.
	if "`level'"==""     local level = 0.05
	if "`eps'"==""       local eps=   0.001
	if "`n2'"==""		 local n2=1e7

/* Save results from the ereturn values of ivreg2 or ivregress.*/
	
	tempname K L S myereturn
		
	*Hold and copy ereturn list
    _estimates hold  `myereturn' , copy
	
	*Save ereturn results after ivreg2
	if `"`e(cmd)'"' == "ivreg2"{
	
		*Number of excluded exogenous regressors.
		mat `K'=`e(exexog_ct)'
		local exexog_num=`K'[1,1]
		
		*Number of included exogenous regressors.	
		mat `L'=`e(inexog_ct)'
		
		*Number of all observations in e(sample). 
		mat `S'=`e(N)'
	
		*Save vce type to be used in "avar" command to compute variance-covariance matrices. 
		*** PH FIX 8/17/18: adjustment for clustering d.f.
		local clustdfadj = 1
		if strpos("`e(vce)'","bw"){
			local vcet=regexr("`e(vce)'","ac bartlett", "")
			local vcet=regexr("`vcet'","h", "")
			local vcet=regexr("`vcet'","=", "(")+")"
		}
		else if strpos("`e(vce)'","cluster"){
			local vcet "robust cluster(`e(clustvar)')"
			local clustdfadj = (`e(N)'/(`e(N)'-1))*(`e(N_clust)'-1)/(`e(N_clust)')
		}
		else {
			local vcet `e(vce)'
		}
	}	
	
	*Save ereturn results after ivregress
	if `"`e(cmd)'"' == "ivregress"{
		
		*Number of included exogenous regressors adjusted for unused variables listed in `e(exogr)'.
		local inexog_num: word count `e(exogr)'
		foreach var in `e(exogr)'{
			if strpos("`var'","b.")|strpos("`var'","o."){
				local inexog_num=`inexog_num'-1
			}
		}
		mat `L'=`inexog_num'
		
		*Number of excluded exogenous regressors adjusted for unused variables listed in `e(insts)'.
		local exexog_num: word count `e(insts)'
		foreach var in `e(insts)'{
			if strpos("`var'","b.")|strpos("`var'","o."){
				local exexog_num=`exexog_num'-1
			}
		}
		local exexog_num=`exexog_num'-`inexog_num'
		mat `K'=`exexog_num'
		
		*Number of all observations in e(sample). 
		mat `S'=`e(N)'
		
	
		*Save vce type to be used in "avar" command to compute variance-covariance matrices. 
		*** PH FIX 8/17/18: adjustment for clustering d.f.
		local clustdfadj = 1
		if strpos("`e(vce)'","unadjusted"){ 
			local vcet ""
		}
		else if strpos("`e(vce)'","robust"){ 
			local vcet "robust "
		}
		else {
			local vcet ""
			if strpos("`e(vce)'","hac bartlett"){
				local lagnum = substr("`e(cmdline)'",strpos("`e(cmdline)'","hac nw")+length("hac nw"),2)
				local lagnum = `lagnum'+1
				local vcet= "`vcet'"+"robust bw"+"("+"`lagnum'"+")"
			}	
			else if strpos("`e(cmdline)'","cluster"){
				local vcet= "`vcet'"+ "cluster (" +"`e(clustvar)'" +")"
				local clustdfadj = (`e(N)'/(`e(N)'-1))*(`e(N_clust)'-1)/(`e(N_clust)')
			}
			else {
				local vcet `e(vce)'
			}
		}
	}

/* Prepare data to compute variance covariance matrix.*/

	* Save the data in temp variables for ivreg2.
	if `"`e(cmd)'"' == "ivreg2"{
		fvrevar `e(exexog)'
		local eexexog  "`r(varlist)'" 
		fvrevar `e(instd)'
		local einstd   "`r(varlist)'"  
		fvrevar `e(depvar)'
		local edepvar  "`r(varlist)'" 
		fvrevar `e(inexog)'

		local einexog "`r(varlist)'"  
		local depvartemp ""
		foreach var in `edepvar' {
			tempvar `var'temp
			qui clonevar  ``var'temp'=`var'
			local depvartemp ``var'temp' `depvartemp'
		}
		
		local endotemp ""
		foreach var in `einstd' {
			tempvar `var'temp
			qui clonevar  ``var'temp'=`var'
			local endotemp ``var'temp' `endotemp'
		}
		
		local exexogtemp ""
		foreach var in `eexexog' {
			tempvar `var'temp
			qui clonevar  ``var'temp'=`var'
			local exexogtemp ``var'temp' `exexogtemp'
		}
		
		local Zincltemp ""
		foreach var in `einexog' {
			tempvar `var'temp
			qui clonevar  ``var'temp'=`var'
			local Zincltemp ``var'temp' `Zincltemp'
		}
	}

	* Save the data in temp variables for ivregress.
	if `"`e(cmd)'"' == "ivregress"{
	
		* Delete the unused variables saved in the ereturn results after ivregress.
		local einsts ""
		foreach var in `e(insts)'{
			if strpos("`var'","b.")==0 & strpos("`var'","o.")==0 {
				local einsts `einsts' `var'
			}
		}
		
		local eexogr ""
		foreach var in `e(exogr)'{
			if strpos("`var'","b.")==0 & strpos("`var'","o.")==0 {
				local eexogr `eexogr' `var'
			}
		}
		if "`eexogr'"!=""{
			local instsname "`einsts'"
			local eexexog: subinstr local instsname  `"`eexogr'"' ""
		}
		else {
			local eexexog "`einsts'"
		}
		
		* Generate temporary variables and save the tempvar names accordingly.
		local depvartemp ""
		foreach var in `e(depvar)' {
			tempvar `var'temp
			qui clonevar  ``var'temp'=`var'
			local depvartemp ``var'temp' `depvartemp'
		}
	
		local endotemp ""
		foreach var in `e(instd)' {
			tempvar `var'temp
			qui clonevar  ``var'temp'=`var'
			local endotemp ``var'temp' `endotemp'
		}
	
		local exexogtemp ""
		foreach var in `eexexog' {
			local tpnm=strtoname("`var'",1)
			tempvar tpv tpnm
			qui gen `tpv'=`var'
			rename `tpv' `tpnm'
			local exexogtemp `exexogtemp' `tpnm'
		}
		
		local Zincltemp ""
		foreach var in `eexogr' {
			local tpnm=strtoname("`var'",1)
			tempvar tpv tpnm
			qui gen `tpv'=`var'
			rename `tpv' `tpnm'
			local Zincltemp `Zincltemp' `tpnm'
		}
	}
	
	* Use observations from e(sample).
	tempvar touse
	gen byte `touse' = e(sample)
   
	* Demean endotemp, exexogtemp and depvartemp variables by regressing them on included exogenous regressors if there is any. 
	*** PH FIX 8/17/18: adjustment for weights
	if `"`e(wtype)'"' != "" {
		local addweight "[`e(wtype)'`e(wexp)']"
        local accuweight = "`addweight'"
        if `"`e(wtype)'"' == "pweight" {
            local accuweight "[aw`e(wexp)']"
        }
	}

	*** PH FIX 8/17/18: adjustment for weights
	foreach var in `depvartemp' `endotemp'{
		tempvar rsd 
		
		*CP 10/28/2020: allow for noconstant
		if ( `noconst_flag'==1) {
		qui reg `var' `Zincltemp' if `touse' `addweight', noconstant
		}
		else {
		qui reg `var' `Zincltemp' if `touse' `addweight'	
		}

		qui predict double `rsd' if `touse' , r
		drop `var'
		rename `rsd' `var'
	}
	
	*** PH FIX 8/17/18: adjustment for weights
	foreach var in `exexogtemp' {
		 tempvar exrsd 
		 
		 *CP 10/28/2020: allow for noconstant
		 if (`noconst_flag'==1) {
		 qui reg `var' `Zincltemp' if `touse' `addweight', noconstant
		 }
		 else{
		 qui reg `var' `Zincltemp' if `touse' `addweight'
		 }
		 qui predict double `exrsd'  if `touse' , r
		 drop `var'
		 rename `exrsd' `var'
	}

	* Orthogonalize the instruments.
	tempname Zs 
	*** PH FIX 8/17/18: adjustment for weights
	qui orthog `exexogtemp' if `touse' `addweight', gen(`Zs'*)

	* Save matricies in mata.
	mata `K'=st_matrix("`K'")
	mata `S'=st_matrix("`S'")
	mata `L'=st_matrix("`L'")

/* Get full W matrix and Omega matrix using residuals from two stage regressions.*/

	tempname W W_1 W_2 W_12 Omega omega_22 res1 res2 F F_eff F_rob W_2mh trW_1s W_12s trW_12s ///
             eigs mineig maxeig pihat dhat btsls sebtsls  W_2inv bgmmf sebgmmf ///
             restsls S_1tsls resgmmf S_1gmmf S_1BGMMf KK ww kappa xhatl bliml ///
             Ommh limlden limlnum limlres Wl sebliml limlres2
	
        
	* Run first stage regression and save residuals as v1 in mata.
	*** PH FIX 8/17/18: adjustment for weights
	qui reg `endotemp' `Zs'* `addweight', noconstant
    mat `pihat' = e(b)'
	qui predict double `res1', r
    mata `pihat' = st_matrix("`pihat'")
    
	* Run second stage regression and save residuals as v2 in mata.
	*** PH FIX 8/17/18: adjustment for weights
	qui reg `depvartemp' `Zs'* `addweight', noconstant
    mat `dhat' = e(b)'
	qui predict double `res2', r
	mata `dhat' = st_matrix("`dhat'") 

        
	* Omega matrix is the unadjusted variance covariance matrix adjusted by the degree of freedom.
	*** PH FIX 8/17/18: adjustment for weights
	qui mat accum `Omega'= `res2' `res1' `accuweight', nocons
	
	mata `Omega'=st_matrix("`Omega'")/`S'
	
    qui mat accum `ww'= `depvartemp' `endotemp' `accuweight', nocons
    mata `ww' = st_matrix("`ww'")/`S'
    mata `Ommh' = matpowersym(`Omega',-0.5)
    mata `kappa' = symeigenvalues(`Ommh'*`ww'*`Ommh')
    mata `kappa' = rowmin(`kappa')
    mata st_numscalar("kappa",`kappa')
    mata st_local("kappa",strofreal(`kappa'))
    
    qui gen `xhatl' = `endotemp'-`kappa'*`res1'
    qui mat accum `limlden'= `xhatl' `endotemp' `depvartemp' `accuweight', nocons
    mat `limlnum' = `limlden'[1,3]
    mat `limlden' = `limlden'[1,2]
    mata `limlnum' = st_matrix("`limlnum'")/`S'
    mata `limlden' = st_matrix("`limlden'")/`S'
    mata `bliml' = `limlnum'/`limlden'
    mata st_numscalar("bliml",`bliml')
    mata st_local("bliml",strofreal(`bliml'))
    
    qui gen `limlres' = `depvartemp'-`endotemp'*bliml
    
*    local vcet ""
*    di strpos("`vcet'","")

    
    if (strpos("`vcet'","")==0) {
        qui mat accum `limlres2' = `limlres' `accuweight', nocons
        mat `sebliml' = `limlres2'*invsym(`S'*`limlden')
        mata `sebliml' = st_matrix("`sebliml'")
        mata `sebliml' = sqrt(`sebliml')
        mata st_numscalar("sebliml",`sebliml')
        mata st_local("sebliml",strofreal(`sebliml'))
    }
    else {
        *qui avar (`limlres') (`xhatl') `addweight', `vcet' noconstant
        qui avar (`limlres') (`Zs'*) `addweight', `vcet' noconstant
        mata `Wl' = st_matrix("r(S)")
        mata `sebliml' = sqrt(invsym(`limlden')*`pihat''`Wl'*`pihat'*invsym(`limlden')/`S')
        mata st_numscalar("sebliml",`sebliml')
        mata st_local("sebliml",strofreal(`sebliml'))
    }

   
  
    *CP 10/28/2020: allow for noconstant
	if ( `noconst_flag'==1) {
	mata `Omega'=st_matrix("`Omega'")/(`S'-`K'-`L')
			}
			else{
	mata `Omega'=st_matrix("`Omega'")/(`S'-`K'-`L'-1)
	}

	* W matrix is the variance-covariance matrix adjusted by the degree of freedom given the specified vce type.
	qui avar (`res2' `res1') (`Zs'*) `addweight', `vcet' noconstant
	
	*CP 10/28/2020: allow for noconstant
	if ( `noconst_flag'==1) {
	   mata `W'=(st_matrix("r(S)")*`S'/(`S'-`K'-`L'))/`clustdfadj'
	}
	else{
	   mata `W'=(st_matrix("r(S)")*`S'/(`S'-`K'-`L'-1))/`clustdfadj'
	}
	
	* omega_22 is the submatrix of omega.
	mata `omega_22'=`Omega'[2,2]
	
	* W_1, W_2 and W_12 are submatrices of W.
	mata `W_1'=`W'[1::`K',1..`K']
	mata `W_12'=`W'[1::`K',`K'+1..2*`K']
	mata `W_2'=`W'[1+`K'::2*`K',1+`K'..2*`K']

    mata `W_2mh' = matpowersym(`W_2',-0.5)
    mata `trW_1s' = trace(`W_2mh'*`W_1'*`W_2mh')
    mata `W_12s' = `W_2mh'*`W_12'*`W_2mh'
    mata `trW_12s' = trace(`W_12s')
    
    mata `eigs' = symeigenvalues((`W_12s'+`W_12s'')/2)
	mata `mineig' = rowmin(`eigs')
	mata `maxeig' = rowmax(`eigs')

    mata `W_2inv' = invsym(`W_2')
    mata `bgmmf' = invsym(`pihat''`W_2inv'*`pihat')*(`pihat''`W_2inv'*`dhat')
    mata st_numscalar("bgmmf",`bgmmf')
    mata st_local("bgmmf",strofreal(`bgmmf'))
    qui gen `resgmmf' = `depvartemp' - `endotemp'*bgmmf
    
    qui avar (`resgmmf') (`Zs'*) `addweight', `vcet' noconstant

/*
    if ( `noconst_flag'==1) {
	   mata `S_1gmmf'=(st_matrix("r(S)")*`S'/(`S'-`K'-`L'))/`clustdfadj'
	}
	else{
	   mata `S_1gmmf'=(st_matrix("r(S)")*`S'/(`S'-`K'-`L'-1))/`clustdfadj'
	}
*/

    mata `S_1gmmf'=st_matrix("r(S)")
    
    mata `sebgmmf' = sqrt(invsym(`pihat''`W_2inv'*`pihat')*(`pihat''`W_2inv'*`S_1gmmf'*`W_2inv''`pihat')*invsym(`pihat''`W_2inv'*`pihat')/`S')
    mata st_numscalar("sebgmmf",`sebgmmf')
    mata st_local("sebgmmf",strofreal(`sebgmmf'))
  
   
    mata `btsls' = invsym(`pihat''`pihat')*(`pihat''`dhat')
    mata st_numscalar("btsls",`btsls')
    mata st_local("btsls",strofreal(`btsls'))
    qui gen `restsls' = `depvartemp' - `endotemp'*btsls
    qui avar (`restsls') (`Zs'*) `addweight', `vcet' noconstant
    mata `S_1tsls' = st_matrix("r(S)")
    mata `sebtsls' = sqrt((invsym(`pihat''`pihat')*(`pihat''`S_1tsls'*`pihat')*invsym(`pihat''`pihat'))/`S')
    mata st_numscalar("sebtsls",`sebtsls')
    mata st_local("sebtsls",strofreal(`sebtsls'))


    mata `F' = `S'*`pihat''`pihat'/(`K'*`omega_22')
    mata st_numscalar("F", `F')
    
    mata `F_rob' = `S'*(`pihat''invsym(`W_2')*`pihat')/`K'
    mata st_numscalar("F_rob", `F_rob')

	
	* Compute effective F estimates and save it as scalar.
	*** PH FIX 8/17/18: adjustment for clustering d.f.
	mata `F_eff'=`S'*`pihat''`pihat'/trace(`W_2')
	mata st_numscalar("F_eff", `F_eff')
	

/*Compute BTSLS and BLIML values.*/

	tempvar BTSLS BLIML BTSLS_start BLIML_start BGMMf_start BGMMf X_GMMf BGMMf_start
	
	*Find initial values for BTSLS maximization by calling program BTSLS_start.
	qui mata `BTSLS_start'=BTSLS_start(`Omega',`W_1',`W_12',`W_2',`eps')
	mata st_local("BTSLS_start", strofreal(`BTSLS_start'*(0,1)'))
	
	*Calculate BTSLS by calling optimization program BTSLS.
	qui mata `BTSLS'=BTSLS(`Omega',`W_1',`W_12',`W_2',`BTSLS_start')
	mata st_local("BTSLS", strofreal(`BTSLS'))
	
	*Find initial values for BLIML maximization by calling program BLIML_start.
	qui mata `BLIML_start'=BLIML_start(`W_1',`W_12',`W_2',`Omega',`eps')
	mata st_local("BLIML_start", strofreal(`BLIML_start'*(0,1)'))
	
	*Calculate BLIML by calling optimization function BLIML.
	qui mata `BLIML'=BLIML(`W_1',`W_2',`W_12',`Omega',`BLIML_start')
	mata st_local("BLIML", strofreal(`BLIML'))

    mata `BGMMf_start' = `bgmmf'
     
  	*Calculate BGMMf by calling optimization function BGMMf.
    mata `BGMMf' = J(4,1,0)
	cap qui mata `BGMMf'[1]=BGMMf(`Omega',`trW_1s',`trW_12s',`mineig',`K',"max",`BGMMf_start')
    cap qui mata `BGMMf'[2]=BGMMf(`Omega',`trW_1s',`trW_12s',`maxeig',`K',"max",`BGMMf_start')
    cap qui mata `BGMMf'[4]=BGMMf(`Omega',`trW_1s',`trW_12s',`maxeig',`K',"min",`BGMMf_start')
    cap qui mata `BGMMf'[3]=BGMMf(`Omega',`trW_1s',`trW_12s',`mineig',`K',"min",`BGMMf_start')     
    
    mata `BGMMf' = colmax(`BGMMf') 
    mata st_local("BGMMf", strofreal(`BGMMf'))
  
  
	
/*Compute and save critical values and effective degree of freedom.*/

	tempvar CVNrt 
	
/*Simplified TSLS*/
	
	*x=20   tau=5%
	mata  `CVNrt'=cpatnaikgen(`W_2', `level' , 20,`n2' )
	mata st_local("c_simple_5", strofreal(`CVNrt'*(1,0)')) 
	mata st_local("EK_simp_5",  strofreal(`CVNrt'*(0,1)'))
	
	*x=10   tau=10%
	mata  `CVNrt'=cpatnaikgen(`W_2', `level' , 10,`n2' )
	mata st_local("c_simple_10", strofreal(`CVNrt'*(1,0)'))
	mata st_local("EK_simp_10",  strofreal(`CVNrt'*(0,1)'))
	
	*x=5   tau=20%
	mata  `CVNrt'=cpatnaikgen(`W_2', `level' , 5,`n2' )
	mata st_local("c_simple_20", strofreal(`CVNrt'*(1,0)'))
	mata st_local("EK_simp_20",  strofreal(`CVNrt'*(0,1)'))
	
	*x=3.33  tau=30%
	mata  `CVNrt'=cpatnaikgen(`W_2', `level' , 3.33 ,`n2')
	mata st_local("c_simple_30", strofreal(`CVNrt'*(1,0)'))
	mata st_local("EK_simp_30",  strofreal(`CVNrt'*(0,1)'))

/*Non-simplified TSLS critical values.*/

	tempvar X_TSLS CVNrt_T 
	
	*X_TSLS is the input noncentrality parameter for computing critical value.
	*x=20  tau=5%
	mata `X_TSLS'=20*`BTSLS'                               
	mata  `CVNrt_T'=cpatnaikgen(`W_2', `level' , `X_TSLS',`n2' )
	mata st_local("X_TSLS_5", strofreal(`X_TSLS'))           
	mata st_local("c_TSLS_5", strofreal(`CVNrt_T'*(1,0)'))
	mata st_local("EK_TSLS_5", strofreal(`CVNrt_T'*(0,1)'))
	
	*x=10 tau=10%  
	mata `X_TSLS'=10*`BTSLS' 
	mata  `CVNrt_T'=cpatnaikgen(`W_2', `level' , `X_TSLS',`n2')
	mata st_local("X_TSLS_10", strofreal(`X_TSLS'))
	mata st_local("c_TSLS_10", strofreal(`CVNrt_T'*(1,0)'))
	mata st_local("EK_TSLS_10", strofreal(`CVNrt_T'*(0,1)'))
	
	*x=5   tau=20%  
	mata `X_TSLS'=5*`BTSLS' 
	mata  `CVNrt_T'=cpatnaikgen(`W_2', `level' , `X_TSLS',`n2')
	mata st_local("X_TSLS_20", strofreal(`X_TSLS'))
	mata st_local("c_TSLS_20", strofreal(`CVNrt_T'*(1,0)'))
	mata st_local("EK_TSLS_20", strofreal(`CVNrt_T'*(0,1)'))
	
	*x=3.33  tau=30%  
	mata `X_TSLS'=3.33*`BTSLS' 
	mata  `CVNrt_T'=cpatnaikgen(`W_2', `level' , `X_TSLS',`n2')
	mata st_local("X_TSLS_30", strofreal(`X_TSLS'))
	mata st_local("c_TSLS_30", strofreal(`CVNrt_T'*(1,0)'))
	mata st_local("EK_TSLS_30", strofreal(`CVNrt_T'*(0,1)'))

/*Non-simplified LIML critical values.*/

	tempvar X_LIML CVNrt_L
	
	*X_LIML is the input noncentrality parameter for computing critical value.
	*x=20  tau=5%  
	mata `X_LIML'=20*`BLIML'                               
	mata  `CVNrt_L'=cpatnaikgen(`W_2', `level' , `X_LIML' ,`n2')
	mata st_local("X_LIML_5", strofreal(`X_LIML'))
	mata st_local("c_LIML_5", strofreal(`CVNrt_L'*(1,0)'))
	mata st_local("EK_LIML_5", strofreal(`CVNrt_L'*(0,1)'))
	
	*x=10   tau=10%  
	mata `X_LIML'=10*`BLIML' 
	mata  `CVNrt_L'=cpatnaikgen(`W_2', `level' , `X_LIML',`n2')
	mata st_local("X_LIML_10", strofreal(`X_LIML'))
	mata st_local("c_LIML_10", strofreal(`CVNrt_L'*(1,0)'))
	mata st_local("EK_LIML_10", strofreal(`CVNrt_L'*(0,1)'))
	
	*x=5   tau=20%  
	mata `X_LIML'=5*`BLIML' 
	mata st_local("X_LIML_20", strofreal(`X_LIML'))
	mata  `CVNrt_L'=cpatnaikgen(`W_2', `level' , `X_LIML',`n2' )
	mata st_local("c_LIML_20", strofreal(`CVNrt_L'*(1,0)'))
	mata st_local("EK_LIML_20", strofreal(`CVNrt_L'*(0,1)'))
	
	*x=3.33  tau=30%  
	mata `X_LIML'=3.33*`BLIML' 
	mata  `CVNrt_L'=cpatnaikgen(`W_2', `level' , `X_LIML',`n2' )
	mata st_local("X_LIML_30", strofreal(`X_LIML'))
	mata st_local("c_LIML_30", strofreal(`CVNrt_L'*(1,0)'))
	mata st_local("EK_LIML_30", strofreal(`CVNrt_L'*(0,1)'))

/*Non-simplified GMMf critical values.*/

	tempvar X_GMMf CVNrt_G 
	
	*X_GMMf is the input noncentrality parameter for computing critical value.
	*x=20  tau=5%
	mata `X_GMMf' = 20*`BGMMf'*`K'                               
	mata  `CVNrt_G' = invnchi2(`K',`X_GMMf',1-`level')/`K'
	mata st_local("X_GMMf_5", strofreal(`X_GMMf'))           
	mata st_local("c_GMMf_5", strofreal(`CVNrt_G'))
		
	*x=10 tau=10%  
	mata `X_GMMf' = 10*`BGMMf'*`K'                               
	mata  `CVNrt_G' = invnchi2(`K',`X_GMMf',1-`level')/`K'
	mata st_local("X_GMMf_10", strofreal(`X_GMMf'))           
	mata st_local("c_GMMf_10", strofreal(`CVNrt_G'))
	
	*x=5   tau=20%  
	mata `X_GMMf' = 5*`BGMMf'*`K'                               
	mata  `CVNrt_G' = invnchi2(`K',`X_GMMf',1-`level')/`K'
	mata st_local("X_GMMf_20", strofreal(`X_GMMf'))           
	mata st_local("c_GMMf_20", strofreal(`CVNrt_G'))
	
	*x=3.33  tau=30%  
	mata `X_GMMf' = 3.33*`BGMMf'*`K'                               
	mata  `CVNrt_G' = invnchi2(`K',`X_GMMf',1-`level')/`K'
	mata st_local("X_GMMf_30", strofreal(`X_GMMf'))           
	mata st_local("c_GMMf_30", strofreal(`CVNrt_G'))
	


/*Save returned values: Effective degrees of freedom K_eff and noncentrality parameter x.*/

	*Simplified noncentrality parameters
	return scalar K_eff_simp_30= `EK_simp_30'
	return scalar K_eff_simp_20= `EK_simp_20'
	return scalar K_eff_simp_10= `EK_simp_10'
	return scalar K_eff_simp_5=  `EK_simp_5'
	return scalar x_simp_30= 3.33
	return scalar x_simp_20= 5
	return scalar x_simp_10= 10
	return scalar x_simp_5=  20
	
	*Non-simplified LIML parameters
	return scalar K_eff_LIML_30= `EK_LIML_30'
	return scalar K_eff_LIML_20= `EK_LIML_20'
	return scalar K_eff_LIML_10= `EK_LIML_10'
	return scalar K_eff_LIML_5=  `EK_LIML_5'
	return scalar x_LIML_30= `X_LIML_30'
	return scalar x_LIML_20= `X_LIML_20'
	return scalar x_LIML_10= `X_LIML_10'
	return scalar x_LIML_5=  `X_LIML_5'
	
	*Non-simplified TSLS parameters
	return scalar K_eff_TSLS_30= `EK_TSLS_30'
	return scalar K_eff_TSLS_20= `EK_TSLS_20'
	return scalar K_eff_TSLS_10= `EK_TSLS_10'
	return scalar K_eff_TSLS_5=  `EK_TSLS_5'
	return scalar x_TSLS_30= `X_TSLS_30'
	return scalar x_TSLS_20= `X_TSLS_20'
	return scalar x_TSLS_10= `X_TSLS_10'
	return scalar x_TSLS_5=  `X_TSLS_5'

/*Save return values: Critical values.*/

	*Simplified TSLS
	return scalar c_simp_30= `c_simple_30'
	return scalar c_simp_20= `c_simple_20'
	return scalar c_simp_10= `c_simple_10'
	return scalar c_simp_5=  `c_simple_5'
	
	*Non-simplified LIML
	return scalar c_LIML_30= `c_LIML_30'
	return scalar c_LIML_20= `c_LIML_20'
	return scalar c_LIML_10= `c_LIML_10'
	return scalar c_LIML_5=  `c_LIML_5'
	
	*Non-simplified TSLS
	return scalar c_TSLS_30= `c_TSLS_30'
	return scalar c_TSLS_20= `c_TSLS_20'
	return scalar c_TSLS_10= `c_TSLS_10'
	return scalar c_TSLS_5=  `c_TSLS_5'

    *Non-simplified GMMf
	return scalar c_GMMf_30= `c_GMMf_30'
	return scalar c_GMMf_20= `c_GMMf_20'
	return scalar c_GMMf_10= `c_GMMf_10'
	return scalar c_GMMf_5=  `c_GMMf_5'
	

	
	* F_eff and inputs
    return scalar F_rob = F_rob
	return scalar F_eff= F_eff
    return scalar btsls = btsls
    return scalar sebtsls = sebtsls
    return scalar bliml = bliml
    return scalar sebliml = sebliml
    return scalar kappa = kappa
    return scalar bgmmf = bgmmf
  	return scalar sebgmmf = sebgmmf
    return scalar eps= `eps'   		
	return scalar level=`level'
	return scalar n2=`n2'
    
	* Number of observations and numbers of instruments
	return scalar N= `e(N)'
	return scalar K= `exexog_num'
			
/* Generate output table.*/

	local b=`level'*100
	display ""
	display in g "Montiel-Pflueger robust weak instrument test"

	display "{txt}{hline 44}"
   	display in g "btsls:	                " in y  %9.4fc btsls
    display in g "sebtsls:                  " in y  %9.4fc sebtsls
    display in g "bliml:                    " in y  %9.4fc bliml     
    display in g "sebliml:                  " in y  %9.4fc sebliml    
    display in g "kappa:                    " in y  %9.4fc kappa
    display in g "Non-Robust F statistic:   " in y  %9.3fc F
   	display in g "Effective F statistic:	" in y  %9.3fc F_eff
	display in g "Confidence level alpha:	" in y  "       `b'%"
	display "{txt}{hline 44}"

	display ""

	display "{txt}{hline 44}"
	display "Critical Values	" "	TSLS" "	     " "LIML"
	display "{txt}{hline 44}"
	display "% of Worst Case Bias" 	
	display in g "tau=5%		     "  in y %9.3fc  `c_TSLS_5' " " %9.3fc `c_LIML_5'
	display in g "tau=10%		     "  in y %9.3fc `c_TSLS_10' " " %9.3fc `c_LIML_10'
	display in g "tau=20%		     "  in y %9.3fc `c_TSLS_20' " " %9.3fc `c_LIML_20'
	display in g "tau=30%		     "  in y %9.3fc `c_TSLS_30' " " %9.3fc `c_LIML_30'
	display "{txt}{hline 44}"

	display "{txt}{hline 44}"
	display in g "bgmmf:	                " in y  %9.4fc bgmmf
    display in g "sebgmmf:                  " in y  %9.4fc sebgmmf
    display in g "Robust F statistic:	    " in y  %9.3fc F_rob
 	display in g "Confidence level alpha:	" in y  "       `b'%"
	display "{txt}{hline 44}"

	display ""

    display "{txt}{hline 44}"
	display "Critical Values	" "	GMMf" "
	display "{txt}{hline 44}"
	display "% of Worst Case Bias" 	
	display in g "tau=5%		     "  in y %9.3fc  `c_GMMf_5'
	display in g "tau=10%		     "  in y %9.3fc `c_GMMf_10'
	display in g "tau=20%		     "  in y %9.3fc `c_GMMf_20'
	display in g "tau=30%		     "  in y %9.3fc `c_GMMf_30'
	display "{txt}{hline 44}"


/*Unhold ereturn result*/

	_estimates unhold  `myereturn'
	
end

/*MATA CODE PART*/
mata
version 10
	
/*cpatnaikgen function computes critical values from variance-covariance matrix W_2, significance level a, 
noncentrality parameter x and non-central inverse F distribution parameter n2.*/
	real matrix cpatnaikgen(W_2,	/*
						*/	a,  	/*
						*/  x,  	/*
						*/  n2){
		W2=eigenvalues(W_2)
		W2=Re(W2)
		W2=W2/sum(W2)
		W2=sort(W2',1)'
		variance=2*sum(W2:*W2)+4*x*max(W2)
		EK=2*(1+2*x)/variance
		Delta=EK*x
		
        /*Use inverse F function to approximate the non-central F distribution by setting n2 large enough.*/
		if (st_numscalar("c(stata_version)")>=13.1) {
			st_numscalar("r(EK)",EK)
			st_numscalar("r(Delta)",Delta)
			st_numscalar("r(a)",a)
			stata("version `c(stata_version)': local cvalue = invnchi2(r(EK), r(Delta), 1-r(a))/r(EK)")
			cvalue = strtoreal(st_local("cvalue"))
		}
		
		else {
			cvalue=invnFtail(EK,n2,Delta,a)
		}
		RTN=cvalue,EK
		return (RTN)
	}


/*Bmaxfunction computes maximal Nagar bias for given values of beta and variance-covariance matrix W. 
This step of the maximization is analytic.*/
	real matrix Bmaxfunction (real scalar beta, 	/*
						*/	  real matrix Omega,    /*
                        */    real matrix W_1,	    /*
						*/	  real matrix W_12, 	/*
						*/	  real matrix W_2){
		S_2=W_2
		S_12=W_12-beta*W_2
		S_1=W_1-2*beta*W_12+beta^2*W_2
		L=eigenvalues(0.5*S_12+0.5*S_12')
		L=Re(L)
		mineig=rowmin(L)
		maxeig=rowmax(L)
		BM = sqrt((Omega[1,1]-2*beta*Omega[1,2]+beta^2*Omega[2,2])/Omega[2,2])
		B=abs((trace(S_12)-2*mineig)/trace(W_2))/BM
		C=abs((trace(S_12)-2*maxeig)/trace(W_2))/BM
		B=max(B\C)
		X=B,beta
		return (X)
	}

/*BTSLS_start computes the start point of maximization for given variance-covariance matrix W 
and deviation parameter eps.*/
	real matrix BTSLS_start (real matrix Omega,  /*
                          */ real matrix W_1,    /*
						  */ real matrix W_12,   /*
						  */ real matrix W_2,    /*
						  /*eps is the maximal fractional deviation of B from the limit
						  as beta goes to infinity; default value eps=0.001.*/
						  */ eps){  
			  					  
		L=eigenvalues(W_2)
		L=Re(L)
		eigmin=rowmin(L)
		
		/*LimitB is the limit of Bmaxfunction as beta approaches +/- infinity. We then find a range [-betastart,betastart] 
		 so that Bmaxfunction(betastart) and Bmaxfunction(-betastart) are within a fraction eps of LimitB.*/
		LimitB=1-2*eigmin/trace(W_2)
		betastart=0
		
		/*Evaluate Bmaxfunction at points=10000 equally spaced values.*/
		points=10000   
		val=max(abs(abs(Bmaxfunction(betastart, Omega, W_1, W_12, W_2)*(1,0)'/LimitB)-1)\   ///
					 abs(abs(Bmaxfunction(-betastart, Omega, W_1, W_12, W_2)*(1,0)'/LimitB)-1))
					 
		/*If val is within a fraction eps of LimitB, set s=0*/
		if (val<eps | val==eps){   
			epserror=0
			s=0
		}
		
		else {
		/*If Val differs from LimitB by more than a fraction eps, increase betastart until condition is satisfied.*/
			while (val>eps){  
				val=max(abs(abs(Bmaxfunction(betastart, Omega, W_1, W_12, W_2)*(1,0)'/LimitB)-1)\   ///
						abs(abs(Bmaxfunction(-betastart, Omega, W_1, W_12, W_2)*(1,0)'/LimitB)-1))
				betastart=betastart+1
			}
			epserror=0

			t=-betastart 
			s=t
			
			/*Evaluate Bmaxfunction at points equally spaced values and return the maximizing beta.*/
			BplotTSLS=Bmaxfunction(t, Omega, W_1, W_12, W_2)*(1,0)'
			while (t<=betastart){
				t=t+betastart/points
				BTSLS=Bmaxfunction(t, Omega, W_1, W_12, W_2)*(1,0)'
				BplotTSLS=max(BTSLS\BplotTSLS)
				if (BTSLS<BplotTSLS){
					s=s
				}
				else{
					s=t
				}
			}
		}
	x=epserror,s
	return(x)
	}

/*Define objective function myfun to be maximized to compute BTSLS given beta, variance-covariance matrix W,
BTSLS_start and returns value y. todo is required action by evaluator(). g, and h are values to be returned.*/
	void myfun (todo,   		/*
			*/  beta,   		/*
            */  Omega,           /*
			*/	 W_1,   		/*
			*/	 W_2,         	/*
			*/	 W_12,   		/*
			*/	 BTSLS_start,   /*
			*/	 y,             /*
			*/	 g,             /*
			*/	 h)
	{
		S_2=W_2
		S_12=W_12-beta*W_2
		S_1=W_1-2*beta*W_12+beta^2*W_2
		L=eigenvalues(0.5*S_12+0.5*S_12')
		L=Re(L)
		mineig=rowmin(L)
		maxeig=rowmax(L)
        BM = sqrt((Omega[1,1]-2*beta*Omega[1,2]+beta^2*Omega[2,2])/Omega[2,2])
		B=abs((trace(S_12)-2*mineig)/trace(W_2))/BM
		C=abs((trace(S_12)-2*maxeig)/trace(W_2))/BM
		y=max(B\C)
	}

/*BTSLS computes the maximizing value BTSLS given function myfun, variance-covariance matrix W 
and start point BTSLS_Start by Nelder-Mead technique.*/
	real matrix BTSLS (real matrix Omega,   /*
                 */     real matrix W_1, 	/*
				 */		real matrix W_12, 	/*
				 */		real matrix W_2, 	/*
				 */		BTSLS_start){
				 
		/*Define initial values as S*/
		S=optimize_init()
		
		/*Define evaluator myfun.*/
		optimize_init_evaluator(S, &myfun ())
		
		/*Set initial value of S to be BTSLS_start.*/
		optimize_init_params(S,BTSLS_start)
		
		/*Specify initial argument used.*/
		optimize_init_argument(S, 1, Omega)
		optimize_init_argument(S, 2, W_1)
		optimize_init_argument(S, 3, W_2)
		optimize_init_argument(S, 4, W_12)
		optimize_init_argument(S, 5, BTSLS_start)
		
		/*Nelder-Mead method is used as optimization techinique.*/
		optimize_init_technique(S, "nm")
		
		/*Set the value of delta to be 0.5, along with the initial parameter S 
		to build simplex required by technique "nm".*/
		optimize_init_nmsimplexdeltas(S,0.5)
		beta_TSLS=optimize(S)
		
		/*Return optimized value as BTSLS.*/
		BTSLS=optimize_result_value(S)
		return (BTSLS)
	}


/*BmaxLIML computes maximal Nagar bias for a given value of beta,
 variance-covariance matrix W and Omega.*/
	real matrix BmaxLIML (real scalar beta,   /*
					 */    real matrix W_1,   /*
					 */	   real matrix W_12,  /*
					 */	   real matrix W_2,   /*
					 /*Omega is the reduced form variance-covariance matrix.*/
					 */    real matrix Omega){    
		S_2=W_2
		S_12=W_12-beta*W_2
		S_1=W_1-2*beta*W_12+beta^2*W_2
		om_1=Omega[1,1]
		om_12=Omega[1,2]
		om_2=Omega[2,2]
		sig_12=om_12-beta*om_2
		sig_1=om_1-2*beta*om_12+beta^2*om_2
		Matrix=2*S_12-sig_12/sig_1*S_1
		Matrix=0.5*(Matrix+Matrix')
		L=eigenvalues(Matrix)
		L=Re(L)
		mineig=rowmin(L)
		maxeig=rowmax(L)
		B=1/sqrt(trace(S_2)*trace(S_1))*(trace(S_12)-sig_12/sig_1*trace(S_1)-mineig)
		B=abs(B)
		B=max(B\abs(1/sqrt(trace(S_2)*trace(S_1))*(trace(S_12)-sig_12/sig_1*trace(S_1)-maxeig)))
		X=B,beta
		return (X)
	}

/*BLIML_start returns start point of maximization given variance-covariance matrix W , Omega 
and deviation parameter eps using BmaxLIML.*/
	real matrix BLIML_start (real matrix W_1,    /*
						*/     real matrix W_12,   /*
						*/     real matrix W_2,    /*
						*/     real matrix Omega,  /*
						*/     eps){
							
		L=eigenvalues(W_2)
		L=Re(L)
		eigmin=rowmin(L)
		eigmax=rowmax(L)
		
		/*LimitB is the limit of BmaxLIML as beta approaches +/- infinity. */
		LimitB=eigmax/trace(W_2)
		betastart=0
		
		/*Points is the number of equally space values. Default value is 10000*/
		points=10000  
		val=max(abs(BmaxLIML(betastart, W_1, W_12, W_2, Omega)*(1,0)'/LimitB-1)\   ///
				abs(BmaxLIML(-betastart, W_1, W_12, W_2, Omega)*(1,0)'/LimitB-1))
				
		/*If val is within a fraction eps of LimitB, set s=0*/
		if (val<eps |val==eps){ 
			epserror=0
			s=0
		}
		
		else{
		/*If val is outside a fraction eps of LimitB, increase betastart until condition is satisfied.*/
			while (val>eps){  
				val=max(abs(BmaxLIML(betastart, W_1, W_12, W_2, Omega)*(1,0)'/LimitB-1)\   ///
						abs(BmaxLIML(-betastart, W_1, W_12, W_2, Omega)*(1,0)'/LimitB-1))
				betastart=betastart+1
			}
			
			epserror=0
			t=-betastart 
			s=t
			
			/*Evaluate BmaxLIML at points equally spaced values and return the maximizing value.*/
			BplotLIML=BmaxLIML(t, W_1, W_12, W_2, Omega)*(1,0)'
			while (t<=betastart){
				t=t+betastart/points
				BLIML=BmaxLIML(t, W_1, W_12, W_2,Omega)*(1,0)'
				BplotLIML=max(BLIML\BplotLIML)
				if (BLIML<BplotLIML){
					s=s
				}
				else{
					s=t
				}
			}
		}
	x=epserror, s
	return(x)
	}

/*Define objective function myfun2 to maximize for computing BLIML given beta, variance-covariance matrix W and Omega,
BLIML_start and returns value y. todo is required action by evaluator(). g, and h are values to be returned */
	void myfun2(todo,       	/*
		  */     beta,  		/*
		  */     W_1,   		/*
		  */     W_2,   		/*
		  */     W_12,  		/*
		  */     Omega,  		/*
		  */     BLIML_start,   /*
		  */     y,             /*
		  */     g,  			/*
		  */     h)
		{
			S_2=W_2
			S_12=W_12-beta*W_2
			S_1=W_1-2*beta*W_12+beta^2*W_2
			om_1=Omega[1,1]
			om_12=Omega[1,2]
			om_2=Omega[2,2]
			sig_12=om_12-beta*om_2
			sig_1=om_1-2*beta*om_12+beta^2*om_2
			Matrix=2*S_12-sig_12/sig_1*S_1
			Matrix=0.5*(Matrix+Matrix')
			L=eigenvalues(Matrix)
			L=Re(L)
			mineig=rowmin(L)
			maxeig=rowmax(L)
			B=abs(1/sqrt(trace(S_2)*trace(S_1))*(trace(S_12)-sig_12/sig_1*trace(S_1)-mineig))
			y=max(B\abs(1/sqrt(trace(S_2)*trace(S_1))*(trace(S_12)-sig_12/sig_1*trace(S_1)-maxeig)))
		}

/*BLIML computes the maximized value BLIML given function myfun2, variance-covariance matrix W 
and Omega and start point BLIML_Start by using the Nelder-Mead technique.*/
	real matrix BLIML (real matrix W_1,		 /*
				*/     real matrix W_2, 	/*
				*/	   real matrix W_12,	/*
				*/	   real matrix Omega,	/*
				*/	   BLIML_start){
				
		/*Define initial values S.*/
		S=optimize_init()
		
		/*Define evaluator myfun2.*/
		optimize_init_evaluator(S, &myfun2 ())
		
		/*Set initial values of S to be BLIML_start.*/
		optimize_init_params(S, BLIML_start)
		
		/*Specify initial argument used*/
		optimize_init_argument(S, 1, W_1)
		optimize_init_argument(S, 2, W_2)
		optimize_init_argument(S, 3, W_12)
		optimize_init_argument(S, 4, Omega) 
		optimize_init_argument(S, 5, BLIML_start) 
		
		/*Nelder-Mead method is used as optimization techinique.*/
		optimize_init_technique(S, "nm")
		
		/*Set the values of delta to be 0.5, along with the initial parameter S.
		to build simplex required by technique "nm".*/
		optimize_init_nmsimplexdeltas(S,0.5)
		beta_BLIML=optimize(S)
		
		/*Return optimize value as BLIML.*/
		BLIML=optimize_result_value(S)
		return (BLIML)
	}

/*Define objective function myfun3 to be maximized to compute BGMMf given beta, variance-covariance matrix W,
 and returns value y. todo is required action by evaluator(). g, and h are values to be returned.*/
	void myfun3(todo,   		/*
			*/  beta,   		/*
            */  Omega,           /*
			*/	trW_1s,   		/*
			*/	trW_12s,      	/*
            */  eig,            /*
            */  K,              /*
			*/	y,              /*
			*/	g,              /*
			*/	h)
	{
		num = (trW_12s-2*eig-(K-2)*beta)/K
		den = sqrt((Omega[1,1]-2*beta*Omega[1,2]+beta^2*Omega[2,2])/Omega[2,2])
		y = num/den

	}

/*BGMMf computes the maximizing value BGMMf given function myfun3, variance-covariance matrix Ws 
and start point BGMMf_Start*/
	real matrix BGMMf ( real matrix Omega,       /*
                 */     real matrix trW_1s, 	/*
				 */		real  matrix trW_12s, 	/*
				 */		real matrix eig, 	/*					 
				 */		real matrix K,     /*
                 */     maxmin,            /*
                 */     BGMMf_start){
				 
		/*Define initial values as S*/
		S=optimize_init()
		
		/*Define evaluator myfun3.*/
		optimize_init_evaluator(S, &myfun3 ())
		
		/*Set initial value of S to be BGMMf_start.*/
		optimize_init_params(S,BGMMf_start)
		
		/*Specify initial argument used.*/
		optimize_init_argument(S, 1, Omega)
		optimize_init_argument(S, 2, trW_1s)
		optimize_init_argument(S, 3, trW_12s)
		optimize_init_argument(S, 4, eig)
		optimize_init_argument(S, 5, K)

        optimize_init_conv_warning(S, "off") 
		optimize_init_which(S,maxmin)
 				       
		beta_GMMf = optimize(S)

        /*Return optimized value as BGMMf.*/
		BGMMf = abs(optimize_result_value(S))
		return (BGMMf)
	}


end

