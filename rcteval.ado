*! version 0.3  Chris J. Kennedy 2013-04-09.
program define rcteval
version 12.0

/*
* Example syntax: 
rcteval `dv' `assignment', covars(`covars') subgroups(`subgroups') review(`assignment_review') balanceby(`balance_var')


TODO:
- base() -> specify the base outcome when analyzing the assignments.
- For multiple assignment options need to check balance across every possible assignment base.
- Show table of covariate means across assignments.
- Option or one-tailed or two-tailed p-values.
- control() -> specify the assignment variable that is the control arm. presumed to be the lowest value if not specified.

*/

syntax varlist(min=2 max=2 numeric) [if] [in] [, COVARS(varlist fv) SUBGroups(varlist) REView(varlist) BALanceby(varlist) SKIPBALance /* CONtrol(integer) */ MODel(string) CLuster(varlist)]

qui: marksample touse // Exclude observations that do not meet the IF or IN criteria (if specified).

tokenize `varlist'
local dv `1'
local assignment `2'

* Time the analysis
timer clear 83
timer on 83

* Review the distribution of the assignment variable.
qui su `assignment' if `touse'
* Set the control value to be the lowest value of the assignment term.
local control = r(min)

dis "Review assignment breakdown within the test universe."
foreach var in `review' {
	tab `var' `assignment' if `touse', row chi2
}

* Run a balance check unless we were told to skip it.
if "`skipbalance'" == "" {
	dis "Balance check on `covars'"
	* Explicitly set the base value of the mlogit to use the control variable.
	if "`balanceby'" != "" {
		dis "Checking balance across separate values of `balanceby'."
		tab `balanceby' `assignment', row chi2
		local balance_bysort = "bysort `balanceby': "
	}
	`balance_bysort'mlogit `assignment' `covars' if `touse', nolog base(`control')
	* TODO: parse results and report variables that have imbalance, plus overall balance check p-value.
	* TODO: give a warning if any imbalanced variable is not included in the covars list.
}

dis _n(3) "a. Simple breakdown -- "
tab `assignment' `dv' if `touse', row chi2
* Save the number of records for which we have the dependent variable, for future error-checking purposes.
local total_records = r(N)

dis _n(3) "b. Basic analysis -- "
reg `dv' ib`control'.`assignment'  if `touse', cluster(`cluster')

dis _n(3) "c. With covariates (`covars') -- " 
reg `dv' ib`control'.`assignment' `covars' if `touse', cluster(`cluster')

* Check that we didn't lose any records with the covariates.
local analyzed_records = e(N)
if e(N) != `total_records' {
	dis as error "DV defined on `total_records' records, but analysis only includes `analyzed_records' records."
	dis as error "WARNING: Covariates may be missing on a portion of the analyzed records."
}

dis _n(3) "d. Logit with covars -- (`covars')"
logit `dv' ib`control'.`assignment' `covars' if `touse', nolog or cluster(`cluster')
* dis "Marginal effects:"
* prchange `assignment', rest(mean) brief fromto

dis _n(3) "e. Subgroup analysis -- (`subgroups')"
foreach var in `subgroups' {
	* Confirm that there are at least two levels of this subgroup.
	qui tab `var' if `touse'
	if r(r) < 2 {
		dis "Skipping subgroup `var' because it does not have 2 or more levels."
		continue
	}
	/** Skip this part to keep the log concise.
	dis _n "Interaction test for `var', no covariates -- "
	logit `dv' `var'##i.`assignment' if `touse', nolog or cluster(`cluster')
	* We allow the distribution of interaction term values to be left as observed in the dataset rather than assumed to be equal (as balanced).
	contrast `var'##i.`assignment', asobserved
	* Note: we do not use the p-value from this version but we do report the results for posterity.
	*/
	
	dis _n(3) "Interaction test for `var', with covariates -- "
    logit `dv' `var'##ib`control'.`assignment' `covars' if `touse', nolog or cluster(`cluster')
    * estimates store model_interact
    * We allow the distribution of interaction term values to be left as observed in the dataset rather than assumed to be equal (as balanced).
	contrast `var'##i.`assignment', asobserved
	matrix p_values = r(p)
	* Interaction term will be the third column in the first row.
	local interact_test = p_values[1, 3]
	
	* Run a reduced form model for comparison. Make sure that the interaction term is defined for all records.
	* logit `dv' ib`control'.`assignment' `covars' if `touse' & `var' != ., nolog or cluster(`cluster')
	* Also compute a LR test to compare the model with interacting treatment to one without.
	* This test is displayed for comparison to the contrast command that includes covariates, but is not explicitly checked in the algorithm.
	* See Long & Freese (2006). "Regression models for Categorical Indepenent Variables Using Stata". 2nd Edition. Chapter 9, p. 424.
	* lrtest . model_interact -> can't do with clustered SEs.

	* TODO: p-value cut-off should be a parameter and/or we should control for multiple comparisons.
	* TODO: ideally we would also check if the subgroup effect is clinically significant (at least 25% different than main effect).
		* (See http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2813449/ )
	if `interact_test' < 0.05 {
		dis "** Found significant interaction effect for `var' (p = " as result %06.4f round(`interact_test', .0001) as text ")."
		dis "By `var', baseline returns -- "
		bysort `var': tab `assignment' `dv' if `touse', row chi2
		/** Skip this to keep the log concise.
		dis "By `var', basic analysis --"
		bysort `var': reg `dv' ib`control'.`assignment' if `touse', cluster(`cluster')
		*/
		dis "By `var' with covariates (`covars') --"
		bysort `var': reg `dv' ib`control'.`assignment' `covars' if `touse', cluster(`cluster')
	}
	else {
		dis "** Did not find significant interaction effect for `var' (p = " as result %06.4f round(`interact_test', .0001) as text ")."
	}
	dis _n
}
* TODO: print summary of interaction effects found with p-values, and list of all variables tested.

* Report our runtime for this experimental evaluation.
timer off 83
qui timer list 83
dis _n(2) "rtcteval.ado analysis time:"
dis "Seconds: " as result r(t83) as text ". Minutes: " as result round(r(t83) / 60, 0.01) as text "."
dis "-----------------------" _n(3) 

end
