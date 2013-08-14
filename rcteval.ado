*! version 0.3  Chris J. Kennedy 2013-04-09.
program define rcteval
version 12.0

/*
* Example syntax: 
rcteval `dv' `assignment', covars(`covars') subgroups(`subgroups') continuous(`continuous subgroups') review(`assignment_review') balanceby(`balance_var')

tiles: ordinal categorical variables that should be converted to tiles, then have a specific HTE algorithm used.
  - Try a 5-tile first at joint p < 0.1 or any individual p < 0.1, then if any is significant run 3-tile and 7-tile for comparison.

Continuous(varlist): specify subgroups that are continuous variable and should receive special treatment. By default subgroups are presumed to be categorical.

Clean: if specified, delete the auto-created tiles for continuous variables. Otherwise leave them in the dataset.
ExpName: name of the experiment, for use in graph filenames.

TODO:
- base() -> specify the base outcome when analyzing the assignments.
- For multiple assignment options need to check balance across every possible assignment base.
- Show table of covariate means across assignments.
- Option or one-tailed or two-tailed p-values.
- control() -> specify the assignment variable that is the control arm. presumed to be the lowest value if not specified.

*/

syntax varlist(min=2 max=2 numeric) [if] [in] [, COVARS(varlist fv) SUBGroups(varlist) CONTinuous(varlist) REView(varlist) BALanceby(varlist) SKIPBALance /* CONtrol(integer) */ MODel(string) EXPname(string) CLuster(varlist) CLEAN]

* Exclude observations that do not meet the IF or IN criteria (if specified).
qui: marksample touse

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
qui tab `assignment' if `touse'
* Save the number of assignment levels for later usage.
local assignment_levels = r(r)

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



*******
* Pre-process subgroups en masse and remove from current analysis if they don't have two or more levels in this sample.
**
local subgroups_clean ""
local subgroups_skip ""
foreach var in `subgroups' {
	* Confirm that there are at least two levels of this subgroup.
	qui tab `var' if `touse'
	if r(r) < 2 {
		* Skip this subgroup.
		local blank " "
		* Don't add a space if this is the first skipped subgroup.
		if "`subgroups_skip'" == "" {
			local blank = ""
		}
		local subgroups_skip "`subgroups_skip'`blank'`var'"	
	}
	else {
		* Keep this subgroup.
		local blank " "
		* Don't add a space if this is the first clean subgroup.
		if "`subgroups_clean'" == "" {
			local blank = ""
		}
		local subgroups_clean "`subgroups_clean'`blank'`var'"
	}
}

dis _n(3) "e. Subgroup analysis (`subgroups_clean')"
* Only display the skip line if we skipped at least 1 subgroup.
local n_skip: list sizeof subgroups_skip
if `n_skip' > 0 {
	dis "Subgroups skipped (`n_skip'): `subgroups_skip'."
}




foreach var in `subgroups_clean' {
	* Confirm that there are at least two levels of this subgroup.
	qui tab `var' if `touse'
	* Save the number of levels for future usage.
	local group_levels = r(r)
	* Note: this section isn't needed anymore because we have already skipped these variables.
	if r(r) < 2 {
		dis "Skipping subgroup `var' because it does not have 2 or more levels."
		continue
	}
	
	* Check if this is a continuous subgroup.
	local is_cont: list var in continuous
	if `is_cont' == 1 {
		* This is a continuous variable.
		
		* Use the continuous prefix designator as part of Stata's factor variable handling.
		* This will help with the interaction test, et al.
		local cont_prefix = "c."
		
		* Create specific tiles for this experiment. Alternatively we could do this earlier for each continuous variable.
		* Drop our target variable name if it exists.
		cap drop etile_`var'
		xtile etile_`var' = `var' if `touse', n(5) 
	}

	/** Skip this part to keep the log concise.
	dis _n "Interaction test for `var', no covariates -- "
	logit `dv' `var'##i.`assignment' if `touse', nolog or cluster(`cluster')
	* We allow the distribution of interaction term values to be left as observed in the dataset rather than assumed to be equal (as balanced).
	contrast `var'##i.`assignment', asobserved
	* Note: we do not use the p-value from this version but we do report the results for posterity.
	*/
	
	dis _n(3) "Interaction test for `var', with covariates -- "
    logit `dv' `cont_prefix'`var'##ib`control'.`assignment' `covars' if `touse', nolog or cluster(`cluster')
    matrix results = r(table)
    
    * If the subgroup has only two levels or is continuous we can check the p-value directly from the regression.
    if `group_levels' == 2 | `is_cont' == 1 {
    	
 	    * Determine the p-value column to check for the interaction effect.
 	    * Generally this is the third p-value in the regression, for a binary assignment variable and binary interaction term.
	    * OLD -- local check_col = `assignment_levels' + `group_levels' + `assignment_levels' * `group_levels'
	    * Leaving the final 1 in here to make the calculation clearer.
	
	    if `is_cont' == 1 {
	    	local p_num = (`assignment_levels' - 1) + 2
	    }
	    else {
		    local p_num = (`assignment_levels' - 1) + (`group_levels' - 1) + 1
		}
	
	    matrix p_values = results["pvalue", 1...]
	    local num_columns = colsof(p_values)
	    local target_col = 0
	    * Loop through p-values and find the correct non-missing value.
	    forvalues i = 1/`num_columns' {
	   		if p_values[1, `i'] != . {
	   			* This is a valid p-value, so decrement our counter until we eventually hit zero.
				local --p_num
	    	}
	    	if `p_num' == 0 {
	    		* We are at the right p-value so save it.
	    		local interact_test = p_values[1, `i']
	    		* Stop looping.
	    		continue, break
	    	}
	    }
	    dis "Preliminary interaction test p-value: " as result %06.4f round(`interact_test', .0001) as text "."
	}
	else {
		* Subgroup has more than two levels, so first do a pairwise comparison of the coefficients.
		* pwcompare `var'#ib`control'.`assignment', asobs effects
		
		* <Insert loop that checks for significant comparsions.>
		
		* Now run a joint test of interaction.
    
   		* estimates store model_interact
	    * We allow the distribution of interaction term values to be left as observed in the dataset rather than assumed to be equal (as balanced).
		contrast `var'##ib`control'.`assignment', asobserved
		matrix p_values = r(p)
		* Interaction term will be the third column in the first row.
		local interact_test = p_values[1, 3]
	}
	
	* Run additional analysis for continuous subgroups.
	if `is_cont' == 1 {
		dis "Running mfpi analysis."
		* Run mfpi -- skipping fp2(`var') for now, it takes too long to run.
		mfpi, with(`assignment') linear(`var') fp1(`var') showmodel adjust(`covars') gendiff(_mfpi_): logit `dv' if `touse', nolog cluster(`cluster')

		* Loop through the results, assuming the linear model is the best fit.
		local model_type = "fp1 fp2"
		* Skip fp2 for now, it takes too long to run.
		local model_type = "fp1"
		
		local best_model = 1
		local best_aic = r(aiclin)
		local best_p = r(Plin)
		local current_model = 1
		
		foreach type in `model_type' {
			* Increment our model counter.
			local ++current_model
			* Determine the model with the lowest AIC.
			* Also determine if the model is significant or not.
			local aic = r(aic`type')
			local p = r(P`type')
			if `aic' < `best_aic' {
				local best_aic = `aic'
				local best_p = `p'
				local best_model = `current_model'
			}
		}
		dis "Best model: `best_model'. Best p: `best_p'. Best AIC: `best_aic'."
		* Run mfpi_plot on that model and save the result.
		local graph_name = "`expname'_subgroup_`var'_`best_model'"
		mfpi_plot `var' if `touse', vn(`best_model') saving("`graph_name'", replace) name(mfpi_curplot)
		graph export "`graph_name'.png", replace
		graph drop mfpi_curplot
	}
	
	* Run a reduced form model for comparison. Make sure that the interaction term is defined for all records.
	* logit `dv' ib`control'.`assignment' `covars' if `touse' & `var' != ., nolog or cluster(`cluster')
	* Also compute a LR test to compare the model with interacting treatment to one without.
	* This test is displayed for comparison to the contrast command that includes covariates, but is not explicitly checked in the algorithm.
	* See Long & Freese (2006). "Regression models for Categorical Indepenent Variables Using Stata". 2nd Edition. Chapter 9, p. 424.
	* lrtest . model_interact -> can't do with clustered SEs.

	* TODO: p-value cut-off should be a parameter and/or we should control for multiple comparisons.
	* TODO: ideally we would also check if the subgroup effect is clinically significant (at least 25% different than main effect).
		* (See http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2813449/ )
	local subgroup_significance = 0.1
	if (`interact_test' < `subgroup_significance') | (`is_cont' == 1 & `best_p' < `subgroup_significance') {
		dis "** Found significant interaction effect for `var' (p = " as result %06.4f round(`interact_test', .0001) as text ")."
		if `is_cont' == 1 {
			* For continuous variables analyze the results by the tiles, not the individual values.
			local subgroup_var = "etile_`var'"
			* Show distribution of the tiles.
			bysort etile_`var': su `var' if `touse', detail
		}
		else {
			local subgroup_var = "`var'"
		}
		dis "By `subgroup_var', baseline returns -- "
		bysort `subgroup_var': tab `assignment' `dv' if `touse', row chi2
		/** Skip this to keep the log concise.
		dis "By `var', basic analysis --"
		bysort `var': reg `dv' ib`control'.`assignment' if `touse', cluster(`cluster')
		*/
		dis "By `subgroup_var' with covariates (`covars') --"
		bysort `subgroup_var': reg `dv' ib`control'.`assignment' `covars' if `touse', cluster(`cluster')
	}
	else {
		dis "** Did not find significant interaction effect for `var' (p = " as result %06.4f round(`interact_test', .0001) as text ")."
	}
	dis _n
}
* TODO: print summary of interaction effects found with p-values, and list of all variables tested.

***
* Clean up variables

* Delete any mfpi variables that were created for the continuous variable analysis.
cap drop _mfpi_*
* Delete the tiles if they aren't wanted any more.
cap drop etile_*

* Report our runtime for this experimental evaluation.
timer off 83
qui timer list 83
dis _n(2) "rtcteval.ado analysis time:"
dis "Seconds: " as result r(t83) as text ". Minutes: " as result round(r(t83) / 60, 0.01) as text "."
dis "-----------------------" _n(3) 

end
