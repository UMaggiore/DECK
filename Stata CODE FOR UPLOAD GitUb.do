
********************************************************************************
**# START TABLE 1
********************************************************************************
clear
cd "C:\Documenti\Furian\DECK"
use deck_trasv, replace

keep if pz_categoria == 2 | pz_categoria == 4

dtable i.volume_attivitàVIV_centrotx i.volume_attivitàDEC_centrotx ///
     i.don_effettuata i.don_tecnica_nefrectomiaVIV don_eta i.don_sesso ///
	 don_gfr don_quality_score ///
	don_altezza don_peso don_BMI i.don_razza i.don_ipertensione ///
	i.don_diabete i.don_causadecesso i.DON_GRUPPO ///
	i.don_rischio ///
	i.causa_ingressoKPD ric_eta i.ric_sesso i.RIC_GRUPPO ///
	ric_PRA ric_CIT ric_peso i.ric_dialisi ric_mesidialisi i.ric_DGF ///
	i.ric_rigetto i.ric_complicanzatx i.ric_complicanzachir i.ric_graftloss ///
	i.causa_graftloss i.ric_decesso ///
	, ///	
	by(pz_categoria, tests) ///	
	define(meansd = mean sd, delimiter(" ± ")) ///
	define(myiqr = p25 p75, delimiter("-")) ///
	define(myrange = min max, delimiter("-")) ///
	factor(volume_attivitàVIV_centrotx volume_attivitàDEC_centrotx ///
	don_effettuata don_tecnica_nefrectomiaVIV  don_sesso ///
	i.don_razza don_ipertensione ///
	don_diabete don_causadecesso DON_GRUPPO ///
	don_rischio ///
	causa_ingressoKPD  ric_sesso RIC_GRUPPO  ///
	ric_dialisi  ric_DGF ric_rigetto ric_complicanzatx ric_complicanzachir ///
	ric_graftloss causa_graftloss ric_decesso, test(fisher)) ///
    continuous(ric_eta ric_PRA ric_CIT ric_peso ric_mesidialisi ///
	don_eta don_gfr don_altezza don_peso don_BMI  ///
	, stat(meansd) test(kwallis)) ///
   continuous(don_quality_score ric_mesidialisi ///
	, stat(median myrange) test(kwallis)) /// 
	column(by(hide) test(p-value)) ///
	title(Table 1. Recipient Characheristics) ///
	note(Mann-Whitney test for continuous variables ///
	(reported as mean ± standard deviation or median (min - max)).) ///
    note(Fisher's exact test for categorical variables ///
	(reported as number (percentage)).) ///
	note(Baseline characteristics of the study population) ///
	note(BMI, Body Mass Index; CNT, Centro nazionale Trapianti; DEC-K, ///
	DECeased-Kidney-paired-exchange; CEK, Compatible End-Kidney; CIK, ///
	Compatible Initial Kidney; HD, HemoDialysis; ///
	KDPI, Kidney Donor Profile Index; LD, living donor; ///
	LDKT, living donor kidney transplantation; ///
	LKPI, Living Kidney Donor Profile Index; PD, Peritoneal Dialysis; ///
	KRT, kidney replacement therapy; Tx, transplant; VLS, videolaparoscopy) ///
	sformat("%s" sd) ///
	nformat("%3.1f" mean sd median p25 p75) ///
	nformat("%3.1f" min max) ///
	sformat("(%s)" myiqr myrange) ///
    nformat("%3.0f" N count fvfrequency) ///
    nformat("%3.1f" fvpercent ) ///
    nformat("%6.3f" kwallis fisher) ///
	export(table1.html, replace)
collect export table 1.xlsx, replace
collect export table 1.docx, replace
collect export table 1.txt, replace
collect export table 1.html, replace


*******************************************************************************
**# END TABLE 1
*******************************************************************************


********************************************************************************
**# START LONGITUDINAL ANALYSES
********************************************************************************
**# restricted cubic spline crude random coefficient with trajectories (linear part for random slopes)

  clear
cd "C:\Documenti\Furian\DECK"
use deck_long, replace

keep if pz_categoria == 2 | pz_categoria == 4


* Optional: sort data for plotting
sort id mese
* Set up the plot command
local plotcmd ""
* Loop over each individual
levelsof id, local(ids)
foreach i of local ids {
    * Get the category for this individual
    quietly summarize pz_categoria if id == `i'
    local cat = r(mean)
    * Choose color based on category
    local color = cond(`cat' == 4, "stc2", "stc1")
    * Add line for this individual to the plot command
    local plotcmd `plotcmd' || line  eGFR_EPI2021 mese if id == `i', lcolor(`color'%30) lwidth(thin)
}


unique mese
global distinct = r(unique)
makespline rcs mese,  distinct($distinct) 
matrix list  r(knots)
  

mixed eGFR_EPI2021 ib4.pz_categoria##(c._rs_rcs_1 c._rcs_1_*)  ///
|| id: _rs_rcs_1, cov(unstr) reml dfm(kroger)
estimates stat, aicconsistent
lincom _b[2.pz_categoria#c._rs_rcs_1] / 6, small sformat(%3.2f) ///
pformat(%3.2f) cformat(%3.1f)  // 1 is 72 months (6 years)
// different linear + spline trend
qui test _b[2.pz_categoria#c._rs_rcs_1] = 0, small
test _b[2.pz_categoria#c._rcs_1_1]=0, small accum
local p_int = r(p)
local spval_int = string( `p_int', "%3.2f")
di `spval_int'
contrast pz_categoria, small
local p_avg = el(r(p),1,1)
local spval_avg = string( `p_avg', "%3.2f")
di `spval_avg'

local atlist
mkmat _rs_rcs_1 _rcs_1_1 in 1/$distinct,matrix(A)
forvalues i = 1/8 {
        local x1 = A[`i',1]
        local x2 = A[`i',2]
        local atlist `atlist' at(_rs_rcs_1= `x1' _rcs_1_1=`x2')
}
margins, `atlist' over(pz_categoria) saving(marg_egfr_sc, replace)
append using marg_egfr_sc
recode _at (1 = 1) (2 = 3) (3 = 6) (4 = 12) (5 = 24) (6 = 36) (7 = 48) (8 = 60)

sort id _rs_rcs_ 
tw line _margin _at if  _by1 == 2, lcolor(stc1)  || ///
   line _margin _at if  _by1 == 4, lcolor(stc2) ||  ///
   scatter _margin _at if  _by1 == 2, msymbol(i)  mcolor(stc1) || ///
   scatter _margin _at if  _by1 == 4, msymbol(i)  mcolor(stc2) ||  ///
   rarea _ci_ub _ci_lb _at if  _by1 == 2, color(stc1%20) lcolor(white)  || ///
   rarea _ci_ub _ci_lb _at if  _by1 == 4, color(stc2%20) lcolor(white)  ||  ///
  `plotcmd' || ///
    ,  ///
	xtitle("Month since Transplant")  xlab(0 `" " "Tx"' 1 3 6 12 24 36 48 60, ///
	labsize(*1)) xsc(range (-1 60))  ///
    ytitle("eGFR ml/min/1.73m{sup:2}") ylabel(0 15 30 45 60 90 120) ///
    legend(order(1 "CIK recipient" 2 "LD recipient")) ///
	title("") ///
	text(3 36 "P value for different slopes=`spval_int'") ///
	note("Non-Adjusted analysis" "") ///
	name(individual_segfr, replace) 
	
graph export egfr_scont_it_crude.tif, replace
graph export egfr_scont_it_crude.pdf, replace
graph export egfr_scont_it_crude.png, replace
graph export egfr_scont_it_crude.svg, replace


* Conversion to TIFF 900 dpi via integrated Python
python:
from PIL import Image

# Path PNG exported from Stata
input_png = "egfr_scont_it_crude.png"
output_tiff = "FIG2A_900dpi.tif"

# Opern image and convert to RGB
img = Image.open(input_png).convert("RGB")

# Save as TIFF with 900 dpi
img.save(output_tiff, dpi=(900, 900), compression="tiff_deflate")

print("Completed coversion: egfr_scont_it_crude.tif with 900 dpi")
end



**# restricted cubic spline adjusted random coefficient with trajectories (linear part for the random slopes)
  clear
cd "C:\Documenti\Furian\DECK"
use deck_long, replace

keep if pz_categoria == 2 | pz_categoria == 4


* Optional: sort data for plotting
sort id mese
* Set up the plot command
local plotcmd ""
* Loop over each individual
levelsof id, local(ids)
foreach i of local ids {
    * Get the category for this individual
    quietly summarize pz_categoria if id == `i'
    local cat = r(mean)
    * Choose color based on category
    local color = cond(`cat' == 4, "stc2", "stc1")
    * Add line for this individual to the plot command
    local plotcmd `plotcmd' || line  eGFR_EPI2021 mese if id == `i', ///
	lcolor(`color'%30) lwidth(thin)
}


unique mese
global distinct = r(unique)
makespline rcs mese,  distinct($distinct) 
matrix list  r(knots)

foreach var of varlist don_quality_score don_eta  ric_CIT  ric_eta ric_PRA  {
	qui summ `var'
	replace  `var' = (`var' - r(mean)) / r(sd)
	 }

mixed eGFR_EPI2021 ib4.pz_categoria##(c._rs_rcs_1 c._rcs_1_*) ///
	c.don_quality_score c.don_eta c.ric_CIT i.volume_attivitàDEC_centrotx ///
	c.ric_eta i.ric_sesso c.ric_PRA || id: _rs_rcs_1, cov(unstr) reml dfm(kroger)
estimates stat, aicconsistent 
lincom _b[2.pz_categoria#c._rs_rcs_1] / 6, small sformat(%3.2f) pformat(%3.2f) ///
	cformat(%3.1f) // 1 is 72 months (6 years)
// different linear + spline trend
qui test _b[2.pz_categoria#c._rs_rcs_1] = 0, small
test _b[2.pz_categoria#c._rcs_1_1]=0, small accum
local p_int = r(p)
local spval_int = string( `p_int', "%3.2f")
di `spval_int'
contrast pz_categoria, small
local p_avg = el(r(p),1,1)
local spval_avg = string( `p_avg', "%3.2f")
di `spval_avg'
	
local atlist
mkmat _rs_rcs_1 _rcs_1_1 in 1/$distinct, matrix(A)
forvalues i = 1/8 {
        local x1 = A[`i',1]
        local x2 = A[`i',2]
        local atlist `atlist' at(_rs_rcs_1= `x1' _rcs_1_1=`x2')
}
margins, `atlist' over(pz_categoria) saving(marg_egfr_sc, replace)
append using marg_egfr_asc
recode _at (1 = 1) (2 = 3) (3 = 6) (4 = 12) (5 = 24) (6 = 36) (7 = 48) (8 = 60)

sort id _rs_rcs_ 
tw line _margin _at if  _by1 == 2, lcolor(stc1)  || ///
   line _margin _at if  _by1 == 4, lcolor(stc2) ||  ///
   scatter _margin _at if  _by1 == 2, msymbol(i)  mcolor(stc1) || ///
   scatter _margin _at if  _by1 == 4, msymbol(i)  mcolor(stc2) ||  ///
   rarea _ci_ub _ci_lb _at if  _by1 == 2, color(stc1%20) lcolor(white)  || ///
   rarea _ci_ub _ci_lb _at if  _by1 == 4, color(stc2%20) lcolor(white)  ||  ///
  `plotcmd' || ///
    ,  ///
	xtitle("Month since Transplant")  xlab(0 `" " "Tx"' 1 3 6 12 24 36 48 60, ///
	labsize(*1)) xsc(range (-1 60))  ///
    ytitle("eGFR ml/min/1.73m{sup:2}") ylabel(0 15 30 45 60 90 120) ///
    legend(order(1 "CIK recipient" 2 "LD recipient")) ///
	title("") ///
	text(3 40 "P value for different slopes=`spval_int'") ///
	note("Adjusted for donor quality (LDKPI/KDPI), Donor age, Cold ischemia time, Center volume," "Recipient age, sex and cPRA") /// 
	name(individual_asegfr, replace) 
	
graph export egfr_scont_it_adj.tif, replace
graph export egfr_scont_it_adj.pdf, replace
graph export egfr_scont_it_adj.png, replace
graph export egfr_scont_it_adj.svg, replace

* Conversion to TIFF 900 dpi with integrated Python
python:
from PIL import Image

# Path of PNG exported from Stata
input_png = "egfr_scont_it_adj.png"
output_tiff = "FIG2B_900dpi.tif"

# Open image and convert to RGB
img = Image.open(input_png).convert("RGB")

# Save as TIFF with 900 dpi
img.save(output_tiff, dpi=(900, 900), compression="tiff_deflate")

print("Completed coversion: egfr_scont_it_adj.tif with 900 dpi")
end



 **# linear crude random coefficient with trajectories
  clear
cd "C:\Documenti\Furian\DECK"
use deck_long, replace

keep if pz_categoria == 2 | pz_categoria == 4


* Optional: sort data for plotting
sort id mese
* Set up the plot command
local plotcmd ""
* Loop over each individual
levelsof id, local(ids)
foreach i of local ids {
    * Get the category for this individual
    quietly summarize pz_categoria if id == `i'
    local cat = r(mean)
    * Choose color based on category
    local color = cond(`cat' == 4, "stc2", "stc1")
    * Add line for this individual to the plot command
    local plotcmd `plotcmd' || line  eGFR_EPI2021 mese if id == `i', ///
	lcolor(`color'%30) lwidth(thin)
}


mixed eGFR_EPI2021 ib4.pz_categoria##(c.mese)  || id: c.mese, cov(unstr) reml dfm(kroger)
lincom _b[2.pz_categoria#c.mese] * 12, small sformat(%3.2f) pformat(%3.2f) cformat(%3.1f) 
// different linear trend
test _b[2.pz_categoria#c.mese] = 0, small
local p_int = r(p)
local spval_int = string( `p_int', "%3.2f")
di `spval_int'
contrast pz_categoria, small
local p_avg = el(r(p),1,1)
local spval_avg = string( `p_avg', "%3.2f")
di `spval_avg'

margins, at(mese=(1 3 6 12 24 36 48 60)) over(pz_categoria) saving(marg_egfr_c, replace)
append using marg_egfr_c
recode _at (1 = 1) (2 = 3) (3 = 6) (4 = 12) (5 = 24) (6 = 36) (7 = 48) (8 = 60)

sort id mese
tw line _margin _at if  _by1 == 2, lcolor(stc1)  || ///
   line _margin _at if  _by1 == 4, lcolor(stc2) ||  ///
   scatter _margin _at if  _by1 == 2, msymbol(i)  mcolor(stc1) || ///
   scatter _margin _at if  _by1 == 4, msymbol(i)  mcolor(stc2) ||  ///
   rarea _ci_ub _ci_lb _at if  _by1 == 2, color(stc1%20) lcolor(white)  || ///
   rarea _ci_ub _ci_lb _at if  _by1 == 4, color(stc2%20) lcolor(white)  ||  ///
   `plotcmd' || ///
    ,  ///
	xtitle("Month since Transplant")  xlab(0 `" " "Tx"' 1 3 6 12 24 36 48 60, labsize(*1)) xsc(range (-1 60))  ///
    ytitle("eGFR ml/min/1.73m{sup:2}") ylabel(0 15 30 45 60 90 120) ///
    legend(order(1 "CIK recipient" 2 "LD recipient")) ///
	title("") ///
	text(3 40 "P value for different slopes=`spval_int'") ///
	note("Non-Adjusted analysis" "") ///
	name(individual_egfr, replace) 
	
graph export egfr_cont_it_crude.tif, replace
graph export egfr_cont_it_crude.pdf, replace
graph export egfr_cont_it_crude.png, replace
graph export egfr_cont_it_crude.svg, replace


**# linear adjusted random coefficient with trajectories
 
clear
cd "C:\Documenti\Furian\DECK"
use deck_long, replace

keep if pz_categoria == 2 | pz_categoria == 4

* Optional: sort data for plotting
sort id mese
* Set up the plot command
local plotcmd ""
* Loop over each individual
levelsof id, local(ids)
foreach i of local ids {
    * Get the category for this individual
    quietly summarize pz_categoria if id == `i'
    local cat = r(mean)
    * Choose color based on category
    local color = cond(`cat' == 4, "stc2", "stc1")
    * Add line for this individual to the plot command
    local plotcmd `plotcmd' || line  eGFR_EPI2021 mese if id == `i', ///
	lcolor(`color'%30) lwidth(thin)
}



foreach var of varlist don_quality_score don_eta  ric_CIT  ric_eta ric_PRA  {
	qui summ `var'
	replace  `var' = (`var' - r(mean)) / r(sd)
	 }

mixed eGFR_EPI2021 ib4.pz_categoria##(c.mese) ///
	c.don_quality_score c.don_eta c.ric_CIT i.volume_attivitàDEC_centrotx ///
	c.ric_eta i.ric_sesso c.ric_PRA  || id: c.mese, cov(unstr) reml dfm(kroger)
lincom _b[2.pz_categoria#c.mese] * 12, small sformat(%3.2f) pformat(%3.2f) cformat(%3.1f) 
// different linear + spline trend
test _b[2.pz_categoria#c.mese] = 0, small
local p_int = r(p)
local spval_int = string( `p_int', "%3.2f")
di `spval_int'
contrast pz_categoria, small
local p_avg = el(r(p),1,1)
local spval_avg = string( `p_avg', "%3.2f")
di `spval_avg'



margins, at(mese=(1 3 6 12 24 36 48 60)) over(pz_categoria) saving(marg_egfr_c, replace)
append using marg_egfr_c
recode _at (1 = 1) (2 = 3) (3 = 6) (4 = 12) (5 = 24) (6 = 36) (7 = 48) (8 = 60)

sort id mese
tw line _margin _at if  _by1 == 2, lcolor(stc1)  || ///
   line _margin _at if  _by1 == 4, lcolor(stc2) ||  ///
   scatter _margin _at if  _by1 == 2, msymbol(i)  mcolor(stc1) || ///
   scatter _margin _at if  _by1 == 4, msymbol(i)  mcolor(stc2) ||  ///
   rarea _ci_ub _ci_lb _at if  _by1 == 2, color(stc1%20) lcolor(white)  || ///
   rarea _ci_ub _ci_lb _at if  _by1 == 4, color(stc2%20) lcolor(white)  ||  ///
  `plotcmd' || ///
    ,  ///
	xtitle("Month since Transplant")  xlab(0 `" " "Tx"' 1 3 6 12 24 36 48 60, ///
	labsize(*1)) xsc(range (-1 60))  ///
    ytitle("eGFR ml/min/1.73m{sup:2}") ylabel(0 15 30 45 60 90 120) ///
    legend(order(1 "CIK recipient" 2 "LD recipient")) ///
	title("") ///
	text(3 40 "P value for different slopes=`spval_int'") ///
	note("Non-Adjusted analysis" "") ///
	name(individual_egfr, replace) 
	

graph export egfr_cont_it_adj.tif, replace
graph export egfr_cont_it_adj.pdf, replace
graph export egfr_cont_it_adj.png, replace
graph export egfr_cont_it_adj.svg, replace


********************************************************************************
**# END LONGITUDINAL ANALYSES
********************************************************************************




********************************************************************************
**# START SURVIVAL ANALYSES
********************************************************************************

**# Survival curve crude

 clear
cd "C:\Documenti\Furian\DECK"
use deck_trasv, replace

keep if pz_categoria == 2 | pz_categoria == 4

sts test pz_categoria, logrank
local p_logrank = chi2tail(r(df), r(chi2))
local spval_logrank = string( `p_logrank', "%3.2f")
di `spval_logrank'
stcox ib4.pz_categoria
estat phtest
lincom _b[2.pz_categoria], hr cformat(%3.2f) pformat(%3.2f) sformat(%3.2f)

streset, sc(30.4375)
sts graph, by(pz_categoria) risktable ///
 xtitle("Month since Transplant") xsc(titlegap(1)) ///
 xlab(0 `" " "Tx"' 1 3 6 12 24 36 48 60, labsize(*1)) tmax(60) ///
 ytitle("Transplant Survival (%)") ///
 ylab(0 "0" .2 "20" .4 "40" .6 "60" .8 "80" 1 "100") ///
 ci ///
 plot1opts(lwidth(*1.4)) plot2opts(lwidth(*1.4)) ///
 ci1opts(color(stc1%20))  ci2opts(color(stc2%20)) ///
 risktable(, title("N at risk", size(*.8))) ///
 risktable(, color(stc1) group(#1) size(*.7)) ///
 risktable(, color(stc2)  group(#2) size(*.7))  ///
 risktable(, rowtitle("CIK Rec:  ") group(#1) size(*.8)) ///
 risktable(, rowtitle("LD  Rec:  ") group(#2) size(*.8)) /// 
 legend(order(5 "CIK recipient" 6 "LD Recipient") pos(6) rows(1)) ///
 text(.05 40 "Log-rank test,  P value=`spval_logrank'" ) ///
 title("")
 
graph export km_crude.tif, replace
graph export km_crude.pdf, replace
graph export km_crude.png, replace
graph export km_crude.svg, replace





**# Adjusted HR

clear
cd "C:\Documenti\Furian\DECK"
use deck_trasv, replace

keep if pz_categoria == 2 | pz_categoria == 4

foreach var of varlist don_quality_score don_eta  ric_CIT  ric_eta ric_PRA  {
	qui summ `var'
	replace  `var' = (`var' - r(mean)) / r(sd)
	 }

stcox ib4.pz_categoria  c.don_quality_score c.ric_CIT 
lincom  _b[2.pz_categoria], hr cformat(%3.2f) pformat(%3.2f) sformat(%3.2f)


********************************************************************************
**# END SURVIVAL ANALYSES
********************************************************************************


 

