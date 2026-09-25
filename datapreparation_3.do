*-------------------------------------------------------------------------------
* Create mini panel dataset for employment, hours and income
*-------------------------------------------------------------------------------

*--------------------------------------
* Wave 1996, Study 3833
*--------------------------------------

use b960312 b960318 b960277 bcsid b960273 empstat b960273 using ///
    "/Users/user/Documents/datasets/britCohortStudy70/UKDA-3833-stata/stata/stata13/bcs96x.dta", clear

gen wave = 1996

* Rename variables
rename (b960312 b960318 b960277) (netpay_raw netpay_period hours)

* Clean negative or invalid values
replace netpay_raw     = . if netpay_raw < 0
replace netpay_period  = . if netpay_period < 0 | netpay_period >= 6
replace hours          = . if hours < 0

* Create weekly net pay variable
gen net_pay_week = netpay_raw * hours       if netpay_period == 1  // Hourly
replace net_pay_week = netpay_raw * 5       if netpay_period == 2  // Daily
replace net_pay_week = netpay_raw           if netpay_period == 3  // Weekly
replace net_pay_week = (netpay_raw * 12)/52 if netpay_period == 4  // Monthly
replace net_pay_week = netpay_raw / 52      if netpay_period == 5  // Annual

* Copy before trimming
gen net_pay_week_95 = net_pay_week

* Trim top 1% and 5%
summarize net_pay_week, detail
replace net_pay_week     = . if net_pay_week > r(p99)
replace net_pay_week_95  = . if net_pay_week > r(p95)

gen employedFT = empstat == 1 | empstat == 3
gen unemployed = empstat == 5

gen isManagerSupervisor= b960273 == 1
replace isManagerSupervisor =. if b960273 <0

keep bcsid empstat net_pay_week net_pay_week_95 hours wave employedFT unemployed isManagerSupervisor

* Harmonise employment variable

replace empstat = empstat +1 if empstat >6

* Save cleaned dataset
save "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_1996.dta", replace


*--------------------------------------
* Wave 2000, Study 5558
*--------------------------------------

use cnetpay cnetprd chours1 bcsid empstat cjsup2 vote97 votewho votenow prtysupp politint using ///
"/Users/user/Documents/datasets/britCohortStudy70/UKDA-5558-stata/stata/stata13_se/bcs2000.dta", clear

gen wave = 2000

* Rename variables
rename (cnetpay cnetprd chours1) (netpay_raw netpay_period hours)

* Clean negative or invalid values
replace netpay_raw    = . if netpay_raw < 0
replace netpay_period = . if netpay_period < 0 | netpay_period >= 6
replace hours         = . if hours < 0

* Create weekly net pay variable
gen net_pay_week = netpay_raw             if netpay_period == 1  // Hourly
replace net_pay_week = netpay_raw / 2     if netpay_period == 2  // Fortnightly
replace net_pay_week = (netpay_raw * 13) / 52 if netpay_period == 3  // Quarterly
replace net_pay_week = (netpay_raw *12)/52 if netpay_period == 4  // Monthly (12/12 = no change)
replace net_pay_week = netpay_raw / 52    if netpay_period == 5  // Annual


* Copy before trimming
gen net_pay_week_95 = net_pay_week

* Trim top 1% and 5%
summarize net_pay_week, detail
replace net_pay_week     = . if net_pay_week > r(p99)
replace net_pay_week_95  = . if net_pay_week > r(p95)


gen employedFT = empstat == 1 | empstat == 3
gen unemployed = empstat == 5

gen isManagerSupervisor = cjsup2== 1 | cjsup2==2
replace isManagerSupervisor=. if cjsup2 ==.
gen manager = cjsup2== 1 
replace manager=. if cjsup2 ==.

rename votewho votewho97
rename politint interestPolitics97

keep bcsid empstat net_pay_week net_pay_week_95 hours wave employedFT unemployed manager isManagerSupervisor vote97 votewho97 interestPolitics97

* Save cleaned dataset
save "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2000.dta", replace


*--------------------------------------
* Wave 2004, Study 5585
*--------------------------------------

use bd7ecact bcsid b7cnetpd b7cnetpy b7chour* b7otimny b7cjsup2 b7vote01 b7votewo bd7othpa b7polint using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-5585-stata/stata/stata13_se/bcs_2004_followup.dta", clear

gen wave = 2004

forvalues i=1/3{ 
	replace b7chour`i' =. if b7chour`i'<0
	}

egen hours = rowtotal (b7chour1 b7chour2 b7chour3), missing

rename (bd7ecact b7cnetpy b7cnetpd ) (empstat netpay_raw netpay_period)

replace netpay_raw    = . if netpay_raw < 0
replace netpay_period = . if netpay_period < 0 | netpay_period >= 6
replace hours         = . if hours < 0

* Create weekly net pay variable
gen net_pay_week = netpay_raw             if netpay_period == 1  // Hourly
replace net_pay_week = netpay_raw / 2     if netpay_period == 2  // Fortnightly
replace net_pay_week = (netpay_raw * 13) / 52 if netpay_period == 3  // Quarterly
replace net_pay_week = (netpay_raw*12)/52        if netpay_period == 4  // Monthly (12/12 = no change)
replace net_pay_week = netpay_raw / 52    if netpay_period == 5  // Annual


* Copy before trimming
gen net_pay_week_95 = net_pay_week

* Trim top 1% and 5%
summarize net_pay_week, detail
replace net_pay_week     = . if net_pay_week > r(p99)
replace net_pay_week_95  = . if net_pay_week > r(p95)


gen employedFT = empstat == 1 | empstat == 3
gen unemployed = empstat == 5

gen isManagerSupervisor = b7cjsup2== 1 | b7cjsup2==2
replace isManagerSupervisor=. if b7cjsup2 <0
gen manager = b7cjsup2== 1 
replace manager=. if b7cjsup2 <0

rename (b7vote01 b7votewo b7polint) (vote01 votewho01 interestPolitics01) 

keep bcsid empstat net_pay_week net_pay_week_95 hours wave employedFT unemployed isManagerSupervisor manager vote01 votewho01 interestPolitics

* Save cleaned dataset
save "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2004.dta", replace


*--------------------------------------
* Wave 2012, Study 7473 
*--------------------------------------

use  bcsid B9NETA B9NETP B9OTIMNY B9CHOUR1 B9CHOUR2 B9CHOUR3 B9CHOUR4 B9CJSUP B9SCQ6 B9SCQ7 B9SCQ4 using"/Users/user/Documents/datasets/britCohortStudy70/UKDA-7473-stata/stata/stata13/bcs70_2012_flatfile.dta", clear

gen wave = 2012

merge 1:m bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-7473-stata/stata/stata13/bcs70_2012_derived.dta", keepusing(BD9ECACT)

drop _merge



forvalues i=1/3{ 
	replace B9CHOUR`i' =. if B9CHOUR`i'<0
	}

egen hours = rowtotal (B9CHOUR1 B9CHOUR2 B9CHOUR3), missing

rename (BD9ECACT B9NETA B9NETP B9SCQ4 ) (empstat netpay_raw netpay_period interestPolitics12)

replace netpay_raw    = . if netpay_raw < 0
replace netpay_period = . if netpay_period < 0 | netpay_period >= 6
replace hours         = . if hours < 0

* Create weekly net pay variable
gen net_pay_week = netpay_raw             if netpay_period == 1  // Hourly
replace net_pay_week = netpay_raw / 2     if netpay_period == 2  // Fortnightly
replace net_pay_week = (netpay_raw * 13) / 52 if netpay_period == 4  // Quarterly
replace net_pay_week = (netpay_raw*12)/52        if netpay_period == 5  // Monthly (12/12 = no change)
replace net_pay_week = netpay_raw / 3    if netpay_period == 3  // THree weeks


* Copy before trimming
gen net_pay_week_95 = net_pay_week

* Trim top 1% and 5%
summarize net_pay_week, detail
replace net_pay_week     = . if net_pay_week > r(p99)
replace net_pay_week_95  = . if net_pay_week > r(p95)


gen employedFT = empstat == 1 | empstat == 3
gen unemployed = empstat == 5


gen isManagerSupervisor = B9CJSUP== 1 | B9CJSUP==2
replace isManagerSupervisor=. if B9CJSUP <0
gen manager = B9CJSUP== 1 
replace manager=. if B9CJSUP <0

gen vote10 = B9SCQ6 >0 & B9SCQ6 !=16
rename B9SCQ6 votewho10
replace votewho10 = . if votewho10<0 | votewho10>=16

gen vote05 = B9SCQ7 >0 & B9SCQ7 !=16
rename B9SCQ7 votewho05
replace votewho05 = . if votewho05<0 | votewho05>=16


keep bcsid empstat net_pay_week net_pay_week_95 hours wave employedFT unemployed manager isManagerSupervisor vote* votewho* interestPolitics12

* Save cleaned dataset
save "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2012.dta", replace


*--------------------------------------
* Wave 2016, Study 8547 
*--------------------------------------


use B10CHOUR1 B10NETW BD10ECACT B10CJSUPA B10VOTE01 B10VOTEWO1 B10VOTE02 B10VOTEWO2 bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-8547-stata/stata/stata13/bcs_age46_main.dta", clear

gen wave = 2016 

rename (B10CHOUR1 B10NETW BD10ECACT )(hours net_pay_week empstat)

replace net_pay_week =. if net_pay_week <0
* Copy before trimming
gen net_pay_week_95 = net_pay_week

* Trim top 1% and 5%
summarize net_pay_week, detail
replace net_pay_week     = . if net_pay_week > r(p99)
replace net_pay_week_95  = . if net_pay_week > r(p95)


gen employedFT = empstat == 1 | empstat == 3
gen unemployed = empstat == 5

gen isManagerSupervisor = B10CJSUPA== 1 
replace isManagerSupervisor=. if B10CJSUPA <0

rename(B10VOTE01 B10VOTEWO1 B10VOTE02 B10VOTEWO2)(vote15 votewho15 vote17 votewho17)

keep bcsid empstat net_pay_week net_pay_week_95 hours wave employedFT unemployed isManagerSupervisor vote* votewho* 
* Save cleaned dataset
save "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2016.dta", replace



*--------------------------------------
* Wave 2021, Study 9347 
*--------------------------------------

use b11econact2 b11netw b11chours* bcsid b11cjsupa b11vote04 b11votewho4 b11vote02 b11votewho2 b11voteeu b11voteeuway b11polint b11cflisn b11cflisd b11cfani b11cfcor b11cfmis b11cftot b11cfrc b11q25a b11q25b b11q25c b11q25d b11q25e b11q25f b11q25g b11q25h b11q25i b11q25j b11q25k b11q25l b11q25m b11q25n b11q25o b11q25p b11q25q b11q25r b11q25s b11q25t using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-9347-stata/stata/stata13/bcs11_age51_main.dta", clear

gen wave = 2021


foreach i of numlist 1 3 4 5{ 
	replace b11chours`i' =. if b11chours`i'<0
	}

egen hours = rowtotal (b11chours1 b11chours3), missing

rename (b11econact2 b11netw ) (empstat net_pay_week)
replace net_pay_week =. if net_pay_week <0
* Copy before trimming
gen net_pay_week_95 = net_pay_week

* Trim top 1% and 5%
summarize net_pay_week, detail
replace net_pay_week     = . if net_pay_week > r(p99)
replace net_pay_week_95  = . if net_pay_week > r(p95)


gen employedFT = empstat == 1 | empstat == 3
gen unemployed = empstat == 5

gen isManagerSupervisor = b11cjsupa== 1 
replace isManagerSupervisor=. if b11cjsupa <0

rename(b11vote04 b11votewho4 b11voteeu b11voteeuway b11polint)(vote19 votewho19 voteRef voteRefWhat interestPolitics21)

keep bcsid empstat net_pay_week net_pay_week_95 hours wave employedFT unemployed isManagerSupervisor vote*  interestPolitics  b11cflisn b11cflisd b11cfani b11cfcor b11cfmis b11cftot b11cfrc b11q25a b11q25b b11q25c b11q25d b11q25e b11q25f b11q25g b11q25h b11q25i b11q25j b11q25k b11q25l b11q25m b11q25n b11q25o b11q25p b11q25q b11q25r b11q25s b11q25t
* Save cleaned dataset
save "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2021.dta", replace


append using "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_1996.dta""/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2000.dta""/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2004.dta""/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2012.dta" "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2016.dta""/Users/user/Documents/datasets/britCohortStudy70/temporary_files/bcs_2021.dta"


save "/Users/user/Documents/datasets/britCohortStudy70/temporary_files/panel_labour_outcomes.dta", replace
