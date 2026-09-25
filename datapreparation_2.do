
* ------------------------------------------------------------------------------
* Cognitive Tests
*
*	There are several cognitive tests in the study, however only three have 
*	distribution dates: the Edingburgh, Fun math and BAS tests (the latter has
*	a number of sub-sections which tap at differnt aspects of cognition).
*	In addition to these, we find a "vocabulary test (i8-i81)", a dictation
*	exercise (i3864-3815) and a reading exercise i111-i216. Becuase the date of 
*	the tests varied, and because development between ages 10 and 11 is not
*	trivial, the latter tests are ignored -as we are going to age standardise
*	the scores
* ------------------------------------------------------------------------------



* ------------------------------------------------------------------------------
* Test FUN math: We use the derived variables, i4001 to i4072
*
* We score "NS/No response (-3)" not as missing value, but as incorrect answer. 
* Same 
* with "More than 1 answer (-4)". Only -6 (No questinnaire) is taken as proper 
* missing
* ------------------------------------------------------------------------------





local testvars i4001 i4002 i4003 i4004 i4005 i4006 i4007 i4008 i4009 i4010 i4011 i4012 i4013 i4014 i4015 i4016 i4017 i4018 i4019 i4020 i4021 i4022 i4023 i4024 i4025 i4026 i4027 i4028 i4029 i4030 i4031 i4032 i4033 i4034 i4035 i4036 i4037 i4038 i4039 i4040 i4041 i4042 i4043 i4044 i4045 i4046 i4047 i4048 i4049 i4050 i4051 i4052 i4053 i4054 i4055 i4056 i4057 i4058 i4059 i4060 i4061 i4062 i4063 i4064 i4065 i4066 i4067 i4068 i4069 i4070 i4071 i4072

foreach var of varlist `testvars' {
       replace `var' = . if inlist(`var',-6)
	replace `var' = 2 if `var' ==-3 | `var' == -4
	replace `var' = . if `var' <0
}

foreach var of varlist `testvars' {
    gen `var'_bin = `var' == 1
	replace `var'_bin  =. if `var' ==.
}



egen fmt_total = rowtotal(i4001_bin i4002_bin i4003_bin i4004_bin i4005_bin i4006_bin i4007_bin i4008_bin i4009_bin i4010_bin i4011_bin i4012_bin i4013_bin i4014_bin i4015_bin i4016_bin i4017_bin i4018_bin i4019_bin i4020_bin i4021_bin i4022_bin i4023_bin i4024_bin i4025_bin i4026_bin i4027_bin i4028_bin i4029_bin i4030_bin i4031_bin i4032_bin i4033_bin i4034_bin i4035_bin i4036_bin i4037_bin i4038_bin i4039_bin i4040_bin i4041_bin i4042_bin i4043_bin i4044_bin i4045_bin i4046_bin i4047_bin i4048_bin i4049_bin i4050_bin i4051_bin i4052_bin i4053_bin i4054_bin i4055_bin i4056_bin i4057_bin i4058_bin i4059_bin i4060_bin i4061_bin i4062_bin i4063_bin i4064_bin i4065_bin i4066_bin i4067_bin i4068_bin i4069_bin i4070_bin i4071_bin i4072_bin), missing

drop *_bin





* Next, we compute BAS as four separate components
* BAS1 Word definition (i3504-i3540)
* Here, we also set "No response" as incorrect/unacceptable.

local testvars i3504 i3505 i3506 i3507 i3508 i3509 i3510 i3511 i3512 i3513 i3514 i3515 i3516 i3517 i3518 i3519 i3520 i3521 i3522 i3523 i3524 i3525 i3526 i3527 i3528 i3529 i3530 i3531 i3532 i3533 i3534 i3535 i3536 i3537 i3538 i3539 i3540

foreach var of varlist `testvars' {
    replace `var' = . if inlist(`var',-6)
	replace `var' = 2 if `var' ==-3 | `var' == -4 | `var' ==9
	replace `var' = . if `var' <0
}

foreach var of varlist `testvars' {
    gen `var'_bin = `var' == 1
	replace `var'_bin  =. if `var' ==.
}

egen bas_def_total = rowtotal(i3504_bin i3505_bin i3506_bin i3507_bin i3508_bin i3509_bin i3510_bin i3511_bin i3512_bin i3513_bin i3514_bin i3515_bin i3516_bin i3517_bin i3518_bin i3519_bin i3520_bin i3521_bin i3522_bin i3523_bin i3524_bin i3525_bin i3526_bin i3527_bin i3528_bin i3529_bin i3530_bin i3531_bin i3532_bin i3533_bin i3534_bin i3535_bin i3536_bin i3537_bin i3538_bin i3539_bin i3540_bin), missing
drop *_bin


* BAS2, Number recall (i3541 - i3574)
local testvars i3541 i3542 i3543 i3544 i3545 i3546 i3547 i3548 i3549 i3550 i3551 i3552 i3553 i3554 i3555 i3556 i3557 i3558 i3559 i3560 i3561 i3562 i3563 i3564 i3565 i3566 i3567 i3568 i3569 i3570 i3571 i3572 i3573 i3574

foreach var of varlist `testvars' {
    replace `var' = . if inlist(`var',-6)
	replace `var' = 2 if `var' ==-3 | `var' == -4 | `var' ==9
	replace `var' = . if `var' <0
}

foreach var of varlist `testvars' {
    gen `var'_bin = `var' == 1
	replace `var'_bin  =. if `var' ==.
}

egen bas_numrec_total =rowtotal(i3541_bin i3542_bin i3543_bin i3544_bin i3545_bin i3546_bin i3547_bin i3548_bin i3549_bin i3550_bin i3551_bin i3552_bin i3553_bin i3554_bin i3555_bin i3556_bin i3557_bin i3558_bin i3559_bin i3560_bin i3561_bin i3562_bin i3563_bin i3564_bin i3565_bin i3566_bin i3567_bin i3568_bin i3569_bin i3570_bin i3571_bin i3572_bin i3573_bin i3574_bin)
drop *_bin


*BAS3, Similarities (i3575-i3616)
local testvars i3575 i3576 i3577 i3578 i3579 i3580 i3581 i3582 i3583 i3584 i3585 i3586 i3587 i3588 i3589 i3590 i3591 i3592 i3593 i3594 i3595 i3596 i3597 i3598 i3599 i3600 i3601 i3602 i3603 i3604 i3605 i3606 i3607 i3608 i3609 i3610 i3611 i3612 i3613 i3614 i3615 i3616


foreach var of varlist `testvars' {
    replace `var' = . if inlist(`var',-6)
	replace `var' = 2 if `var' ==-3 | `var' == -4 | `var' ==9
	replace `var' = . if `var' <0
}

foreach var of varlist `testvars' {
    gen `var'_bin = `var' == 1
	replace `var'_bin  =. if `var' ==.
}

egen bas_sims_total =rowtotal(i3575_bin i3576_bin i3577_bin i3578_bin i3579_bin i3580_bin i3581_bin i3582_bin i3583_bin i3584_bin i3585_bin i3586_bin i3587_bin i3588_bin i3589_bin i3590_bin i3591_bin i3592_bin i3593_bin i3594_bin i3595_bin i3596_bin i3597_bin i3598_bin i3599_bin i3600_bin i3601_bin i3602_bin i3603_bin i3604_bin i3605_bin i3606_bin i3607_bin i3608_bin i3609_bin i3610_bin i3611_bin i3612_bin i3613_bin i3614_bin i3615_bin i3616_bin), missing
drop *_bin


* BAS4, Matrices
local testvars i3617 i3618 i3619 i3620 i3621 i3622 i3623 i3624 i3625 i3626 i3627 i3628 i3629 i3630 i3631 i3632 i3633 i3634 i3635 i3636 i3637 i3638 i3639 i3640 i3641 i3642 i3643 i3644

foreach var of varlist `testvars' {
    replace `var' = . if inlist(`var',-6)
	replace `var' = 2 if `var' ==-3 | `var' == -4 | `var' ==9
	replace `var' = . if `var' <0
}

foreach var of varlist `testvars' {
    gen `var'_bin = `var' == 1
	replace `var'_bin  =. if `var' ==.
}

egen bas_matrix_total =rowtotal(i3617_bin i3618_bin i3619_bin i3620_bin i3621_bin i3622_bin i3623_bin i3624_bin i3625_bin i3626_bin i3627_bin i3628_bin i3629_bin i3630_bin i3631_bin i3632_bin i3633_bin i3634_bin i3635_bin i3636_bin i3637_bin i3638_bin i3639_bin i3640_bin i3641_bin i3642_bin i3643_bin i3644_bin), missing
drop *_bin

* Next, we age standardise the variables.


replace dobc10 = 1970  // 1970 instead of 70
gen dob =ym(dobc10, dobb10)

* Next, test Fun math
replace i2503y = 1980 if i2503y== 80
replace i2503y = 1981 if i2503y== 81
gen fmt_date = ym(i2503y, i2503m)
gen age_at_fmt = (fmt_date - dob)

* Age at BAS
replace i3503y = 1980 if i3503y== 80
replace i3503y = 1981 if i3503y== 81
gen bas_date = ym(i3503y, i3503m)
gen age_at_bas = (bas_date - dob)

/* 
Most kids were tested in the BAS and FM at the same age, even at the same 
day, however seveal kids were tested with several months difference as can be seen
in these tables which tabulate age (in years) when taking the BAS and hte FMT

           |             check_b
   check_f |         9         10         11 |     Total
-----------+---------------------------------+----------
         9 |     1,557        163          0 |     1,720 
        10 |       239      8,117         19 |     8,375 
        11 |         1         16        145 |       162 
-----------+---------------------------------+----------
     Total |     1,797      8,296        164 |    10,257 


The implication is that we cannot/should use any of the other measures for which 
the date of the test is not explicit (the Edinburgh in particular). So we focus 
on the BAS and FMT (which still give us a rich variety of measures).


For standardisation, we are going to use periods of two months. Ages above 67 
bimesters are  included as 67 (these are 3 and 4 observations in the FMT and BAS 
respectively). Similarly, 58 bimesters are included in the 59 bimesters; this affects 12 
and 14 observations only.
*/

gen bimester_f = floor(age_at_f/2)
replace bimester_f = 59 if bimester_f == 58
replace bimester_f = 67 if bimester_f == 68 | bimester_f == 69

foreach var in fmt_total  {
	egen `var'_sd = std(`var'), by(bimester_f)

	
}



gen bimester_b = floor(age_at_b/2)
replace bimester_b = 59 if bimester_b == 58
replace bimester_b = 67 if bimester_b == 68 | bimester_b == 69

foreach var in bas_def_total bas_matrix_total bas_numrec_total bas_sims_total {
	egen `var'_sd = std(`var'), by(bimester_b)
	
}

pca *_sd ,  mineigen(1)

predict pc1


* Next, perserverance score. We age standardise this measure as well.
replace j009c = 1980 if j009c== 80
replace j009c = 1981 if j009c== 81
gen edq_date = ym(j009c, j009b)
gen age_at_edq = (edq_date - dob)

gen bimester_edq = floor(age_at_edq/2)
replace bimester_edq = 59 if bimester_edq == 58
replace bimester_edq = 67 if bimester_edq == 68 | bimester_edq == 69
replace j087 =. if j087 <0
replace j139 =. if j139 <0
egen perserverance_1_sd = std(j087), by(bimester_edq)
egen perserverance_2_sd = std(j139), by(bimester_edq)

replace j085=. if j085<0
egen c_disruptive_in_class = std(j085), by(bimester_edq)

foreach var in  j127 j128 j129 j130 j131 j132 j133 j134 j135 j136 j137 j138  j140 j141 j142 j143 j144 j145 j146 j147 j148 j149 j150 j151 j152 j153 j154 j155 j156 j157 j158 j159 j160 j161 j162 j163 j164 j165 j166 j167 j168 j169 j170 j171 j172 j173 j174 j175 j176 j177 {
	replace `var'=. if `var'<0
	egen c_`var'_sd = std(`var'), by(bimester_edq)

	
}

* Perseverance is coming from the Teacher's questionnaire. It is measured alongside
* another 49 markers of overall "classroom behaviour". It turns out that
* perserverance correlates strongly with some of these. So waht I'm going to do 
* is to separate the effect of perserverance and "all the rest" by using a 
* principal component. Specifically, I pass all the other measures through PCA. 

pca c_j127_sd c_j128_sd c_j129_sd c_j130_sd c_j131_sd c_j132_sd c_j133_sd c_j134_sd c_j135_sd c_j136_sd c_j137_sd c_j138_sd c_j140_sd c_j141_sd c_j142_sd c_j143_sd c_j144_sd c_j145_sd c_j146_sd c_j147_sd c_j148_sd c_j149_sd c_j150_sd c_j151_sd c_j152_sd c_j153_sd c_j154_sd c_j155_sd c_j156_sd c_j157_sd c_j158_sd c_j159_sd c_j160_sd c_j161_sd c_j162_sd c_j163_sd c_j164_sd c_j165_sd c_j166_sd c_j167_sd c_j168_sd c_j169_sd c_j170_sd c_j171_sd c_j172_sd c_j173_sd c_j174_sd c_j175_sd c_j176_sd c_j177_sd

predict pcA pcB pcC pcD pcE pcF pcG pcH

correlate perserverance_1_sd perserverance_2_sd pcA pcB pcC pcD pcE pcF pcG pcH

* The measures of perserverance is stronlgy negatively correlated with PC-A and
* strongly positively correlated with PC-C. Therefore in our regressions we include
* principal components B, D, E, F, G and H
* Including the PC's drops the sample size from 9066 to 6306, although the 
* results remain the same. 

* 


* Now, child's selfreport on how good at different subjects.

foreach var of varlist k036 k037 k038 k039 k040 k041 k042 k043{
	
	gen good_at_`var' = `var' == 1
	replace good_at_`var' =. if `var' ==-4
}

foreach var of varlist k080 k088 k094 k075{
	gen luckvar_`var' = `var' ==1
	replace luckvar_`var' =. if `var' <0
}

drop _merge

* Outcomes at age 30

merge 1:1 bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-5558-stata/stata/stata13_se/bcs2000.dta", keepusing(bencod* seg sc cjorg cnetpay cnetprd chours1 cjperm jdemand* skill* econact  actagel2 cjsup2 )
drop _merge

merge 1:1 bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-5558-stata/stata/stata13_se/bcs6derived.dta", keepusing(HIACA00)

drop _merge
merge 1:1 bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-8547-stata/stata/stata13/bcs_age46_main.dta", keepusing(BD10ECACT B10INCP4 B10INCW B10FINNOW B10IASI01 B10IASI02 B10IASI03 B10IASI04 B10IASI05 B10IASI06 B10IASI07 B10IASI08 B10IASI09 B10IASI10 B10DEBT01 B10DEBT02 B10DEBT03 B10DEBT04 B10DEBT05 B10DEBT06 B10DEBT07 B10DEBT08 B10DEBT09 B10DEBT10 B10DEBT11 B10DEBT12 B10JBATIS B10INCC01 B10CJSUPA B10SOC3 B10SIC3 B10NSSECSB B10NSSECOP B10NSSECAN)

drop _merge
merge 1:1 bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-9347-stata/stata/stata13/bcs11_age51_main.dta", keepusing(b11iasi01 b11iasi02 b11iasi03 b11iasi04 b11iasi05 b11iasi06 b11iasi07 b11iasi08 b11savtotc b11debt01 b11debt02 b11debt03 b11debt04 b11debt05 b11debt06 b11debt07 b11debt08 b11debt09 b11debt10 b11debt11 b11expectret b11uslw b11netw b11econact2 b11beenvac b11numvac)

drop _merge
merge 1:1 bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-9347-stata/stata/stata13/savings.dta", keepusing(logTotalSavingsCapped51)
 
* age 51 
*===========================


*===========================

gen employed = econact == 1 |econact == 3 // Full time employed.
gen employedAt46 = BD10ECACT == 1 | BD10ECACT == 3
gen employedAt51 =  b11econact2 == 1 |  b11econact2 ==3

gen unemployedAnytime = econact== 5 | BD10ECACT ==5 | b11econact2== 5
 
gen profORmanagerialClass = sc == 1 | sc == 2

gen ismanager = cjsup2 == 1


gen ismanagerAt46 = B10CJSUPA == 1 
replace ismanagerAt46=. if B10CJSUPA <0

gen weeklyIncAt46 =B10INCW 
replace weeklyIncAt46 = . if B10INCW <0
gen logWeeklyInc46 = log(weeklyIncAt46 )

gen weeklyIncAt51 = b11netw
replace weeklyIncAt51=. if b11netw <0
gen logWeeklyInc51 = log(weeklyIncAt51)



gen noLiquidSavingsAt46 = B10IASI02==1 | B10IASI03==1| B10IASI04==1| B10IASI05==1| B10IASI06==1
gen noLiquidSavingsAt51 = b11iasi02 ==1 | b11iasi03 ==1 |  b11iasi04 ==1 |  b11iasi05 ==1 |  b11iasi06 ==1 

gen personalLoan46= B10DEBT04 ==1 | B10DEBT07 ==1
gen personalLoan51 = b11debt05 ==1 

gen noDebtAt46 = B10DEBT10==1


foreach var of varlist skill1a skill2a skill3a skill4a skill5a skill6a skill7a skill8a skill9a {
	
	gen skill_`var' = (`var' == 1 )
	replace skill_`var' =. if `var' >=8
}
drop _merge


gen noqualifications = HIACA00 == 0
gen degree = HIACA00 == 7 | HIACA00 ==8

gen interferes = jdemand3== 1
replace interferes = . if jdemand3==9

gen monthInc = cnetpay if cnetprd==4
replace monthInc = cnetpay /12 if cnetprd==5
replace monthInc = cnetpay*52 / 12 if cnetprd == 1
replace monthInc = cnetpay*26 / 12 if cnetprd == 2
replace monthInc = cnetpay*13 / 12 if cnetprd == 3
replace monthInc =. if cnetprd == 6
replace monthInc =. if monthInc >4000 | monthInc <0
gen loginc = log(monthInc)

tab sc, gen(social_class_job)

gen manual_worker = seg >=6



global set7 pcB pcD pcE pcF pcG pcH


// This wasn't used in the end 
/*
use "/Users/user/Documents/datasets/britCohortStudy70/UKDA-9347-stata/stata/stata13/bcs11_age51_savings_longf.dta", clear
replace b11savtot =. if b11savtot<0

collapse (sum) b11savtot, by(bcsid)
gen totalSavingsCapped51 = b11savtot
replace totalSavingsCapped51 = 300000 if totalSavingsCapped51>300000
gen logTotalSavingsCapped51= log(totalSavingsCapped51)

save "/Users/user/Documents/datasets/britCohortStudy70/UKDA-9347-stata/stata/stata13/savings.dta",replace

*/
