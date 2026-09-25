use "/Users/user/Documents/datasets/britCohortStudy70/UKDA-2666-stata/stata/stata13/bcs7072a.dta", clear

merge 1:1 bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-2666-stata/stata/stata13/bcs1derived.dta"

keep a0012 a0012a a0005a a0009 a0018 a0255 a0307 a0308  bcsid  a0010 a0013 a0014 a0017  a0043b a0164 a0242 a0278 a0195b BD1REGN BD1TEENM BD1MAGE BD1FAGE a0163 a0166 a0168 a0169 a0170 a0171

replace a0005a=. if a0005a <0

* Education of the mother. We first check that the Butler Act
* worked as expected.

gen mum_dob = 1970 - a0005a

preserve
gen left14 = a0009 == 14
gen left15 = a0009 ==15
collapse left* ,by(mum_d)
twoway(scatter left15 mu, xline(1934))
restore

gen mum_lftAfterSLA = (mum_dob <=1933 & a0009 >14 | mum_dob >1933 & a0009 >15)
replace mum_lftAfterSLA=1 if a0009 == 97 & a0005a >15


replace BD1FAGE = . if BD1FAGE <0
gen dad_dob = 1970 - BD1FAGE 

preserve
gen left14 = a0010 == 14
gen left15 = a0010 ==15
collapse left* ,by(dad_dob)
twoway(scatter left14 dad, xline(1934))
restore


gen dad_lftAfterSLA = (dad_dob <=1933 & a0010 >14 | dad_dob >1933 & a0010 >15)
replace dad_lftAfterSLA=1 if a0010 == 97 & BD1FAGE >15

* Region of residence

sum BD1REGN

* Father/Mother's social class.

replace a0014 = . if a0014<0
replace a0018 = . if a0018 <0

* Smoking during pregnancy
replace a0043b = . if a0043b <0 
gen mum_smoked = a0043b >=4

* Gestation age; we use post 37 pre and 42 post
gen normal_gestation = a0195b >37 & a0195b<42

* Gender baby 
replace a0255 = . if a0255 < 0

* Weight baby
replace a0278 = . if a0278 < 0

* Premarital conception
replace a0012a = . if a0012a <0

* Marital status at birth
replace a0012 = .  if a0012 < 0

* Parity
replace a0166 = . if a0166 <0
replace a0168 = . if a0168 <0


* Covariates

gen c_delivery_age_mum = a0005a
gen c_delivery_age_dad = BD1FAGE
gen c_stayedEd_mum = mum_lftAfterSLA
gen c_stayedEd_dad =dad_lftAfterSLA
gen c_region_birth = BD1REGN
gen c_socclass_mum = a0018 // Job class of mother before pregnancy
gen c_socclsss_dad = a0014 // job class of father at birth
gen c_mumsmoked = mum_smoked
gen c_abnormal_gest = normal_gestation ==0
gen c_sexbirth = a0255
gen c_birthweight = a0278
gen c_premarconception = a0012a == 1
replace c_premarconception = . if a0012a ==.
gen c_singleconception = a0012 != 2
replace c_singleconception = . if a0012 == .
gen c_parity = a0166
gen c_prevDeathUnder7 = a0170 >0
gen c_prevDeathAfter7 = a0171 >0
gen c_prevStillbirth = a0169 >0
gen c_prevMiscarries = a0168 >0

gen isInPerinatalWave = 1

keep bcsid c_*  isInPerinatalWave


// Covariates from 1975 (age 5). NOTE FOR FUTURE: IN f699a there are some intersting
// Variables on mother's attitude (authoritarian attitudes). THese can be used to evaluate
// future behaviour (including political behaviour). 


merge 1:1 bcsid using  "/Users/user/Documents/datasets/britCohortStudy70/UKDA-2699-stata/stata/stata13/f699b.dta"
 

keep bcsid e006 e007 e020 e038 e075 e076 e104 e105 e130 e189a e189b e191 e192 e193 e194 e197 e204b e204a e205 e206 e208 e209a e209b e217 e218 e220 e228a e228b e232 e234 e230 e229 e245 e246a e246b c_* isInPerinatalWave

* breast feeding

gen c_breastfed = e020 ==2 |e020 ==3

* ever inpatient

gen c_inpatientAt5 = e038 ==1

* Present or past hearing, sight or speech difficulty

gen c_hearSeeSpeak_issueAt5 = (e075 == 1 | e075 == 2 |e075 == 3 |e075 == 4) | ///
							(e076 == 1 | e076 == 2 |e076 == 3 |e076 == 4) | ///
							(e104 == 1 | e104 == 2 |e104 == 3) | ///
							(e105 == 1 | e105 == 2 |e105 == 3) 

* Child read at by father or mother

gen c_childReadAt_5 = (e130 == 1 | e130 == 2) 

* QUalifications

gen c_mum_nonqual = e189a == 1
gen c_mum_secondary = e189a >1 & e189a<6
gen c_mum_degree = e189a == 7

gen c_dad_nonqual = e189b == 1
gen c_dad_secondary = e189b >1 & e189b<6
gen c_dad_degree = e189b == 7

* MOther's occupation at 5
gen c_mumWorksFTAt5 = e208 >=20
gen c_mumWorksPTAt5 = e208 >0 & e208<20 


* Accommodation
replace e218 = . if e218<0
tab e218, gen(c_typeAcommodationAt5)

replace e220 =. if e220 <0
tab e220, gen(c_accTenureAt5)

replace e228a =. if e228a <0
gen c_roomNumAt5 = e228a

* Ethnicity 
replace e245 = . if e245< 0
gen c_nonwhite = e245 != 1

* number siblings
replace e006 =. if e006 < 0
replace e007 =. if e007 < 0
egen c_noSiblingsAt5 = rowtotal(e006 e007)

keep bcsid c_* isInPerinatalWave


merge 1:1 bcsid using "/Users/user/Documents/datasets/britCohortStudy70/UKDA-3723-stata/stata/stata13_se/sn3723.dta"

* Covariates from Parent questionnarie
*----------------------------------------------------------------

* Children in the house (in addition to CM)

replace a4a_42=. if a4a_42<0
gen c_childHhldAt10 = a4a_42
 
 * Child admitted to hospital
gen c_inpatientBefore10 = b16_1 == 1

* Fathers qualification
gen c_fatherNoQual10 = c1_9==1
gen c_fatherDegree10 = c1_6==1

* MOthers qualification
gen c_motherNoQual10 = c1_20==1 
gen c_motherDegree10 = c1_17==1

* Father unemployed
gen c_fatherUnemployed10 = c2_2 ==1 |c2_3==1 | c2_4==1 
gen c_fatherStayhome10 = c2_4==1
* Mother unemployed 

gen c_motherUnemployed10 = c2_10==1 | c2_11 ==1 |c2_13==1 
gen c_motherStayhome10 = c2_12==1

* Benefits: family income supple, supplementary benefit, sickeness, disablement, attendance allowance or unemployemnt 


gen c_benefit = c8_3==1 | c8_4==1 | c8_7 ==1 |c8_8 ==1 |c8_9==1 | c8_10==1 

* Selfreported income (in intervals, use middle of interval)
gen c_hholdIncomePW10 = 17.5
replace c_hholdIncomePW10 = 42 if c9_2==1
replace c_hholdIncomePW10=74.5 if c9_3 == 1
replace c_hholdIncomePW10=124.5 if c9_4 ==1
replace c_hholdIncomePW10=174.5 if c9_5 ==1
replace c_hholdIncomePW10=224.5 if c9_6 ==1 
replace c_hholdIncomePW10=274.5 if c9_7 ==1 

tab c_hholdIncomePW10, gen(c_incband)

* household tenure
tab d2, gen(c_accTenureAt10)

replace d5_1 =. if d5_1 <0
gen c_roomNumAt10 = d5_1

* Father BMI 
gen c_fatherBMIat10 = e2_2 / ((e2_1/100)^2)
replace c_fatherBMIat10 =. if e2_2 <0 | e2_1 <0

* Mother BMI 
gen c_motherBMIat10 = e1_2 / ((e1_1/100)^2)
replace c_motherBMIat10 =. if e1_2 <0 | e1_1 <0

* Mother drunk during pregnancy
gen c_motherDrunk = e4_2 <=3 

* Mother smokes 
gen c_motherSmokesat10 = e9_1== 1 | e9_1==2

* Father smokes
gen c_fatherSmokesat10 = e11_1 == 1 | e11_1==2


* Medical questionnarie
*----------------------------------------------------------------

* Age of child at exam examdatb examdatc

gen dob_s = ym( 1970, 4)
replace examdatc=. if examdatc<0
replace examdatb=. if examdatb<0
replace examdatc = 1980 if examdatc==80
replace examdatc = 1981 if examdatc==81
gen exam_date = ym(examdatc,examdatb)

gen bimester_f = floor(exam_date/2)



* Eyesight: some interfering defect
gen c_eyesightIssue10 = meb11_1 >=1 & meb11_1<7

* Hearing issues
gen c_hearingIssue10 = meb12_15 >=1 & meb12_15<=2

* Speech
gen c_speechIssue10 = meb15_1 ==3 | meb15_1 ==4

* Height to age

replace meb17 =. if meb17<0
egen c_heightAt10_sd = std(meb17), by(bimester_f)

* Weight to age (too many missing values)

replace meb19_1 =. if meb19_1<0
egen c_weightAt10_sd = std(meb19_1), by(bimester_f)




* Mother's questionnarie
*----------------------------------------------------------------

* Family co-activities

forvalues i =107/113{
	
	gen act`i'= m`i' == 1
}
egen c_disengagedAt10 = rowtotal(act107 act108 act109 act110 act111 act112 act113)


* Teacher questionnarie
*----------------------------------------------------------------

*school type

gen c_schoolMaintainedAt10 = j251==1
gen c_schoolVolCtrlAt10 = j252==1
gen c_schoolVolAidedAt10 = j253==1
gen c_schoolDirectGrantAt10 = j254==1
gen c_schoolIndepAt10 = j255==1

drop bimester_f



/* Define sets of covariates.
*----------------------------------------------------------------
Set 1: Key variables

	Region of birth
	Delivery age of the mother
	Mother smoked during pregnancy
	Abnormal gestation (too short or too long)
	Sex
	Weight at birth
	

Set 2: of perinatal information:

From the perinatal questionnarie:
	
	Mum stayed in education beyond the SLA
	Dad stayed in education beyond the SLA
	Parity
	previous death of baby under 7 days
	previous death of baby over 7 days
	previous still births
	previous miscarriages
	
Set 3: Age 5
	
	Breastfed baby
	Inpatient event
	Hear/see/speech issues
	child read at by mum/dad
	mother no qualifications
	mother secondary education
	mother degree 
	Idem father
	Mother works full time
	mother works part time
	Non-white
	Number of siblings
	

Set 4: Age 10

	Impatient episode
	Father/mother unemployed
	Father/mother stay home
	family receives benefit
	household income per week 
	Tenure
	number of rooms
	mother smokes 
	father smokes

Set 5: Medical
	Issues with sight, hearing or speech.
	height standardised
	weight standardised
	
Set 6: Teacher's questionnaire
	type of school

Set 7 Mothers questionnaire 

	Disengaged

*/

gen c_delivery_age_mum_2 = c_delivery_age_mum^2
global set1 c_delivery_age_mum c_delivery_age_mum_2 i.c_region_birth c_mumsmoked c_abnormal_gest i.c_sexbirth c_birthweight 

global set2 c_stayedEd_mum c_stayedEd_dad c_parity c_prevDeathUnder7 ///
	c_prevDeathAfter7 c_prevStillbirth c_prevMiscarries 
	
global set3   c_breastfed c_inpatientAt5 c_hearSeeSpeak_issueAt5 c_childReadAt_5 c_mum_nonqual c_mum_secondary c_mum_degree c_dad_nonqual c_dad_secondary c_dad_degree c_mumWorksFTAt5 c_mumWorksPTAt5   c_nonwhite c_noSiblingsAt5


global set4 c_childHhldAt10 c_inpatientBefore10  c_fatherUnemployed10 c_fatherStayhome10 c_motherUnemployed10 c_motherStayhome10 c_benefit c_incband* c_accTenureAt101 c_accTenureAt102 c_accTenureAt103 c_accTenureAt104 c_accTenureAt105 c_accTenureAt106 c_accTenureAt107 c_roomNumAt10   c_motherSmokesat10 c_fatherSmokesat10

global set5 c_eyesightIssue10 c_hearingIssue10 c_speechIssue10 c_heightAt10_sd c_weightAt10_sd

global set6 c_disengagedAt10


* Generate an adversity index.



sum c_birthweight, det
gen belowMedianBirthweight = c_birthweight < r(p50)

gen incomeDisadvantage = c_incband1 ==1 | c_incband2==1 | c_incband3 ==1 


egen shitindex = rowtotal(c_abnormal_gest   c_inpatientAt5 c_hearSeeSpeak_issueAt5 c_mum_nonqual c_dad_nonqual c_inpatientBefore10 c_fatherUnemployed10 c_motherUnemployed10 c_benefit c_eyesightIssue10 c_hearingIssue10 c_speechIssue10 c_disengagedAt10 belowMedianBirthweight incomeDisadvantage), missing

sum shitindex, det

gen shitlife = shitindex >=5  // Above median amount of shit life.
