*Calculating vulnerability based on weekly information: Counting Approach.
*This version also includes the dimensional contribution
*Keep cutoffs the same as, 0, 0.5. 0.65, 0.80, 0.95
*COnsidered three different ways of calculating vulnerability
*Separated identification and aggregation
*HHsize taken in to acccount

global path  "PLEASE USE YOUR GLOBAL PATH"

use "${path}\PLEASE PUT YOUR DATA FILE NAME", clear



set more off
preserve
drop ref 
gen weekid=week(date)
gen yearid=year(date)
gen monthid=month(date)
gen qid=quarter(date)

*xtset hhid time0

drop if yearid==2019
drop if hhid==6  // very few observations
replace hhid=hhid-1 if hhid>6
drop if hhid==69 // very few observations
replace hhid= hhid-1 if hhid>69

replace hhsize=round(hhsize)  //frequency weights do not take decimals
 

*creating income variable
gen expd =0 // expenditure on legitimate business 
replace expd=-(takaout) if class=="EEM"|class=="ESH"|class=="EFA"|class=="EBU"|class=="EVE"|class=="EWO"
*note that takaout already is in negative in the data
*EEM	emoney stock purchases
*ESH	shop stock purchases
*EFA	farm inputs
*EBU	business asset purchase
*EVE	vehicle repair costs
*EWO	work tools and materials

gen income =  takain - expd // taking legitimate business expenditure out of income
replace income =0 if income ==.
replace income = income/hhsize // taking household size in to consideration


*Construct weekly per capita income 
egen weekincR=sum(income), by(hhid yearid weekid)
replace weekincR=1 if weekincR<=0
gen weekincLR=ln(weekincR)


egen tag1=tag(hhid yearid weekid)

bysort hhid: gen time = sum(tag1)

**** Number of Weeks

*number of weeks of observation for each household in each year
bysort hhid yearid weekid:  gen wkcount = _n == 1 
by hhid yearid: gen wkcount1 = sum(wkcount)
by hhid yearid: replace wkcount = sum(wkcount)
by hhid yearid: replace wkcount = wkcount[_N]  
by hhid yearid: replace wkcount = 0 if wkcount == .

*table hhid yearid, stat(max wkcount) nototals


 ****Summary Statistic based on the whole sample
 
gen yid=1 if yearid==2015
replace yid=2 if yearid==2016
replace yid=3 if yearid==2017
replace yid=4 if yearid==2018
*replace yid=5 if yearid==2019

**************summary of full sample******************

mat HINCA= J(6,4,0)

forval j=1/4 {
	
									
					quietly distinct hhid if  yid==`j' 
					scalar n`j'= r(ndistinct) // number of households each year
					
					quietly su weekincR  if  tag1==1 & yid==`j' [fw=hhsize], d
					mat HINCA[1,`j'] = r(mean) // average income per week per person
					mat HINCA[2,`j']= r(sd) // standard deviation
					mat HINCA[3,`j']= r(p50) // median
					
					quietly su weekincR  if  tag1==1 & yid==`j'
					scalar a`j'=r(N)/n`j' // avg observations per household per year
					mat HINCA[4,`j']= a`j'
				    
					
					mat HINCA[5,`j'] = n`j' // number of households
					
	
					quietly su hhsize if tag1==1 & yid==`j'
					mat HINCA[6,`j'] =r(mean) // avg household size
					
				}			

*** Income is per capita per week
mat colnames HINCA = 2015 2016 2017 2018 
mat rownames HINCA = MeanInc SDIncome MedianInc AvgweeklyObs NHH AvgHHSize

**** creating matrix for weeks
quietly su hhid if tag1==1
scalar nh1=r(max)

mat HWK= J(`=nh1', 4, 0)

forval i=1/`=nh1'{
			
			forval j=1/4{
				 quietly su wkcount  if  tag1==1 & yid==`j' & hhid==`i'
				 mat HWK[`i',`j'] = r(max)
						}
			}

mat colnames HWK = 2015 2016 2017 2018 			
*matlist HWK
*/
**********************Calculating the the contribution of shocks for 2017 DISTRIBUTION APPROACH******************************

**** Dropping observations that has less than 50 weeks by 2017

drop if hhid>49 

*
*****************Summary for sample in 2017****************
mat HINC= J(6,2,0)
	
									
					quietly distinct hhid if  yearid<=2017 
					scalar n17= r(ndistinct) // number of households each year
					
					quietly su weekincR  if  tag1==1 & yearid<=2017  [fw=hhsize], d
					mat HINC[1,1] = r(mean) // average income per week per person
					mat HINC[2,1]= r(sd) // standard deviation
					mat HINC[3,1]= r(p50) // median
					
					quietly su weekincR  if  tag1==1 & yearid<=2017 
					mat HINC[4,1]= r(N)/n17 // avg observations per household per year
				    
					
					mat HINC[5,1] = n17 // number of households
					
	
					quietly su hhsize if tag1==1 & yearid<=2017 
					mat HINC[6,1] =r(mean) // avg household size
					
							

*** Income is per capita per week
mat colnames HINC = 2017 2018 
mat rownames HINC = MeanInc SDIncome MedianInc AvgweeklyObs NHH AvgHHSize
*/


sort hhid yearid weekid 

keep if tag1==1
				
quietly su hhid if yearid<=2017 & tag1==1
scalar nh=r(max)


mat PCoeff = J(`=nh', 5,0)

tsset hhid time, weekly


forval i=1/`=nh' {
					
			
					eststo m`i': quietly prais weekincLR i.qid time if hhid==`i' & yearid<=2017 & tag1==1, rho(regress) vce(robust) 
				
					*eststo m`i':quietly reg weekincLR L.weekincLR if hhid==`i' & yearid<=2017 & tag1==1, vce(robust) // use this for AR1 results presented in the appendix
					
					
					mat PCoeff[`i',1]=_b[_cons]
					mat PCoeff[`i',2]=_se[_cons]
					mat PCoeff[`i',3]=_b[time]
					mat PCoeff[`i',4]=_se[time]
					*mat PCoeff[`i',3]=_b[L.weekincLR] // relevant for AR1
					*mat PCoeff[`i',4]=_se[L.weekincLR] // relevant for AR1
					mat PCoeff[`i',5]= e(N)
					
					scalar rmse`i'=e(rmse)
					
					predict mB`i' if hhid==`i' & yearid<=2017 & tag1==1, xb 
					quietly su mB`i' if hhid==`i' & yearid<=2017 & tag1==1 
					scalar meanB`i' = r(mean)
					
					predict sdB`i' if hhid==`i' & yearid<=2017 & tag1==1, stdp
					quietly su sdB`i' if hhid==`i' & yearid<=2017 & tag1==1 
					scalar meansdB`i'=r(mean)
					
					**** Take both variance of e and variance of mean forcast (see Woolridge Chp6.6)***
					scalar stddevB`i' = ((meansdB`i')^2+ (rmse`i')^2)^(0.5)
					
					set obs 600000
					set seed 12345
				
					gen weekincB`i' =rnormal(meanB`i',stddevB`i')
					replace weekincB`i' = exp(weekincB`i') //because we are generating log income, hence the income in taka needs exponential
					replace weekincB`i' = 0 if weekincB`i' ==1 // this is because earlier we adjusted any zero income as 1 to facilitate taking log values
	
				}
				

matlist PCoeff				
				
*esttab m* using "${path}\results.csv", keep(_cons time) varlabels($VARNAMES) starlevels(* 0.05 ** 0.01 *** 0.001) b(%5.3f) se(%5.3f) replace

*esttab m* using "${path}\results.csv", keep(_cons L1.) varlabels($VARNAMES) starlevels(* 0.05 ** 0.01 *** 0.001) b(%5.3f) se(%5.3f) replace							

mat P1=J(`=nh', 4, 0) // stores probability of different states
mat P2=J(`=nh', 5, 0) // average income in different states and household size
mat P3=J(`=nh', 4, 0) // probabibility weighted deprivation in different states
mat P4=J(`=nh', 4, 0) // probability weighted squared deprivation
mat P5=J(`=nh', 4, 0) // probability conditional on deprivation being positive
mat P6=J(`=nh', 4, 0)
mat P7=J(`=nh', 4, 0) // mainly to store total deprivation to calculate 


*mat Pn=J(50,1,0)

	forval i= 1/`=nh' {
						
							*Calculating probabilities

							quietly su weekincB`i'
							scalar t`i'=r(N)
							quietly count if weekincB`i'<=482.84
							scalar x1`i'=r(N)
							quietly count if weekincB`i'>482.84 & weekincB`i'<=813.22 
							scalar x2`i'=r(N)
							quietly count if weekincB`i'>813.22 & weekincB`i'<=1397.72 
							scalar x3`i'=r(N)
							quietly count if weekincB`i'>1397.72 
							scalar x4`i'=r(N)
							matrix P1[`i',1]=x1`i'/t`i'
							matrix P1[`i',2]=x2`i'/t`i'				
							matrix P1[`i',3]=x3`i'/t`i'
							matrix P1[`i',4]=x4`i'/t`i'
		****P1 above calculates the frequency of different states for different households
			*** calculating average income in each of the shocks				
							quietly sum weekincB`i' if weekincB`i'<=482.84, meanonly
							scalar y1`i'=r(mean)
							if (r(mean)<0) {
								scalar drop y1`i'
								scalar y1`i'=0 // ensuring that lowest income from shocks is zero                         
								}
							
							else if (r(mean)==.) {
								scalar drop y1`i'
								scalar y1`i'=0 
										}
							
							quietly sum weekincB`i' if weekincB`i'>482.84 & weekincB`i'<813.22, meanonly
							scalar y2`i' =r(mean)
							if (r(mean)==.) {
								scalar drop y2`i'
								scalar y2`i'=0 
										}
							quietly sum weekincB`i' if weekincB`i'>813.22 & weekincB`i'<=1397.72 , meanonly
							scalar y3`i' =r(mean)
							if (r(mean)==.) {
								scalar drop y3`i'
								scalar y3`i'=0 
										}
							quietly sum weekincB`i' if weekincB`i'>1397.72, meanonly
							scalar y4`i' =r(mean)
								if (r(mean)==.) {
								scalar drop y4`i'
								scalar y4`i'=0 
										}
				
							matrix P2[`i',1]=y1`i'
							matrix P2[`i',2]=y2`i'
							matrix P2[`i',3]=y3`i'
							matrix P2[`i',4]=y4`i'
							
						   quietly sum hhsize if hhid==`i' // household size
								mat P2[`i',5]=r(mean)	
							
				*** probability weighted deprivations by states		
							
							if (y1`i'>=0){ // extreme shock zero income is high deprivation
								mat P3[`i',1]=((1397.72-y1`i')/1397.72)*P1[`i',1] 
										}
				
							if (y2`i'>0){
								mat P3[`i',2]=((1397.72-y2`i')/1397.72)*P1[`i',2] 
										}
										
							if (y3`i'>0){
								mat P3[`i',3]=((1397.72-y3`i')/1397.72)*P1[`i',3]
										}
				
								mat P3[`i',4]= P3[`i',1]+P3[`i',2]+P3[`i',3] // probability weighted total deprivation
							*
							
					**** squared poverty gap
							if (y1`i'>=0){
								mat P4[`i',1]=((1397.72-y1`i')/1397.72)^2*P1[`i',1] 
										}
				
							if (y2`i'>0){
								mat P4[`i',2]=((1397.72-y2`i')/1397.72)^2*P1[`i',2] 
										}
										
							if (y3`i'>0){
								mat P4[`i',3]=((1397.72-y3`i')/1397.72)^2*P1[`i',3]
										}
				
							mat P4[`i',4]= P4[`i',1]+P4[`i',2]+P4[`i',3] // probability weighted total squared deprivation
							
							
						***** probability conditional on deprivation being positive
						
							if (y1`i'>=0){
								mat P5[`i',1]=P1[`i',1] 
										}
				
							if (y2`i'>0){
								mat P5[`i',2]=P1[`i',2] 
										}
										
							if (y3`i'>0){
								mat P5[`i',3]=P1[`i',3]
										}
				
							mat P5[`i',4]= P5[`i',1]+P5[`i',2]+P5[`i',3] //total probability of falling into poverty
							
						***** P7 for the standard calculations***********
							
							if P3[`i',4]>0.5 {
							
										mat P7[`i',1] = 1
										mat P7[`i',2] = P3[`i',4]  //poverty gap
										mat P7[`i',3] = P4[`i',4]  //squared poverty gap
										
								}		
							
										mat P7[`i',4]= P2[`i',5]  // hhsize
						
							
				}
	

mat VulA=J(6,8,0) // calculate vulnerability poverty gap in 2017 and 2018
matrix colnames VulA = N HC Vul17 Intens17 N HC Vul18 Intens18
matrix  rownames VulA = R0 R50 R65 R80 R95 StdM

mat Shk=J(5,6,0) // calculate shocks
matrix colnames Shk = Ext17 Mod17 Mild17 Ext18 Mod18 Mild18
matrix  rownames Shk = R0 R50 R65 R80 R95	

mat VulB=J(6,8,0) // calculate vulnerability squared in 2017 and 2018
matrix colnames VulB = HC Vul17 Intens17 Med17 HC Vul18 Intens18 Med18
matrix  rownames VulB = R0 R50 R65 R80 R95 StdM
	

	
	forval k=5(-1)1 {
			    
				*scalar theta`k'= 1.25-`k'*0.25  // k times 0.25 deducted from 1 give the cut off values of 0, 0.25. 0.5, 0.75, 1
				
				if `k'==5 {
					scalar theta`k'= 0
				}
				else if `k'==4 {
					scalar theta`k'= 0.5
				}
				else if `k'==3 {
					scalar theta`k'= 0.65
				}
				else if `k'==2 {
					scalar theta`k'= 0.8
				}
				else if `k'==1 {
					scalar theta`k'= 0.95
				}
				
		
				mat P6`k'=J(`=nh', 4, 0)
				mat P6S`k'=J(`=nh', 5, 0)
				forval i= 1/`=nh' {	
						
								mat P6`k'[`i',4]= P2[`i',5]  // ensuring household size if part of P6 for weighting purposes
								mat P6S`k'[`i',5]= P2[`i',5]	// hhsize included for weighting under summary 
									
									if P5[`i',4]>=theta`k'{ 
												
												mat P6`k'[`i',1]=1 //counting those who are deprived
												mat P6`k'[`i',2]=P3[`i',4] //poverty gap
												mat P6`k'[`i',3]=P4[`i',4] // sq pov gap
												
												mat P6S`k'[`i',1]=P3[`i',1]  //extreme shock
												mat P6S`k'[`i',2]=P3[`i',2] //moderate shock
												mat P6S`k'[`i',3]=P3[`i',3]  // mild shock
												mat P6S`k'[`i',4]=P3[`i',4]  // total shock
														}
							}
					
				scalar r`k'= 6-`k' // to ensure that results are stacked from the top to the bottom since k start from 5
				
				 *** calculating for vul gap and squared vul gap
				svmat P6`k', names(c`k')
				
				quietly summ c`k'1  [fw=c`k'4], meanonly
				mat VulA[r`k', 1] =r(sum) //number of individuals who are vuln.
				mat VulA[r`k', 2]= r(mean) // headcount
				mat VulB[r`k', 1]= r(mean) // headcount
				
				quietly summ c`k'2 [fw=c`k'4], d
				mat VulA[r`k', 3] =r(mean) //avg vul
				mat VulA[r`k', 4] =r(sum)/VulA[r`k', 1] // intensity
				*mat VulA[r`k', 5]= r(p50) // median vulnerability
				
				
				quietly summ c`k'3 [fw=c`k'4], d
				mat VulB[r`k', 2] =r(mean) //avg vul
				mat VulB[r`k', 3] =r(sum)/VulA[r`k', 1] // intensity
				mat VulB[r`k', 4]= r(p50) // median vulnerability
				
				
			********************Calculating the shocks for 2017**********	
				
				svmat P6S`k', names(s`k')
				
				quietly summ s`k'4  [fw=s`k'5]  // total deprivation from all shocks
				scalar td`k'= r(sum)
				
				quietly summ s`k'1  [fw=s`k'5]
				mat Shk[r`k', 1] =(r(sum)/td`k')*100 //percentage of shock from extreme
				
				
				quietly summ s`k'2 [fw=s`k'5]
				mat Shk[r`k', 2] =(r(sum)/td`k')*100 //percentage of shock from moderate
				
				quietly summ s`k'3 [fw=s`k'5]
				mat Shk[r`k', 3] =(r(sum)/td`k')*100 //percentage of shock from mild
						
				}

********Calculating Standard case for 2017**********************
				svmat P7, names(w)
				
				quietly summ w1  [fw=w4], meanonly
				mat VulA[6, 1] =r(sum) //number of individuals who are vuln.
				mat VulA[6, 2]= r(mean) // headcount
				mat VulB[6, 1]= r(mean) 
				
				quietly summ w2 [fw=w4], d
				mat VulA[6, 3] =r(mean) //avg vul
				mat VulA[6, 4] =r(sum)/VulA[6, 1] // intensity
				*mat VulA[6, 5]= r(p50) // median vulnerability
				
				
				quietly summ w3 [fw=w4], d
				mat VulB[6, 2] =r(mean) //avg vul
				mat VulB[6, 3] =r(sum)/VulA[6, 1] // intensity
				mat VulB[6, 4]= r(p50) // median vulnerability
										
						

scalar drop _all

mat drop P1 P2 P3 P4 P5 P7

mata

 mata drop  P6*  ///

 end

restore


**********************Calculating for 2018 based on the DISTRIBUTION based approach**************************

preserve
drop ref  
gen weekid=week(date)
gen yearid=year(date)
gen monthid=month(date)
gen qid=quarter(date)

**** Dropping done previously for 2017 same issues hold for 2018
drop if yearid==2019

drop if hhid==6
replace hhid=hhid-1 if hhid>6
drop if hhid==69
replace hhid= hhid-1 if hhid>69


***** We are dropping households for which we do not have information in 2018

drop if hhid==4 // obs 4 no data for 2018
replace hhid=hhid-1 if hhid>4
drop if hhid==7  //obs 8 no data for 2018
replace hhid= hhid-1 if hhid>7
drop if hhid==12 //obs 14 no data for 2018
replace hhid= hhid-1 if hhid>12
drop if hhid==21 //obs 24 no data for 2018
replace hhid= hhid-1 if hhid>21
drop if hhid==23 //obs 27 no data for 2018
replace hhid= hhid-1 if hhid>23
drop if hhid==31 // obs 36 no data for 2018
replace hhid= hhid-1 if hhid>31
drop if hhid==34 //obs 40 no data for 2018
replace hhid= hhid-1 if hhid>34
drop if hhid==36 // obs 43 no data for 2018
replace hhid= hhid-1 if hhid>36
drop if hhid==39 // obs 47 no data for 2018
replace hhid= hhid-1 if hhid>39
drop if hhid==39 //obs 48 no data for 2018
replace hhid= hhid-1 if hhid>39
drop if hhid==57 //obs 67 deleted because less than 50 weeks info
replace hhid= hhid-1 if hhid>57



replace hhsize=round(hhsize)  //frequency weights do not take decimals
 

*creating income variable
gen expd =0 // expenditure on legitimate business 
replace expd=-(takaout) if class=="EEM"|class=="ESH"|class=="EFA"|class=="EBU"|class=="EVE"|class=="EWO"

gen income =  takain - expd // taking legitimate business expenditure out of income
replace income = income/hhsize // taking household size in to consideration


*Construct weekly per capita income 
egen weekincR=sum(income), by(hhid yearid weekid)
replace weekincR=1 if weekincR<=0
gen weekincLR=ln(weekincR)


egen tag1=tag(hhid yearid weekid)


*********for the time trend********

bysort hhid: gen time = sum(tag1)

***********summary of sample in 2018**************
									
					quietly distinct hhid if  yearid<=2018 
					scalar n18= r(ndistinct) // number of households each year
					
					quietly su weekincR  if  tag1==1 & yearid<=2018  [fw=hhsize], d
					mat HINC[1,2] = r(mean) // average income per week per person
					mat HINC[2,2]= r(sd) // standard deviation
					mat HINC[3,2]= r(p50) // median
					
					quietly su weekincR  if  tag1==1 & yearid<=2018 
					mat HINC[4,2]= r(N)/n18  // avg observations per household per year
				    
					
					mat HINC[5,2] = n18 // number of households
					
	
					quietly su hhsize if tag1==1 & yearid<=2018 
					mat HINC[6,2] =r(mean) // avg household size
					

sort hhid yearid weekid 

keep if tag1==1
				
quietly su hhid if yearid<=2018 & tag1==1
scalar nh=r(max)

mat PCoeff1 = J(`=nh', 5,0)

tsset hhid time, weekly

forval i=1/`=nh' {
					
					 quietly prais weekincLR i.qid time if hhid==`i' & yearid<=2018 & tag1==1, rho(regress) vce(robust)
					
					*eststo m`i':quietly reg weekincLR L.weekincLR if hhid==`i' & yearid<=2018 & tag1==1, vce(robust) // use this for the AR1 results in the appendix
				
					scalar rmse`i'=e(rmse)
					
					
					mat PCoeff1[`i',1]=_b[_cons]
					mat PCoeff1[`i',2]=_se[_cons]
					mat PCoeff1[`i',3]=_b[time]
					mat PCoeff1[`i',4]=_se[time]
					*mat PCoeff1[`i',3]=_b[L.weekincLR] //use this for AR1
					*mat PCoeff1[`i',4]=_se[L.weekincLR] //use this for AR1
					mat PCoeff1[`i',5]= e(N)
					
					
					
					predict mB`i' if hhid==`i' & yearid<=2018 & tag1==1, xb 
					quietly su mB`i' if hhid==`i' & yearid<=2018 & tag1==1 
					scalar meanB`i' = r(mean)
					
					***** we use stdf instead of stdp becuase stdp is just error of the mean rather than y..see Woolridge 6.4 (wisconsin notes) ****
					predict sdB`i' if hhid==`i' & yearid<=2018 & tag1==1, stdp
					quietly su sdB`i' if hhid==`i' & yearid<=2018 & tag1==1 
					scalar meansdB`i' = r(mean)
					
					scalar stddevB`i' = ((meansdB`i')^2+ (rmse`i')^2)^(1/2)
					
					set obs 600000
					set seed 12345
					
					gen weekincB`i' =rnormal(meanB`i',stddevB`i')
					replace weekincB`i' = exp(weekincB`i') //because we are generating log income, hence the income in taka needs exponential
					replace weekincB`i' = 0 if weekincB`i' ==1 // this is because earlier we adjusted any zero income as 1 to facilitate taking log values
				
				}

mat P1=J(`=nh', 4, 0) // stores probability of different states
mat P2=J(`=nh', 5, 0) // average income in different states and household size
mat P3=J(`=nh', 4, 0) // probabibility weighted deprivation in different states
mat P4=J(`=nh', 4, 0) // probability weighted squared deprivation
mat P5=J(`=nh', 4, 0) // probability conditional on deprivation being positive
*mat P6=J(`=nh', 4, 0)
mat P7=J(`=nh', 4, 0) // mainly to store total deprivation to calculate 


*mat Pn=J(50,1,0)

	forval i= 1/`=nh' {
						
							*Calculating probabilities
							*quietly count if weekincB`i'>=0 
							quietly su weekincB`i'
							scalar t`i'=r(N)
							*quietly count if weekincB`i'<=509.62 & weekincB`i'>=0
							quietly count if weekincB`i'<=509.62 
							scalar x1`i'=r(N)
							quietly count if weekincB`i'>509.62 & weekincB`i'<=858.30 
							scalar x2`i'=r(N)
							quietly count if weekincB`i'>858.30 & weekincB`i'<=1475.21 
							scalar x3`i'=r(N)
							quietly count if weekincB`i'>1475.21 
							scalar x4`i'=r(N)
							matrix P1[`i',1]=x1`i'/t`i'
							matrix P1[`i',2]=x2`i'/t`i'				
							matrix P1[`i',3]=x3`i'/t`i'
							matrix P1[`i',4]=x4`i'/t`i'
		****P1 above calculates the frequency of different states for different households
			*** calculating average income in each of the shocks				
							quietly sum weekincB`i' if weekincB`i'<=509.62, meanonly
							scalar y1`i'=r(mean)
							if (r(mean)<0) {
								scalar drop y1`i'
								scalar y1`i'=0 // ensuring that lowest income from shocks is zero                         
								}
							
							else if (r(mean)==.) {
								scalar drop y1`i'
								scalar y1`i'=0 
										}
							
							quietly sum weekincB`i' if weekincB`i'>509.62 & weekincB`i'<=858.30, meanonly
							scalar y2`i' =r(mean)
							if (r(mean)==.) {
								scalar drop y2`i'
								scalar y2`i'=0 
										}
							quietly sum weekincB`i' if weekincB`i'>858.30 & weekincB`i'<=1475.21 , meanonly
							scalar y3`i' =r(mean)
							if (r(mean)==.) {
								scalar drop y3`i'
								scalar y3`i'=0 
										}
							quietly sum weekincB`i' if weekincB`i'>1475.21, meanonly
							scalar y4`i' =r(mean)
								if (r(mean)==.) {
								scalar drop y4`i'
								scalar y4`i'=0 
										}
				
							matrix P2[`i',1]=y1`i'
							matrix P2[`i',2]=y2`i'
							matrix P2[`i',3]=y3`i'
							matrix P2[`i',4]=y4`i'
							
						   quietly sum hhsize if hhid==`i' // household size
								mat P2[`i',5]=r(mean)	
							
				*** probability weighted deprivations by states		
							
							if (y1`i'>=0){ // extreme shock zero income is high deprivation
								mat P3[`i',1]=((1475.21-y1`i')/1475.21)*P1[`i',1] 
										}
				
							if (y2`i'>0){
								mat P3[`i',2]=((1475.21-y2`i')/1475.21)*P1[`i',2] 
										}
										
							if (y3`i'>0){
								mat P3[`i',3]=((1475.21-y3`i')/1475.21)*P1[`i',3]
										}
				
								mat P3[`i',4]= P3[`i',1]+P3[`i',2]+P3[`i',3] // probability weighted total deprivation
							*
							
					**** squared poverty gap
							if (y1`i'>=0){
								mat P4[`i',1]=((1475.21-y1`i')/1475.21)^2*P1[`i',1] 
										}
				
							if (y2`i'>0){
								mat P4[`i',2]=((1475.21-y2`i')/1475.21)^2*P1[`i',2] 
										}
										
							if (y3`i'>0){
								mat P4[`i',3]=((1475.21-y3`i')/1475.21)^2*P1[`i',3]
										}
				
							mat P4[`i',4]= P4[`i',1]+P4[`i',2]+P4[`i',3] // probability weighted total squared deprivation
							
							
						***** probability conditional on deprivation being positive
						
							if (y1`i'>=0){
								mat P5[`i',1]=P1[`i',1] 
										}
				
							if (y2`i'>0){
								mat P5[`i',2]=P1[`i',2] 
										}
										
							if (y3`i'>0){
								mat P5[`i',3]=P1[`i',3]
										}
				
							mat P5[`i',4]= P5[`i',1]+P5[`i',2]+P5[`i',3] //total probability of falling into poverty
							
						***** P7 for the standard calculations***********
							
							if P3[`i',4]>0.5 {
							
										mat P7[`i',1] = 1  // indicates which households are free
										mat P7[`i',2] = P3[`i',4]  //poverty gap
										mat P7[`i',3] = P4[`i',4]  //squared poverty gap
										
								}	
							
										mat P7[`i',4]= P2[`i',5]  // hhsize
						
							
				}
	

	forval k=5(-1)1 {
			    
				*scalar theta`k'= 1.25-`k'*0.25  // k times 0.25 deducted from 1 give the cut off values of 0, 0.5. 0.65, 0.80, 0.95
				
				if `k'==5 {
					scalar theta`k'= 0
				}
				else if `k'==4 {
					scalar theta`k'= 0.5
				}
				else if `k'==3 {
					scalar theta`k'= 0.65
				}
				else if `k'==2 {
					scalar theta`k'= 0.80
				}
				else if `k'==1 {
					scalar theta`k'= 0.95
				}
				
				
				
				
				
				mat P6`k'=J(`=nh', 4, 0)
				mat P6S`k'=J(`=nh', 5, 0)
				forval i= 1/`=nh' {	
						
								mat P6`k'[`i',4]= P2[`i',5]  // ensuring household size if part of P6 for weighting purposes
								mat P6S`k'[`i',5]= P2[`i',5]	// hhsize included for weighting under summary 
									
									if P5[`i',4]>=theta`k'{ 
												
												mat P6`k'[`i',1]=1 //counting those who are deprived
												mat P6`k'[`i',2]=P3[`i',4] //poverty gap
												mat P6`k'[`i',3]=P4[`i',4] // sq pov gap
												
												mat P6S`k'[`i',1]=P3[`i',1]  //extreme shock
												mat P6S`k'[`i',2]=P3[`i',2] //moderate shock
												mat P6S`k'[`i',3]=P3[`i',3]  // mild shock
												mat P6S`k'[`i',4]=P3[`i',4]  // total shock
														}
							}
					
				scalar r`k'= 6-`k' // to ensure that results are stacked from the top to the bottom since k start from 5
				
 **** calculating for vul gap and squared vul gap for 2018************
 
				svmat P6`k', names(c`k')
				
				quietly summ c`k'1  [fw=c`k'4], meanonly
				mat VulA[r`k', 5] =r(sum) //number of individuals who are vuln.
				mat VulA[r`k', 6]= r(mean) // headcount
				mat VulB[r`k', 5]= r(mean) // headcount
				
				quietly summ c`k'2 [fw=c`k'4], d
				mat VulA[r`k', 7] =r(mean) //avg vul
				mat VulA[r`k', 8] =r(sum)/VulA[r`k', 5] // intensity
				*mat VulB[r`k', 5]= r(p50) // median vulnerability
				
				
				quietly summ c`k'3 [fw=c`k'4], d
				mat VulB[r`k', 6] =r(mean) //avg vulSq
				mat VulB[r`k', 7] =r(sum)/VulA[r`k', 5] // intensitySq
				mat VulB[r`k', 8]= r(p50) // median vulnerabilitySq
				
				
			********************Calculating the shocks for 2018**********	
				
				svmat P6S`k', names(s`k')
				
				quietly summ s`k'4  [fw=s`k'5]  // total deprivation from all shocks
				scalar td`k'= r(sum)
				
				quietly summ s`k'1  [fw=s`k'5]
				mat Shk[r`k', 4] =(r(sum)/td`k')*100 //percentage of shock from extreme
				
				
				quietly summ s`k'2 [fw=s`k'5]
				mat Shk[r`k', 5] =(r(sum)/td`k')*100 //percentage of shock from moderate
				
				quietly summ s`k'3 [fw=s`k'5]
				mat Shk[r`k', 6] =(r(sum)/td`k')*100 //percentage of shock from mild
						
				}


********Calculating Standard case for 2018**********************
				svmat P7, names(w)
				
				quietly summ w1  [fw=w4], meanonly
				mat VulA[6, 5] =r(sum) //number of individuals who are vuln.
				mat VulA[6, 6]= r(mean) // headcount
				mat VulB[6, 5]= r(mean) // headcount
				
				quietly summ w2 [fw=w4], d
				mat VulA[6, 7] =r(mean) //avg vul
				mat VulA[6, 8] =r(sum)/VulA[6, 5] // intensity 
				*mat VulB[6, 5]= r(p50) // median vulnerability
				
				
				quietly summ w3 [fw=w4], d
				mat VulB[6, 6] =r(mean) //avg vul
				mat VulB[6, 7] =r(sum)/VulA[6, 5] // intensity
				mat VulB[6, 8]= r(p50) // median vulnerability
										
							





scalar drop _all

mat drop P1 P2 P3 P4 P5 P7

mata

 mata drop  P6*  ///

 end
 
restore



**************************Final Vulnerability Result*****************************

*** Vulnerability of 2017 - Table 1
matlist VulA  

*** Contribution of Shocks alpha =1 in 2017/2018 - Table 2
matlist Shk

*** Data Summary By Year Sample- Table C1 Appendix C
matlist HINC

**** Coefficiet of 2017 -  Table C2 Appendix C
matlist PCoeff

**** Coefficient of 2018 - Table C3 Appendix C
matlist PCoeff1

*** Table C4 and C5 is equivalent to Table 1 and Table 2 (above) for AR1 specificiation

**** Data Summary All- Table D1 Appendix D
matlist HINCA

*** Weekcount - Table D2 Appendix D
matlist HWK






