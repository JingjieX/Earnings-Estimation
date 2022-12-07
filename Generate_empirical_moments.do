clear all

set maxvar 32767

use "/Users/jingjiexu/Desktop/509 earnings estimation/Earnings estimation/Ready_newdata/ready_newdata.dta"
ssc install tabout, replace


/*---------------------------------Data Cleaning------------------------------*/
/*------------------Set initial conditions----------------*/
gen Tinit=67
gen Tlast=96
gen ageinit=20
gen agelast=64
gen minyrs=20



/*--------IMPUTING MISSING GRADES VARIABLE IN THE YEARS 70-74-------------*/
local i=94
while `i'<=97 {
     rename upedu`i'h grade`i'
     local i=`i'+1
}

replace grade72=edcn72
replace grade72=. if grade72>25
replace grade75=edcn75
replace grade75=. if grade75>25

gen grade69=0
replace grade69=grade68 if seqno68==1 & seqno69==1
replace grade69=grade72 if grade69==0 & seqno69==1 & seqno72==1

replace grade70=grade68 if seqno68==1 & seqno70==1
replace grade70=grade72 if grade69==0 & seqno70==1 & seqno72==1

replace grade71=grade70 if seqno70==1 & seqno71==1
replace grade71=grade72 if grade71==0 & seqno71==1 & seqno72==1

replace grade73=grade72 if seqno72==1 & seqno73==1
replace grade73=grade75 if grade73==0 & seqno73==1 & seqno75==1

replace grade74=grade72 if seqno72==1 & seqno74==1
replace grade74=grade75 if grade74==0 & seqno74==1 & seqno75==1


/*---------REMOVING MISSING OBSERVATIONS IN GRADES-----------*/
local i=Tinit+1
while `i'<=Tlast+1 {
     replace grade`i'=. if grade`i'>30
     local i=`i'+1
}


/*---------Input average nominal wage rates-----------*/
gen awg67=2.85
gen awg68=3.02
gen awg69=3.22
gen awg70=3.40
gen awg71=3.63
gen awg72=3.90
gen awg73=4.14
gen awg74=4.43
gen awg75=4.73
gen awg76=5.06
gen awg77=5.44
gen awg78=5.87
gen awg79=6.33
gen awg80=6.84
gen awg81=7.43
gen awg82=7.86
gen awg83=8.19
gen awg84=8.48
gen awg85=8.73
gen awg86=8.92
gen awg87=9.13
gen awg88=9.43
gen awg89=9.80
gen awg90=10.19
gen awg91=10.50
gen awg92=10.76
gen awg93=11.03
gen awg94=11.32
gen awg95=11.64
gen awg96=12.03

/*Note: awg96==rawg96 by definition of 96 as base year*/



/*--------Generating dummy variables for families with head aged AGEINIT to AGELAST, with male head of family-----------*/
/*This will have at least "minyrs" years of datapoint for each surviving person, but not necessarily consecutive 20 years.*/


local i = Tinit-1
while `i' < Tlast {
       local i=`i'+1
       gen rawg`i'=awg`i'/prc`i'
       replace rhdlbin`i'=. if rhdwg`i'<=(2*rawg`i'/awg96) | rhdwg`i'>(400*rawg`i'/awg96) | hwkhrs`i'>5110 | hwkhrs`i'<520
       replace rhdlbin`i'=. if rhdlbin`i'==0 & hwkhrs`i'>0
       replace rhdlbin`i'=. if rhdlbin`i'>0 & hwkhrs`i'==0
}  





gen kept=0
local i = Tinit /*+1*/
while `i'<= Tlast {
      local i = `i'+1
      local ii = `i'-1
      /*-----------Identify people in the age groups--------*/
      gen agedum`i'=0
      replace agedum`i'=1 if agehd`i'>=ageinit & agehd`i'<=agelast
      /*--------- Identify head of households------------- */
      gen seqdum`i'=0
      replace seqdum`i'=1 if seqno`i'==1
      /*------------------Identify sex of the head of household---------- */
      gen sexdum`i'=0
      replace sexdum`i'=1 if sexhd`i'==1
      /* -------------Identify if labor income is positive in that year------------*/
      gen labdum`ii'=0
      replace labdum`ii'=1 if rhdlbin`ii'>0 & rhdlbin`ii'~=.

       replace kept=kept+agedum`i'*seqdum`i'*sexdum`i'*labdum`ii'

      gen lnrhdlbin`ii'=ln(rhdlbin`ii')
      gen expr`i'=agehd`i'-max(grade`i',12)-6

}

/*----------------KEEP PEOPLE WITH AT LEAST "MINYRS" YEARS OF DATA POINTS -------------*/

drop if kept < minyrs

/*-------------DROP ANY REMAINING UNUSABLE VARIABLES----------------*/
drop educ* prc* rawg* awg* upedu*
aorder

/*-----------The sequence no. of each observation--------*/
gen idno = _n
/*-----------Number of observations left--------*/
gen k = _N 
display k

/* -------------CONVERT DATA INTO LONG FORM FOR PANEL IMPLEMENTATION-----------*/
reshape long age agehd agedum expr hdedcn hdlbin hdwg rhdwg grade hwkhrs id labdum rhdlbin lnrhdlbin numfam relh seqdum seqno sexdum sexhd, i(idno) j(year)

drop if year < Tinit | year==67
drop if year > Tlast

tsset idno year



save "/Users/jingjiexu/Desktop/509 earnings estimation/Earnings estimation/ready_newdata_cleaned.dta"







/*---------------A Time-invariant model specification-------------------------*/



use "/Users/jingjiexu/Desktop/509 earnings estimation/Earnings estimation/ready_newdata_cleaned.dta"



/*---------------GENERATE RESIDUALS-------------------------*/


gen residual=.
local j=ageinit
	 while `j'<= agelast{
		disp "agehd=`j'"
		quietly reg lnrhdlbin agehd expr if seqno*labdum*agedum*sexdum==1 & agehd==`j'
		predict resid`j', resid
		replace residual=resid`j' if seqno*labdum*agedum*sexdum==1 & agehd==`j'		
	 local j = `j'+1
}



/*------------------OBTAIN SAMPLE MOMENTS--------------*/

local j=ageinit
	while `j'<= agelast{
		local n=0
		while `n'<= agelast -`j'{
			
*			local `j1'=`j'+1

			local jn=`j'+`n'  //`jn' is used to label y_i,j+n as yi`jn'
			
*			gen yi`j'_yi`j1'=.
*			replace yi`j'_yi`j1'=resid`j'*resid`j1'    //y_ij*y_i,j+1


			gen yi`j'_yi`jn'=.
			replace yi`j'_yi`jn'=resid`j'*resid`jn'  //the product y_j*y_j+n
						
*			egen m_`j'_0=mean(yi`j'_yi`j') 
*			egen m_`j'_1=mean(yi`j'_yi`j1')  
			
			egen m`j'_`jn'=mean(yi`j'_yi`jn')  /*m`j'_`jn': the covariance of residual earnings between age j and j+n individuals*/
			
			if `n'!=0{
				gen m`jn'_`j'=m`j'_`jn'
			}
			
			local n=`n'+1
		}		

		local j=`j'+1
		
}


keep m*_* 
duplicates drop
*outfile using empiricalmoments, nolabel replace wide
export excel using "/Users/jingjiexu/Desktop/509 earnings estimation/Earnings estimation/empiricalmoments.xlsx", sheet("empiricalmoments") firstrow(variables)
save "/Users/jingjiexu/Desktop/509 earnings estimation/Earnings estimation/empiricalmoments.dta"




*Next, turn to matlab to construct matric M and M_hat
		
		
/*		
/*----------CONSTRUCT A MATRIX CONTAINING THEORETICAL MOMENTS----------------*/
	
/*------------------OBTAIN THEORETICAL MOMENTS--------------*/

gen rho=.
gen var_a=.
gen var_e=.
gen var_v=.



local j=ageinit
	while `j'<= agelast{
		local n=0
		while `n'<= agelast -`j'{
			
			local jn=`j'+`n'  
			gen h_m`j'_`jn'=.
			
			if `n'==0 {
				replace h_m`j'_`jn'=var_a+var_e+var_v
			}
			
			else if `n'~=0 {
				replace h_m`j'_`jn'=var_a+rho^n*var_v
				gen h_m`jn'_`j'=h_m`j'_`jn'
			}
			
			local n=`n'+1
		}		

		local j=`j'+1
		
}
*
*/




/*
/*----------SAVE EACH ROW AS A MATRIX-----------*/			
mkmat m20_*, matrix(R20) 
mkmat m21_*, matrix(R21) 
mkmat m22_*, matrix(R22)
mkmat m23_*, matrix(R23)
mkmat m24_*, matrix(R24)
mkmat m25_*, matrix(R25)
mkmat m26_*, matrix(R26)
mkmat m27_*, matrix(R27)
mkmat m28_*, matrix(R28)
mkmat m29_*, matrix(R29)

mkmat m30_*, matrix(R30) 
mkmat m31_*, matrix(R31) 
mkmat m32_*, matrix(R32) 
mkmat m33_*, matrix(R33)
mkmat m34_*, matrix(R34)
mkmat m35_*, matrix(R35)
mkmat m36_*, matrix(R36)
mkmat m37_*, matrix(R37)
mkmat m38_*, matrix(R38)
mkmat m39_*, matrix(R39)

mkmat m40_*, matrix(R40) 
mkmat m41_*, matrix(R41) 
mkmat m42_*, matrix(R42) 
mkmat m43_*, matrix(R43)
mkmat m44_*, matrix(R44)
mkmat m45_*, matrix(R45)
mkmat m46_*, matrix(R46)
mkmat m47_*, matrix(R47)
mkmat m48_*, matrix(R48)
mkmat m49_*, matrix(R49)

mkmat m50_*, matrix(R50) 
mkmat m51_*, matrix(R51) 
mkmat m52_*, matrix(R52) 
mkmat m53_*, matrix(R53)
mkmat m54_*, matrix(R54)
mkmat m55_*, matrix(R55)
mkmat m56_*, matrix(R56)
mkmat m57_*, matrix(R57)
mkmat m58_*, matrix(R58)
mkmat m59_*, matrix(R59)

mkmat m60_*, matrix(R60) 
mkmat m61_*, matrix(R61) 
mkmat m62_*, matrix(R62) 
mkmat m63_*, matrix(R63)
mkmat m64_*, matrix(R64)

*/



/*

/*----------CONVERT THE EMPIRICAL MOMENTS TO MATRIX-----------*/

/*----------SAVE EACH ROW AS A MATRIX-----------*/
local j=20
while (`j'<= 64) {
	mkmat m`j'_*, matrix(R`j')
		
	local `j'=`j'+1	
}
* Now we have 45 （1*65）matrices R20 R21 R22 …… R64. 
* Eg: R20 = [m20_20, m20_21, m20_22, m20_23 …… m20_64]


/*------Horitontally conbine the matrices to obtain the matrix containing empirical moments----*/
*mat M=[R20;R21;R22……;R64]

/*----------VECTORIZE THE MATRIX CONTAINING EMPIRICAL MOMENTS----------*/
mat M = vec(M)

*For simplicity, the weighting matrix W here is set to the identity matrix.
matrix W=I(65)


/*----------CONVERT THE THEORETICAL MOMENTS TO MATRIX-----------*/	

*/

* A 4*4 sample


*Simplicity, here consider two 4*4 matrics cantaining empirical and theoretical moments to solve the problem



/*----------A SUBMATRIX CONTAINING EMPIRICAL MOMENTS-----------*/
keep m20_20 m20_21 m20_22 m20_23 m21_20 m21_21 m21_22 m21_23 m22_20 m22_21 m22_22 m22_23 m23_20 m23_21 m23_22 m23_23

/*----------SAVE EACH ROW AS A MATRIX-----------*/
mkmat m20_*, matrix(R20) 
mkmat m21_*, matrix(R21) 
mkmat m22_*, matrix(R22)
mkmat m23_*, matrix(R23)

/*------Horitontally combine the matrices to obtain the matrix containing empirical moments----*/
matrix M = [R20\R21\R22\R23]
matrix list M

svmat M
outfile using M, nolabel replace wide

/*----------VECTORIZE THE MATRIX CONTAINING EMPIRICAL MOMENTS----------*/
mat M_vectorized = vec(M)
mat l M_vectorized

svmat M_vectorized
outfile using M_vectorized, nolabel replace wide

keep M_vectorized

export excel using "/Users/jingjiexu/Desktop/509 earnings estimation/Earnings estimation/M_vectorized.xlsx", sheet("M_vectorized") firstrow(variables)

















/*---------------A Time-variant model specification---------------------------*/
use "/Users/jingjiexu/Desktop/509 earnings estimation/Earnings estimation/ready_newdata_cleaned.dta", clear


gen agehdsq=(expr^2)/100
gen agehdcu=(expr^3)/1000
gen agehdqr=(expr^4)/10000


gen residual=.
local i=Tinit+1
      while `i'<= Tlast{
          disp "year=`i'"
          quietly reg lnrhdlbin agehd agehdsq agehdcu if seqno*labdum*agedum*sexdum==1 & year==`i'
          replace alphage=_b[agehd] if  seqno*labdum*agedum*sexdum==1 & year==`i'
          replace alphagesq=_b[agehdsq] if  seqno*labdum*agedum*sexdum==1 & year==`i'
          replace alphagecu=_b[agehdcu] if  seqno*labdum*agedum*sexdum==1 & year==`i'
          replace alphacons=_b[_cons] if  seqno*labdum*agedum*sexdum==1 & year==`i'
       *  replace alphanumfam=_b[numfam] if  seqno*labdum*agedum*sexdum==1 & year==`i'
          predict resid`i', resid
          replace residual=resid`i' if seqno*labdum*agedum*sexdum==1 & year==`i'
          drop resid`i'
      local i = `i'+1
}



