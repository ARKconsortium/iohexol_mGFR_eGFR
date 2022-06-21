************************************************************************** 
  Study ID       : ARK GSK-CKD study
  Program name   : ARKprog.sas
  Description    : derive iGFR results from iohexol measurements alongside
				   eGFR estimates
				   The program accompanies 
  Input          : ARKinput.sas7bdat
  Output         : ARKoutput.sas7bdat
**************************************************************************;

*change the below filepath to point to the data files that you have downloaded onto your drive;
Libname mydata 'C:\ARK study\datarepos';

*Read in iohexol concentration data from Uganda, Malawi and South Africa from drive;
data ARKinput;
	set mydata.ARKinput_final;
run;

*sort dataset by country and sex;
proc sort data=ARKinput out=fitcurve_;
	by country sex;
run;

*import population creatinine values needen in FAS equations;
data fitcurve;
	merge fitcurve_
		  mydata.sumcreat;
	by country sex; 
run;
proc sort data=fitcurve;by idno;run;

*expand dataset into 'long' format: one row per measurement per person;
data long;
	set fitcurve(in=a keep=idno t1_iohexol t1 rename=(t1_iohexol=conc t1=time))
		fitcurve(in=b keep=idno t2_iohexol t2 rename=(t2_iohexol=conc t2=time))
		fitcurve(in=c keep=idno t3_iohexol t3 rename=(t3_iohexol=conc t3=time))
		fitcurve(in=d keep=idno t4_iohexol t4 rename=(t4_iohexol=conc t4=time));
	if a then timepoint = 1;
	if b then timepoint = 2;
	if c then timepoint = 3;
	if d then timepoint = 4;
run;

*sort dataset by participant id and measurement time point;
proc sort data=long;
	by idno time;
run;

*only measurements from time points 2 (2hrs), 3 (3hrs) and 4 (4hrs) are used in the analysis;
data final;
	set long;
	*remove the first measurements from the dataset;
	if timepoint ne 1;

	*log transform all iohexol concentration measurements;
	logy = log(conc);
run;

*If needed, direct SAS output and log to external files in designated output folder;
*change the below filepaths to point to a folder on your drive;
PROC PRINTTO log="C:\ARK study\Output\logoutput.txt" new;
RUN;
PROC PRINTTO print="C:\ARK study\Output\output.lst" new;
RUN;

*for each participant, estimate line of best fit through the three measurements points (time,logy).;
*save the intercept and slope in dataset sols;
proc mixed data=final;
	by idno;
	model logy=time/s;
	ods output SolutionF=sols;
run;

*for each participant, calculate the r-statistic for the line of best fit;
proc corr data=final;
	by idno;
	var logy time;
	ods output PearsonCorr = pearson_exp;
run;

*Calculate R-squared as the square of the Pearson correlation coefficient;
data rsquared_exp;
	set pearson_exp(where=(variable="logy") keep=idno time variable);
	rsquared_exp=time**2;
	keep idno rsquared_exp;
run;

*Direct SAS output back to terminal;
PROC PRINTTO ;
RUN;

*combine estimates of intercept and slope on sigle line and exponentiate to;
*get constant and rate parameter for the exponential curve;
data combine_exp;
	merge	sols(in=a where=(effect="Intercept") keep=idno effect estimate rename=(estimate=lconst))
			sols(in=b where=(effect="time") 	 keep=idno effect estimate rename=(estimate=lrate ));
	by idno;
	drop effect;
	expconst = exp(lconst);
	exprate  = exp(lrate );
run;

proc sort data=fitcurve; by idno; run;

*append exponential curve and R-squared onto original dataset;
data broad;
	merge	fitcurve
			rsquared_exp
			combine_exp;
	by idno;
	r = sqrt(rsquared_exp);
run;

data mGFR;
	set broad;

	*calculate volume of iohexol injected. This is done by converting the weight of the syringe 
	*content (in grams) to mL by multiplying by the specific gravity for iohexol, 1.406;
	ml_injected = (ioh_weight/1.406);

	*calculate total dose (in mg) of iohexol administered. 
	*The iohexol solution injected contains 755 mg of iohexol per mL;
	dose = 755*ml_injected;
	
	*calculate the volume of distribution, Vd;
	Vd 	 = dose / expconst;

	*calculate the slope-intercept GFR;
	SI_GFR = -lrate*Vd*1000;

	*calculate body surface area using the Haycock formula;
	BSA = 0.024265*(weight**0.5378)*height**0.3964;

	*calculate BSA corrected SI_GFR;
	SI_GFR_BSA = SI_GFR*(1.73/BSA);

	*calculate the Brochner-Mortensen GFR from the BSA corrected SI_GFR;
	GFR_BM = 0.990778*SI_GFR_BSA-0.001218*(SI_GFR_BSA**2);

	** Define 2 sensitivity analysis subsets from the full dataset based on, 
	* (1) an r value of at least 0.985, (2) a volume of distribution in the normal range;
	sens_anal1 = (r>=0.985);
	sens_anal2 = (sex = "Male" AND (13<=vd<=20)) OR (sex = "Female" AND (11<=vd<=17));
run;
 
proc sort data=mGFR; by country idno; run;

data eGFR1;
	set mGFR;

	*convert creatinine measured in umol/L to mg/dL;
	creat_mg_dL = creat*0.0113;

	*derive MDRD eGFR estimates;
	sexfactor = (sex EQ "Female")*0.742+(sex NE "Female")*1;
	eGFR_MDRD 		= 175*creat_mg_dL**(-1.154)*age**(-0.203)*sexfactor;

	*apply ethnicity factor of 1.212 ;
	eGFR_MDRD_ethn  = 175*creat_mg_dL**(-1.154)*age**(-0.203)*sexfactor*1.212;

	*derive the creatinine-based 2009 and 2021 CKD-EPI eGFR estimates;
	if sex="Female" then do;
		CKD_EPI_2009   = 141*min(creat/61.9,1)**(-0.329)*max(creat/61.9,1)**(-1.209)*0.9930**age*1.018;
		CKD_EPI_2021   = 142*min(creat/61.9,1)**(-0.241)*max(creat/61.9,1)**(-1.200)*0.9938**age*1.012;
	end;
	if sex="Male" then do;
		CKD_EPI_2009   = 141*min(creat/79.6,1)**(-0.411)*max(creat/79.6,1)**(-1.209)*0.9930**age;
		CKD_EPI_2021   = 142*min(creat/79.6,1)**(-0.302)*max(creat/79.6,1)**(-1.200)*0.9938**age;
	end;

	*apply ethnicity factor to the CKD-EPI 2009 eGFR estimates;
	CKD_EPI_2009_ethn = CKD_EPI_2009*1.159;

	* derive Cockcroft-Gault eGFR estimates with and without BSA correction;
	CG_eGFR = (140-age)*weight* (1.04*(sex EQ "Female") + 1.23*(sex NE "Female")) / creat;
	CGbsa_eGFR = CG_eGFR*(1.73/BSA);

	* derive Lund-Malmo eGFR estimates;
	if sex eq "Female" then do;
		if creat  < 150 then X = 2.50 + 0.0121*(150-creat);
		if creat >= 150 then X = 2.50 - 0.9260*log(creat/150);
	end;
	if sex eq "Male" then do;
		if creat  < 180 then X = 2.56 + 0.00968*(180-creat);
		if creat >= 180 then X = 2.56 - 0.92600*log(creat/180);
	end;
	LM_eGFR = exp(X-0.0158*age+0.438*log(age));
run;

Data eGFR2;
	set eGFR1;

	*derive eGFR estimates from FAS equation;
	Q = popcreat_median;
	FAS = 107.3/(creat/Q);
	if age>40 then FAS = FAS*(0.988**(age-40));
	
	*derive eGFR estimates from the CKD-EPI cystatin and cystatin-creatinine equations;
	if cystatin ne . then do;
		if 		cystatin <= 0.8 then do;
			eGFR_cys 		= 133*((cystatin/0.8)**(-0.499))*0.996**age*(0.932*(sex EQ "Female")+(sex EQ "Male"));
			if sex EQ "Female" then do;
				if creat_mg_dL<=0.7 then do;
					eGFR_cys_creat = 130.815*((creat_mg_dL/0.7)**(-0.248))*((cystatin/0.8)**(-0.375)) *0.995**age;
				end;
				if creat_mg_dL>0.7  then do;
					eGFR_cys_creat = 130.815*((creat_mg_dL/0.7)**(-0.601))*((cystatin/0.8)**(-0.375)) *0.995**age;
				end;
			end; 
			if sex EQ "Male" then do;
				if creat_mg_dL<=0.9 then do;
					eGFR_cys_creat = 135*((creat_mg_dL/0.9)**(-0.207))*((cystatin/0.8)**(-0.375)) *0.995**age;
				end;
				if creat_mg_dL>0.9  then do;
					eGFR_cys_creat = 135*((creat_mg_dL/0.9)**(-0.601))*((cystatin/0.8)**(-0.375)) *0.995**age;
				end;
			end; 		
		end;
		else if cystatin >  0.8 then do;
			eGFR_cys 		= 133*((cystatin/0.8)**(-1.328))*0.996**age*(0.932*(sex EQ "Female")+(sex EQ "Male"));
			if sex EQ "Female" then do;
				if creat_mg_dL<=0.7 then do;
					eGFR_cys_creat = 130.815*((creat_mg_dL/0.7)**(-0.248))*((cystatin/0.8)**(-0.711)) *0.995**age;
				end;
				if creat_mg_dL>0.7  then do;
					eGFR_cys_creat = 130.815*((creat_mg_dL/0.7)**(-0.601))*((cystatin/0.8)**(-0.711)) *0.995**age;
				end;
			end; 
			if sex EQ "Male" then do;
				if creat_mg_dL<=0.9 then do;
					eGFR_cys_creat = 135*((creat_mg_dL/0.9)**(-0.207))*((cystatin/0.8)**(-0.711)) *0.995**age;
				end;
				if creat_mg_dL>0.9  then do;
					eGFR_cys_creat = 135*((creat_mg_dL/0.9)**(-0.601))*((cystatin/0.8)**(-0.711)) *0.995**age;
				end;
			end; 		
		end;
	end;

	*derive eGFR estimates from the 2021 CKD-EPI cystatin-creatinine equation;
	if cystatin ne . then do;
		if 		cystatin <= 0.8 then do;
			if sex EQ "Female" then do;
				if creat_mg_dL<=0.7 then do;
					eGFR_cys_creat_2021 = 130.005*((creat_mg_dL/0.7)**(-0.219))*((cystatin/0.8)**(-0.323)) *0.9961**age;
				end;
				if creat_mg_dL>0.7  then do;
					eGFR_cys_creat_2021 = 130.005*((creat_mg_dL/0.7)**(-0.544))*((cystatin/0.8)**(-0.323)) *0.9961**age;
				end;
			end; 
			if sex EQ "Male" then do;
				if creat_mg_dL<=0.9 then do;
					eGFR_cys_creat_2021 = 135*((creat_mg_dL/0.9)**(-0.144))*((cystatin/0.8)**(-0.323)) *0.9961**age;
				end;
				if creat_mg_dL>0.9  then do;
					eGFR_cys_creat_2021 = 135*((creat_mg_dL/0.9)**(-0.544))*((cystatin/0.8)**(-0.323)) *0.9961**age;
				end;
			end; 		
		end;
		else if cystatin >  0.8 then do;
			if sex EQ "Female" then do;
				if creat_mg_dL<=0.7 then do;
					eGFR_cys_creat_2021 = 130.005*((creat_mg_dL/0.7)**(-0.219))*((cystatin/0.8)**(-0.778)) *0.9961**age;
				end;
				if creat_mg_dL>0.7  then do;
					eGFR_cys_creat_2021 = 130.005*((creat_mg_dL/0.7)**(-0.544))*((cystatin/0.8)**(-0.778)) *0.9961**age;
				end;
			end; 
			if sex EQ "Male" then do;
				if creat_mg_dL<=0.9 then do;
					eGFR_cys_creat_2021 = 135*((creat_mg_dL/0.9)**(-0.144))*((cystatin/0.8)**(-0.778)) *0.9961**age;
				end;
				if creat_mg_dL>0.9  then do;
					eGFR_cys_creat_2021 = 135*((creat_mg_dL/0.9)**(-0.544))*((cystatin/0.8)**(-0.778)) *0.9961**age;
				end;
			end; 		
		end;
	end;
run;

Data eGFR3;
	set eGFR2;

	*derive eGFR estimates from the new FAS equation;
	*step 1: determine Q;
	if age <= 25 then do;
		if sex EQ "Female" then ln_Q = 3.080 + 0.177*age - 0.223*log(age) - (0.00596*(age**2)) + (0.0000686*(age**3));
		if sex EQ "Male"   then ln_Q = 3.200 + 0.259*age - 0.543*log(age) - (0.00763*(age**2)) + (0.0000790*(age**3));
		Q = exp(ln_Q);
	end;
	if age >  25 then do;
		if sex EQ "Female" 	then Q = 62;
		if sex EQ "Male" 	then Q = 80;
	end;

	*step 2: evaluate FAS;
	if age <= 40 then do;
		if ((creat/Q)< 1) then FAS2021 = 107.3*((creat/Q)**-0.322);
		if ((creat/Q)>=1) then FAS2021 = 107.3*((creat/Q)**-1.132);
	end;
	if age >  40 then do;
		if ((creat/Q)< 1) then FAS2021 = 107.3 * ((creat/Q)**-0.322) * (0.990**(age-40));
		if ((creat/Q)>=1) then FAS2021 = 107.3 * ((creat/Q)**-1.132) * (0.990**(age-40));		
	end;
run;

*derive eGFR estimates from the ARK equations;
data eGFR4;
	set eGFR3;
	BMI = weight/(height*height/10000);

	* ARK Model #1: separate BMI coefficient for women and men;
	if sex="Female" then ARKM1 = 126.2797899*(creat_mg_dL/0.818730753)**(-0.3443)*0.993084**age*0.992498279**bmi;
	if sex="Male" 	then ARKM1 = 126.2797899*min(creat_mg_dL/0.818730753,1)**(-0.3443)*max(creat_mg_dL/0.82,1)**(-0.5708)*0.993084**age*0.99911**bmi;

	* ARK Model #2: common BMI coefficient for women and men;
	if sex="Female" then ARKM2 = 120.9769462*(creat_mg_dL/0.818730753)**(-0.3395)*0.993054**age*0.994117**bmi;
	if sex="Male" 	then ARKM2 = 141.6557757*min(creat_mg_dL/0.818730753,1)**(-0.3395)*max(creat_mg_dL/0.82,1)**(-0.5585)*0.993054**age*0.994117**bmi;

	* ARK Model #3: without BMI;
	if sex="Female" then ARKM3 = 102.7192974*(creat_mg_dL/0.818730753)**(-0.3393)*0.993094**age;
	if sex="Male" 	then ARKM3 = 123.7174084*min(creat_mg_dL/0.818730753,1)**(-0.3393)*max(creat_mg_dL/0.82,1)**(-0.5742)*0.993094**age;
run;

*output final dataset to drive;
data mydata.ARKoutput;
	set eGFR4;
run;
