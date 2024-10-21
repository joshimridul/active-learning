
/*==========================================================================

Title: create_aggregate_effects.do 
Author: Mridul Joshi
Date: Thu Oct 17 15:13:24 2024

Description: Create aggregate effects from data extracted from studies

===========================================================================*/

clear all
version 16.0
set more off
pause off
cap log close 

local icc 0.20

* Adair et al

/*

hedges g - 0.306
se - 0.149
se_corrected - 0.647

*/

input q m1 sd1 n1 m2 sd2 n2

1 	0.03	0.13	91	0.03	0.14	90
2 	0.26	0.15	91	0.24	0.15	90
3 	0.36	0.16	91	0.22	0.15	90
4 	0.44	0.19	91	0.26	0.18	90
5 	0.49	0.17	91	0.34	0.23	90
end

egen wgt_mean1 = mean(m1)  // all questions answered by all respondents
egen wgt_mean2 = mean(m2)

egen tot_n1 = total(n1)
egen tot_n2 = total(n1)

gen within_term_1 = (sd1^2)*(n1 - 1)
gen within_term_2 = (sd2^2)*(n2 - 1)

egen ss_within_1 = total(within_term_1)
egen ss_within_2 = total(within_term_2)

gen m1_sq = m1*m1
gen m2_sq = m2*m2

egen tot_m1_sq = total(m1_sq)
egen tot_m2_sq = total(m2_sq)

gen ss_between_1 = n1*tot_m1_sq - ((n1^2)*tot_m1_sq)/tot_n1
gen ss_between_2 = n1*tot_m2_sq - ((n2^2)*tot_m2_sq)/tot_n2

gen overall_sd_1 = sqrt((ss_between_1 + ss_within_1)/(tot_n1 - 1))

gen overall_sd_2 = sqrt((ss_between_2 + ss_within_2)/(tot_n2 - 1))

gen pooled_sd =  sqrt(((n1 - 1)*overall_sd_1^2 + (n2 - 1)*overall_sd_2^2)/(n1+n2-1))

gen cohen_d = (wgt_mean1 - wgt_mean2)/pooled_sd

gen var_d =  (n1 + n2)/(n1*n2) + ( cohen_d^2/(2*(n1+n2)) )

gen se_d = sqrt(var_d)

gen j = 1 - (3/(4*(n1+n2-2)-1) )

gen hedges_g = cohen_d * (j)

gen se_g = se_d*j

gen deff = 1 + (((n1 + n2)/2) - 1 )*`icc'

gen se_d_corr = se_d*sqrt(deff)

gen se_g_corr = se_g*sqrt(deff)

gen t_stat_d = cohen_d/se_d

gen t_stat_d_corr = cohen_d/se_d_corr

gen t_stat_g = hedges_g/se_g

gen t_stat_g_corr = hedges_g/se_g_corr

gen sig_d = abs(t_stat_d_corr) > 1.96

gen sig_g = abs(t_stat_g_corr) > 1.96


* Bowman et al 

/*

hedges g -  -0.119
se - 0.051
se_corrected - 0.051

*/


clear 

input coefficient std_err n_t n_c n_t_clus n_c_clus
	  -0.133	0.057	 796	 734	49	 47
end


gen pooled_sd = std_err /sqrt((1/n_t) + (1/n_c))

gen p = `icc'
gen n = n_t + n_c
gen m = n_t_clus + n_c_clus
gen n_by_m = n/m

* Calculate h

gen h_num =  ( (n-2) - 2*(n_by_m -1)*p )^2

gen h_den = (n-2)*(1-p)^2 + n_by_m*(n - 2*n_by_m)*p^2 + 2*(n - 2*n_by_m)*p*(1-p)

gen h = h_num/h_den

gen j = 1 - ( 3 /(4*h -1))

gen lambda = 1 - (2*n_by_m*p/(n-1))

gen cohen_d = (coefficient/pooled_sd)*sqrt(lambda)

gen hedges_g = j*cohen_d*sqrt(lambda)


gen se_d = sqrt(lambda*(std_err/pooled_sd)^2 + (cohen_d^2/(2*h)) )

gen se_g = sqrt((j^2)*lambda*(std_err/pooled_sd)^2 + (hedges_g^2/(2*h)) )


gen t_stat_d = cohen_d/se_d

gen t_stat_g = hedges_g/se_g

gen sig_d = abs(t_stat_d) > 1.96

gen sig_g = abs(t_stat_g) > 1.96


* Chan and Bauer

/*

hedges g -  0.0098
se - 0.094
se_corrected - 0.454

*/


clear 

input year exam score_c score_t sd n
	  2008 1 66 68 15 297
	  2008 2 61 62 15 297
	  2008 3 58 60 15 297
	  2008 4 63 64 15 297
	  2009 1 72 72 18 150
	  2009 2 75 72 18 150
	  2009 3 65 61 18 150
	  2009 4 64 61 18 150
end

collapse score_c score_t sd n, by(year)

/*
year	score_c	score_t	sd	n
2008	62		63.5	15	297
2009	69		66.5	18	150
*/

gen m_c =  ((62*297) + (69*150))/(297+150)
gen m_t =  ((63.5*297) + (66.5*150))/(297+150)

gen within_term = (sd^2)*(n-1)
egen ss_within = total(within_term)

gen between_term1_c = (m_c^2)*n 
egen sum_between_term1_c = total(between_term1_c)

gen between_term1_t = (m_t^2)*n 
egen sum_between_term1_t = total(between_term1_t)


gen between_term2_c = (m_c)*n 
egen sum_between_term2_c = total(between_term2_c)

gen between_term2_t = (m_t)*n 
egen sum_between_term2_t = total(between_term2_t)

gen sum_between_term2_c_sq = sum_between_term2_c*sum_between_term2_c
gen sum_between_term2_t_sq = sum_between_term2_t*sum_between_term2_t

gen ss_between_c = sum_between_term1_c - sum_between_term2_c_sq/(297+150)

gen ss_between_t = sum_between_term1_t - sum_between_term2_t_sq/(297+150)

gen sd_c = sqrt((ss_within + ss_between_c)/(297+150-1))
gen sd_t = sqrt((ss_within + ss_between_t)/(297+150-1))


gen pooled_sd =  sqrt(((223.5 - 1)*sd_c^2 + (223.5 - 1)*sd_t^2)/(297+150-1))

gen cohen_d = (m_t - m_c)/pooled_sd

gen var_d =  (223.5 + 223.5)/(223.5*223.5) + ( cohen_d^2/(2*(223.5+223.5)) )

gen se_d = sqrt(var_d)

gen j = 1 - (3/(4*(223.5 + 223.5 -2)-1) )

gen hedges_g = cohen_d * (j)

gen se_g = se_d*j

gen deff = 1 + (111.75 - 1)*`icc'

gen se_d_corr = se_d*sqrt(deff)

gen se_g_corr = se_g*sqrt(deff)

gen t_stat_d = cohen_d/se_d

gen t_stat_d_corr = cohen_d/se_d_corr

gen t_stat_g = hedges_g/se_g

gen t_stat_g_corr = hedges_g/se_g_corr

gen sig_d = abs(t_stat_d_corr) > 1.96

gen sig_g = abs(t_stat_g_corr) > 1.96


* Gruner 

/*

hedges g -  -0.0009
se - 0.111
se_corrected - 0.36

*/


clear

input period m_c sd_c n_c m_t sd_t n_t
		  1 1179.74 99.1 16 1161.75 99.1 55
		  2 1119.59 97 99 1118.61 97 220
end


gen mn_c = m_c*n_c
gen mn_t = m_t*n_t

egen tot_mn_t = total(mn_t)
egen tot_mn_c = total(mn_c)


egen tot_n_t = total(n_t)
egen tot_n_c = total(n_c)

gen overall_m_c = tot_mn_c/tot_n_c
gen overall_m_t = tot_mn_t/tot_n_t

gen within_term_c = (sd_c^2)*(n_c-1)
egen ss_within_c = total(within_term_c)

gen within_term_t = (sd_t^2)*(n_t-1)
egen ss_within_t = total(within_term_t)

gen between_term1_c = (m_c^2)*n_c 
egen sum_between_term1_c = total(between_term1_c)

gen between_term1_t = (m_t^2)*n_t 
egen sum_between_term1_t = total(between_term1_t)


gen between_term2_c = (m_c)*n_c 
egen sum_between_term2_c = total(between_term2_c)

gen between_term2_t = (m_t)*n_c 
egen sum_between_term2_t = total(between_term2_t)

gen sum_between_term2_c_sq = sum_between_term2_c*sum_between_term2_c
gen sum_between_term2_t_sq = sum_between_term2_t*sum_between_term2_t

gen ss_between_c = sum_between_term1_c - sum_between_term2_c_sq/(99+16)

gen ss_between_t = sum_between_term1_t - sum_between_term2_t_sq/(220+55)

gen overall_sd_c = sqrt((ss_within_c + ss_between_c)/(390-1))
gen overall_sd_t = sqrt((ss_within_t + ss_between_t)/(390-1))


gen pooled_sd =  sqrt(((99 + 16 - 1)*overall_sd_c^2 + (220 + 55 - 1)*overall_sd_t^2)/(390-1))

gen cohen_d = (overall_m_t - overall_m_c)/pooled_sd

gen var_d =  (390)/(275*115) + ( cohen_d^2/(2*(390)) )

gen se_d = sqrt(var_d)

gen j = 1 - (3/(4*(390 -2)-1) )

gen hedges_g = cohen_d * (j)

gen se_g = se_d*j

gen deff = 1 + (48.75 - 1)*`icc'

gen se_d_corr = se_d*sqrt(deff)

gen se_g_corr = se_g*sqrt(deff)

gen t_stat_d = cohen_d/se_d

gen t_stat_d_corr = cohen_d/se_d_corr

gen t_stat_g = hedges_g/se_g

gen t_stat_g_corr = hedges_g/se_g_corr

gen sig_d = abs(t_stat_d_corr) > 1.96

gen sig_g = abs(t_stat_g_corr) > 1.96


* Kazemi 

/*

hedges g -  0.133
se - 0.218
se_corrected - 0.657

*/


clear 

input m_c sd_c n_c m_t sd_t n_t
	  15.07 8.25 42 16.11 7.31 41			 	
end

gen overall_m_t = m_t
gen overall_m_c = m_c 

gen overall_sd_c = sd_c
gen overall_sd_t = sd_t

gen pooled_sd =  sqrt(((42 - 1)*overall_sd_c^2 + (41- 1)*overall_sd_t^2)/(83-1))

gen cohen_d = (overall_m_t - overall_m_c)/pooled_sd

gen var_d =  (83)/(42*41) + ( cohen_d^2/(2*(83)) )

gen se_d = sqrt(var_d)

gen j = 1 - (3/(4*(83-2)-1) )

gen hedges_g = cohen_d * (j)

gen se_g = se_d*j

gen deff = 1 + (41.5 - 1)*`icc'

gen se_d_corr = se_d*sqrt(deff)

gen se_g_corr = se_g*sqrt(deff)

gen t_stat_d = cohen_d/se_d

gen t_stat_d_corr = cohen_d/se_d_corr

gen t_stat_g = hedges_g/se_g

gen t_stat_g_corr = hedges_g/se_g_corr

gen sig_d = abs(t_stat_d_corr) > 1.96

gen sig_g = abs(t_stat_g_corr) > 1.96


* Kramer

/*

hedges g - 0.771
se - 0.154
se_corrected - 0.154

*/

clear 

input cohen_d ci_l ci_u n    n_clus
	  0.771 0.468 1.073 8111 32	 	   // supplementary material
end

gen se_d = (0.771 -0.468)/1.96

gen j = 1 - (3/(4*(811 -2)-1) )

gen hedges_g = cohen_d * (j)

gen se_g = se_d*j

gen se_d_corr = se_d

gen se_g_corr = se_g

gen t_stat_d = cohen_d/se_d

gen t_stat_d_corr = cohen_d/se_d_corr

gen t_stat_g = hedges_g/se_g

gen t_stat_g_corr = hedges_g/se_g_corr

gen sig_d = abs(t_stat_d_corr) > 1.96

gen sig_g = abs(t_stat_g_corr) > 1.96


* Lopez-Fernandez 

/*

hedges g -  -0.182
se - 0.179
se_corrected - 0.473

*/


clear 

input case m0_c sd0_c m1_c sd1_c n_c  m0_t sd0_t m1_t sd1_t n_t
		1  5.2 1.4 7.4 1.3 37 4.8 1.1 6.9 1.2 38
		2  2.4 1.0 7.3 2.1 24 2.6 1.3 7.1 2.1 25 
end

gen mn0_c = m0_c*n_c
gen mn0_t = m0_t*n_t

gen mn1_c = m1_c*n_c
gen mn1_t = m1_t*n_t

egen tot_mn0_t = total(mn0_t)
egen tot_mn0_c = total(mn0_c)

egen tot_mn1_t = total(mn1_t)
egen tot_mn1_c = total(mn1_c)

egen tot_n0_t = total(n_t)
egen tot_n0_c = total(n_c)

egen tot_n1_t = total(n_t)
egen tot_n1_c = total(n_c)

gen overall_m0_c = tot_mn0_c/tot_n0_c
gen overall_m0_t = tot_mn0_t/tot_n0_t

gen overall_m1_c = tot_mn1_c/tot_n1_c
gen overall_m1_t = tot_mn1_t/tot_n1_t

gen within_term0_c = (sd0_c^2)*(n_c-1)
egen ss_within0_c = total(within_term0_c)

gen within_term1_c = (sd1_c^2)*(n_c-1)
egen ss_within1_c = total(within_term1_c)

gen within_term0_t = (sd0_t^2)*(n_t-1)
egen ss_within0_t = total(within_term0_t)

gen within_term1_t = (sd1_t^2)*(n_t-1)
egen ss_within1_t = total(within_term1_t)


gen between_term1_0_c = (m0_c^2)*n_c 
egen sum_between_term1_0_c = total(between_term1_0_c)

gen between_term1_1_c = (m1_c^2)*n_c 
egen sum_between_term1_1_c = total(between_term1_1_c)

gen between_term1_0_t = (m0_t^2)*n_t 
egen sum_between_term1_0_t = total(between_term1_0_t)

gen between_term1_1_t = (m1_t^2)*n_t 
egen sum_between_term1_1_t = total(between_term1_1_t)


gen between_term2_0_c = (m0_c)*n_c 
egen sum_between_term2_0_c = total(between_term2_0_c)

gen between_term2_1_c = (m1_c)*n_c 
egen sum_between_term2_1_c = total(between_term2_1_c)

gen between_term2_0_t = (m0_t)*n_t 
egen sum_between_term2_0_t = total(between_term2_0_t)

gen between_term2_1_t = (m1_t)*n_t 
egen sum_between_term2_1_t = total(between_term2_1_t)


gen sum_between_term2_0_c_sq = sum_between_term2_0_c*sum_between_term2_0_c
gen sum_between_term2_1_c_sq = sum_between_term2_1_c*sum_between_term2_1_c

gen sum_between_term2_0_t_sq = sum_between_term2_0_t*sum_between_term2_0_t
gen sum_between_term2_1_t_sq = sum_between_term2_1_t*sum_between_term2_1_t


gen ss_between0_c = sum_between_term1_0_c - sum_between_term2_0_c_sq/(37+24)
gen ss_between1_c = sum_between_term1_1_c - sum_between_term2_1_c_sq/(37+24)

gen ss_between0_t = sum_between_term1_0_t - sum_between_term2_0_t_sq/(38+25)
gen ss_between1_t = sum_between_term1_1_t - sum_between_term2_1_t_sq/(38+25)

gen overall_sd0_c = sqrt((ss_within0_c + ss_between0_c)/(124-1))
gen overall_sd0_t = sqrt((ss_within0_t + ss_between0_t)/(124-1))

gen overall_sd1_c = sqrt((ss_within1_c + ss_between1_c)/(124-1))
gen overall_sd1_t = sqrt((ss_within1_t + ss_between1_t)/(124-1))


gen pooled_sd = sqrt(((overall_sd1_t^2)*(n_t-1) + (overall_sd1_c^2)*(n_c-1)) / (n_c + n_t - 2))

gen cohen_d = ((overall_m1_t - overall_m0_t) -  (overall_m1_c - overall_m0_c))/pooled_sd

gen var_d =  (124)/(61*63) + ( cohen_d^2/(2*(124)) )

gen se_d = sqrt(var_d)

gen j = 1 - (3/(4*(124-2)-1) )

gen hedges_g = cohen_d * (j)

gen se_g = se_d*j

gen deff = 1 + (31 - 1)*`icc'

gen se_d_corr = se_d*sqrt(deff)

gen se_g_corr = se_g*sqrt(deff)

gen t_stat_d = cohen_d/se_d

gen t_stat_d_corr = cohen_d/se_d_corr

gen t_stat_g = hedges_g/se_g

gen t_stat_g_corr = hedges_g/se_g_corr

gen sig_d = abs(t_stat_d_corr) > 1.96

gen sig_g = abs(t_stat_g_corr) > 1.96



* Lovelace 

/*

hedges g -  0.634
se - 0.330
se_corrected - 0.700

*/

clear 

input tstat n_c n_t
	  1.97 19 18	   // supplementary material
end

gen cohen_d = (2*tstat)/sqrt(n_c + n_t)

gen var_d =  (37)/(19*18) + ( cohen_d^2/(2*(37)) )

gen se_d = sqrt(var_d)

gen j = 1 - (3/(4*(37-2)-1) )

gen hedges_g = cohen_d * (j)

gen se_g = se_d*j

gen deff = 1 + (18.5 - 1)*`icc'

gen se_d_corr = se_d*sqrt(deff)

gen se_g_corr = se_g*sqrt(deff)

gen t_stat_d = cohen_d/se_d

gen t_stat_d_corr = cohen_d/se_d_corr

gen t_stat_g = hedges_g/se_g

gen t_stat_g_corr = hedges_g/se_g_corr

gen sig_d = abs(t_stat_d_corr) > 1.96

gen sig_g = abs(t_stat_g_corr) > 1.96

 


* Miller


/*

hedges g -  0.178
se - 0.093
se_corrected - 0.547

*/


clear 

input m_c sd_c n_c m_t sd_t n_t
	  112.20 48.04 177 120.66 47.28 320		 	
end



gen overall_m_t = m_t
gen overall_m_c = m_c 

gen overall_sd_c = sd_c
gen overall_sd_t = sd_t

gen pooled_sd =  sqrt(((n_c - 1)*overall_sd_c^2 + (n_t- 1)*overall_sd_t^2)/(n_c + n_t-1))

gen cohen_d = (overall_m_t - overall_m_c)/pooled_sd

gen var_d =  (n_c + n_t)/(n_c*n_t) + ( cohen_d^2/(2*(n_c+n_t)) )

gen se_d = sqrt(var_d)

gen j = 1 - (3/(4*(n_c+n_t-2)-1) )

gen hedges_g = cohen_d * (j)

gen se_g = se_d*j

gen deff = 1 + (165.66 - 1)*`icc'

gen se_d_corr = se_d*sqrt(deff)

gen se_g_corr = se_g*sqrt(deff)

gen t_stat_d = cohen_d/se_d

gen t_stat_d_corr = cohen_d/se_d_corr

gen t_stat_g = hedges_g/se_g

gen t_stat_g_corr = hedges_g/se_g_corr

gen sig_d = abs(t_stat_d_corr) > 1.96

gen sig_g = abs(t_stat_g_corr) > 1.96


* Pando-cerra 

/*

hedges g -  0.4697
se - 0.2057
se_corrected - 0.663

*/


clear 

input m_c sd_c n_c m_t sd_t n_t
	  56.95 20.64 45 66.38 19.46 51 	 	
end



gen overall_m_t = m_t
gen overall_m_c = m_c 

gen overall_sd_c = sd_c
gen overall_sd_t = sd_t

gen pooled_sd =  sqrt(((n_c - 1)*overall_sd_c^2 + (n_t- 1)*overall_sd_t^2)/(n_c + n_t-1))

gen cohen_d = (overall_m_t - overall_m_c)/pooled_sd

gen var_d =  (n_c + n_t)/(n_c*n_t) + ( cohen_d^2/(2*(n_c+n_t)) )

gen se_d = sqrt(var_d)

gen j = 1 - (3/(4*(n_c+n_t-2)-1) )

gen hedges_g = cohen_d * (j)

gen se_g = se_d*j

gen deff = 1 + (48 - 1)*`icc'

gen se_d_corr = se_d*sqrt(deff)

gen se_g_corr = se_g*sqrt(deff)

gen t_stat_d = cohen_d/se_d

gen t_stat_d_corr = cohen_d/se_d_corr

gen t_stat_g = hedges_g/se_g

gen t_stat_g_corr = hedges_g/se_g_corr

gen sig_d = abs(t_stat_d_corr) > 1.96

gen sig_g = abs(t_stat_g_corr) > 1.96


* Roselli 

/*

hedges g -  0.1288
se - 0.2127
se_corrected - 0.6557

*/

clear

input cohen_d n_c n_t
	  0.13 44 43
end

gen var_d =  (n_c + n_t)/(n_c*n_t) + ( cohen_d^2/(2*(n_c+n_t)) )

gen se_d = sqrt(var_d)

gen j = 1 - (3/(4*(n_c+n_t-2)-1) )

gen hedges_g = cohen_d * (j)

gen se_g = se_d*j

gen deff = 1 + (43.5 - 1)*`icc'

gen se_d_corr = se_d*sqrt(deff)

gen se_g_corr = se_g*sqrt(deff)

gen t_stat_d = cohen_d/se_d

gen t_stat_d_corr = cohen_d/se_d_corr

gen t_stat_g = hedges_g/se_g

gen t_stat_g_corr = hedges_g/se_g_corr

gen sig_d = abs(t_stat_d_corr) > 1.96

gen sig_g = abs(t_stat_g_corr) > 1.96




