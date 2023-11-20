
*Manually creating instruments for each panel and product (diesel j=1, gasoline j=2) and region (indexed by "i"). 

estimates clear		

forvalues j=1/2{
drop p`j'
}																													

quietly reg promedio_bandera  max_comp_REGION  premium tarjeta_puntos  competencia_mayorista_provincia  eess_prop_tot_prov urbanizacion  ln_parque_automotor if (anio!=2015 & identificador==1 & utilidad_media>0 & canal_reclasificado!="EESS propias" & gasoleos==1 & max_comp_REGION>0 & max_comp_REGION!=missing(max_comp_REGION)), vce(robust)
predict p1, xb	
replace p1=0 if missing(p1)																														

quietly reg promedio_bandera promedio_otros_canales refineria premium tarjeta_puntos competencia_mayorista_provincia  eess_prop_tot_prov urbanizacion ln_parque_automotor if (anio!=2015 & identificador==1 & utilidad_media>0 & canal_reclasificado!="EESS propias" & naftas==1 & promedio_otros_canales>0), vce(robust)
predict p2, xb																															
replace p2=0 if missing(p2)	


matrix drop _all
set matsize 9000

*Creating one group for each product (diesel j=1, gasoline j=2) and region (indexed by "i"):

forvalues j=1/2{
forvalues i=1/5{
drop grupo_`j'`i'
gen grupo_`j'`i'=0
}
}

*Filter observations by taking only those with available instrument (also leaving out those obs. that do not belong to the panel):
forvalues i=1/5{
replace grupo_1`i'=1 if (identificador==1 & utilidad_media>0 & gasoleos==1 & REGIONES==`i' & anio!=2015  & canal_reclasificado!="EESS propias"  & p1>0)
replace grupo_2`i'=1 if (identificador==1 & utilidad_media>0 & naftas==1 & REGIONES==`i' & anio!=2015  & canal_reclasificado!="EESS propias"   & p2>0)
}

*Creating matrices of mean utilities for diesel (j=1) and gasoline (j=2) for every region, only for those observations that belong to groups previously defined:
forvalues j=1/2{
forvalues i=1/5{
mkmat utilidad_media if grupo_`j'`i'==1, mat(um_`j'`i')
}
}
*Creating matrices with regresors for diesel (j=1) and gasoline (j=2) for every region, only for those observations that belong to groups previously defined::
forvalues i=1/5{
mkmat p1  premium tarjeta_puntos competencia_mayorista_provincia eess_prop_tot_prov urbanizacion ln_parque_automotor constante if grupo_1`i'==1, mat(x_1`i')
mkmat p2  premium tarjeta_puntos competencia_mayorista_provincia eess_prop_tot_prov urbanizacion  ln_parque_automotor constante if grupo_2`i'==1, mat(x_2`i')
}


*General coefficients:
*IV estimations for every group are generated to latter use:

*Diesel:
forvalues i=1/5{
quietly ivregress gmm utilidad_media (promedio_bandera=  p1) premium tarjeta_puntos competencia_mayorista_provincia eess_prop_tot_prov urbanizacion ln_parque_automotor ///
if (identificador==1 & utilidad_media>0 & gasoleos==1  & anio!=2015 & canal_reclasificado!="EESS propias" & REGIONES==`i' & p1>0)
matrix b_1_`i'= e(b) 
matrix var_b_1_`i'= e(V)
scalar n_1_`i'= e(N)
scalar rss_1_`i'= e(rss)
}


*Gasoline:
forvalues i=1/5{
quietly ivregress gmm utilidad_media (promedio_bandera=  p2) premium tarjeta_puntos competencia_mayorista_provincia eess_prop_tot_prov urbanizacion ln_parque_automotor ///
if (identificador==1 & utilidad_media>0 & naftas==1  & anio!=2015 & canal_reclasificado!="EESS propias" & REGIONES==`i' & p2>0)
matrix b_2_`i'= e(b) 
matrix var_b_2_`i'= e(V)
scalar n_2_`i'= e(N)
scalar rss_2_`i'= e(rss)
}


*Coefficients by product (diesel/gasoline)

forvalues j=1/2{
mat define beta_barra_`j'= 1/5*(b_`j'_1 + b_`j'_2 +b_`j'_3 +b_`j'_4+ b_`j'_5)
*mat list beta_barra_`j'
}
 
 
*Standard deviation for every product and panel:

*Diesel:
forvalues j=1/1{
forvalues i=1/5{
scalar sigma_`j'_`i'`i'= rss_`j'_`i'/ (n_`j'_`i' - 8- 8) 
}
}

*Gasoline:
forvalues j=2/2{
forvalues i=1/5{
scalar sigma_`j'_`i'`i'= rss_`j'_`i'/ (n_`j'_`i' - 8- 9) 
}
}

*Auxiliary calculations

forvalues j=1/2{
forvalues i=1/5{
mat define aux_`j'_`i'= (b_`j'_`i')' * b_`j'_`i'
}
mat define aux_suma_`j'= aux_`j'_1 + aux_`j'_2 + aux_`j'_3 + aux_`j'_4 + aux_`j'_5
mat define SIGMA_`j'= [aux_suma_`j'- 5* (beta_barra_`j')'*beta_barra_`j' ]/5
}

forvalues j=1/2{
mat define inv_SIGMA_`j'=inv(SIGMA_`j')
}


*Following Swamy paper, estimation of A_i for every j (j=1 diesel, j=2 gasoline), A_ji:

forvalues j=1/2{
forvalues i=1/5{
mat define aux_xx_`j'_`i'= (x_`j'`i')' * (x_`j'`i')
mat define A`j'_`i'= inv(inv_SIGMA_`j' + (1/sigma_`j'_`i'`i') * aux_xx_`j'_`i') * inv_SIGMA_`j'
}
}

*Auxiliary calculations and variance calculation for beta's (for every product: diesel j=1 and gasoline j=2):

forvalues j=1/2{
forvalues i=1/5{
mat define aux_var_beta_`j'_`i'=(sigma_`j'_`i'`i') * inv(aux_xx_`j'_`i')
}
}

forvalues j=1/2{
mat define var_beta_`j'= inv(inv( SIGMA_`j' + aux_var_beta_`j'_1 ) + inv(SIGMA_`j' + aux_var_beta_`j'_2 ) +inv(SIGMA_`j' + aux_var_beta_`j'_3 ) + inv(SIGMA_`j' + aux_var_beta_`j'_4 ) + inv(SIGMA_`j' + aux_var_beta_`j'_5 ))
}


*ASY- VAR(beta_i)

mat define  A_id=I(8)

forvalues j=1/2{
forvalues i=1/5{
mat define var_beta_`j'_`i'= var_beta_`j'+ (A_id - A`j'_`i') * (var_b_`j'_`i' - var_beta_`j') * (A_id - A`j'_`i')'
}
}


*Variance:
forvalues j=1/2{
forvalues i=1/5{
mat define var_beta_`j'_`i'_final= vecdiag( var_beta_`j'_`i')
}
}

forvalues j=1/2{
mat define var_beta_`j'_final= vecdiag(var_beta_`j')
}


*T statistic for each  beta_j_i (j=1 product diesel, j=2 product gasoline, with "i" indexing regions).
*Using mata to auxiliary calculations:

mata

var_beta_1_final = st_matrix("var_beta_1_final") 
var_beta_2_final = st_matrix("var_beta_2_final") 
var_beta_11_final = st_matrix("var_beta_1_1_final") 
var_beta_12_final = st_matrix("var_beta_1_2_final") 
var_beta_13_final = st_matrix("var_beta_1_3_final") 
var_beta_14_final = st_matrix("var_beta_1_4_final") 
var_beta_15_final = st_matrix("var_beta_1_5_final") 
var_beta_21_final = st_matrix("var_beta_2_1_final") 
var_beta_22_final = st_matrix("var_beta_2_2_final") 
var_beta_23_final = st_matrix("var_beta_2_3_final") 
var_beta_24_final = st_matrix("var_beta_2_4_final") 
var_beta_25_final = st_matrix("var_beta_2_5_final") 


*se
se_1_final = sqrt(var_beta_1_final) 
se_2_final = sqrt(var_beta_2_final) 
se_11_final = sqrt(var_beta_11_final) 
se_12_final = sqrt(var_beta_12_final) 
se_13_final = sqrt(var_beta_13_final) 
se_14_final = sqrt(var_beta_14_final) 
se_15_final = sqrt(var_beta_15_final) 
se_21_final = sqrt(var_beta_21_final) 
se_22_final = sqrt(var_beta_22_final) 
se_23_final = sqrt(var_beta_23_final) 
se_24_final = sqrt(var_beta_24_final) 
se_25_final = sqrt(var_beta_25_final) 

st_matrix("se_1_final", se_1_final )
st_matrix("se_2_final", se_2_final )
st_matrix("se_11_final",se_11_final)
st_matrix("se_12_final",se_12_final)
st_matrix("se_13_final",se_13_final)
st_matrix("se_14_final",se_14_final)
st_matrix("se_15_final",se_15_final) 
st_matrix("se_21_final" ,se_21_final) 
st_matrix("se_22_final" ,se_22_final) 
st_matrix("se_23_final" ,se_23_final) 
st_matrix("se_24_final" ,se_24_final) 
st_matrix("se_25_final" ,se_25_final) 

end

*Auxiliary calculations: 

forvalues j=1/2{
forvalues i=1/5{
mat define  aux_W_`j'`i'= inv( SIGMA_`j' + aux_var_beta_`j'_`i' ) 
}
}


forvalues j=1/2{
forvalues i=1/5{
mat define  W_`j'`i'=  var_beta_`j'* aux_W_`j'`i'
}
}


*Controlling matrix (to be sure that the sum of assigned weights equals one:
forvalues j=1/2{
mat define control_W_`j'= W_`j'1+W_`j'2+W_`j'3+W_`j'4+W_`j'5
mat list control_W_`j'
}

*General beta coefficients for every product :  beta_sombrero_swammy_`j'
forvalues j=1/2{
forvalues i=1/5{
matrix define b_`j'_`i'= (b_`j'_`i')'
}
}

forvalues j=1/2{
mat define beta_sombrero_swammy_`j'=  W_`j'1 * b_`j'_1 +  W_`j'2* b_`j'_2+   W_`j'3* b_`j'_3  +   W_`j'4* b_`j'_4 +   W_`j'5* b_`j'_5
}

*Beta coefficients by product (j)and regions (i):beta_swammy_`j'`i'
forvalues j=1/2{
forvalues i=1/5{
mat define aux_beta_`j'`i'= (aux_xx_`j'_`i' * b_`j'_`i')/sigma_`j'_`i'`i' + inv_SIGMA_`j'*beta_sombrero_swammy_`j'
}
}

forvalues j=1/2{
forvalues i=1/5{
mat define beta_swammy_`j'`i'= inv(inv_SIGMA_`j' + (1/sigma_`j'_`i'`i') * aux_xx_`j'_`i') * aux_beta_`j'`i'
}
}


*Auxiliary calculations:
forvalues j=1/2{
mat define beta_sombrero_swammy_`j'=(beta_sombrero_swammy_`j')'
forvalues i=1/5{
mat define beta_swammy_`j'`i'=(beta_swammy_`j'`i')'
}
}

mata

beta_sombrero_swammy_1 = st_matrix("beta_sombrero_swammy_1") 
beta_sombrero_swammy_2 = st_matrix("beta_sombrero_swammy_2") 

beta_swammy_11 = st_matrix("beta_swammy_11") 
beta_swammy_12 = st_matrix("beta_swammy_12") 
beta_swammy_13 = st_matrix("beta_swammy_13") 
beta_swammy_14 = st_matrix("beta_swammy_14") 
beta_swammy_15 = st_matrix("beta_swammy_15") 

beta_swammy_21 = st_matrix("beta_swammy_21") 
beta_swammy_22 = st_matrix("beta_swammy_22") 
beta_swammy_23 = st_matrix("beta_swammy_23") 
beta_swammy_24 = st_matrix("beta_swammy_24") 
beta_swammy_25 = st_matrix("beta_swammy_25") 

t_valor_1=beta_sombrero_swammy_1 :/ se_1_final
t_valor_2=beta_sombrero_swammy_2 :/ se_2_final
t_valor_11=beta_swammy_11 :/ se_11_final
t_valor_12=beta_swammy_12 :/ se_12_final
t_valor_13=beta_swammy_13 :/ se_13_final
t_valor_14=beta_swammy_14 :/ se_14_final
t_valor_15=beta_swammy_15 :/ se_15_final
t_valor_21=beta_swammy_21 :/ se_21_final
t_valor_22=beta_swammy_22 :/ se_22_final
t_valor_23=beta_swammy_23 :/ se_23_final
t_valor_24=beta_swammy_24 :/ se_24_final
t_valor_25=beta_swammy_25 :/ se_25_final

 
*To stata:

t_valor_1=st_matrix("t_valor_1",t_valor_1)
t_valor_2=st_matrix("t_valor_2",t_valor_2)
t_valor_11=st_matrix("t_valor_11",t_valor_11)
t_valor_12=st_matrix("t_valor_12",t_valor_12)
t_valor_13=st_matrix("t_valor_13",t_valor_13)
t_valor_14=st_matrix("t_valor_14",t_valor_14)
t_valor_15=st_matrix("t_valor_15",t_valor_15)
t_valor_21=st_matrix("t_valor_21",t_valor_21)
t_valor_22=st_matrix("t_valor_22",t_valor_22)
t_valor_23=st_matrix("t_valor_23",t_valor_23)
t_valor_24=st_matrix("t_valor_24",t_valor_24)
t_valor_25=st_matrix("t_valor_25",t_valor_25)
 
 end
 
forvalues j=1/2{
forvalues i=1/5{
mat list t_valor_`j'`i'
}
}

*Degrees of fredom by product (j) and region:
 
*Diesel
forvalues j=1/1{
forvalues i=1/5{
matrix define df_`j'`i'= (n_`j'_`i' - 8 - 8, n_`j'_`i' - 8 - 8, n_`j'_`i' - 8 - 8, n_`j'_`i' - 8 - 8, ///
n_`j'_`i' - 8 - 8, n_`j'_`i' - 8 - 8, n_`j'_`i' - 8 - 8, n_`j'_`i' - 8 - 8)
}
}
*Gasoline:
forvalues j=2/2{
forvalues i=1/5{
matrix define df_2`i'= (n_`j'_`i' - 8 - 9, n_`j'_`i' - 8 - 9, n_`j'_`i' - 8 - 9, n_`j'_`i' - 8 - 9, ///
n_`j'_`i' - 8 - 9, n_`j'_`i' - 8 - 9, n_`j'_`i' - 8 - 9, n_`j'_`i' - 8 - 9)
}
}


sum  p1 if (identificador==1 & utilidad_media>0 & gasoleos==1 & anio!=2015  & canal_reclasificado!="EESS propias"   & p1>0)
sum  p2 if (identificador==1 & utilidad_media>0 & naftas==1 & anio!=2015  & canal_reclasificado!="EESS propias"  & p2>0)

/*Degrees of freedom for general betas beta_j: observations minus degrees of fredom lost by estimating individual coefficients of every
panel (8*5) and general coefficients for the whole database (8) minus degrees od fredom lost by instrumenting price in every panel (8*5) for 
diesel productos and (9*5) for gasoline products:
*/

 mat define df_1= (16130, 16130, 16130, 16130,16130, 16130,16130,16130)
 mat define df_2= (10240, 10240,10240, 10240, 10240, 10240, 10240, 10240)   

*Auxiliary calculations:

mata
df_1=st_matrix("df_1")
df_2=st_matrix("df_2")
df_11=st_matrix("df_11")
df_12=st_matrix("df_12")
df_13=st_matrix("df_13")
df_14=st_matrix("df_14")
df_15=st_matrix("df_15")
df_21=st_matrix("df_21")
df_22=st_matrix("df_22")
df_23=st_matrix("df_23")
df_24=st_matrix("df_24")
df_25=st_matrix("df_25")


*p-valor: ttail(df, t_valor)

t_valor_1=beta_sombrero_swammy_1 :/ se_1_final
t_valor_2=beta_sombrero_swammy_2 :/ se_2_final
t_valor_11=beta_swammy_11 :/ se_11_final
t_valor_12=beta_swammy_12 :/ se_12_final
t_valor_13=beta_swammy_13 :/ se_13_final
t_valor_14=beta_swammy_14 :/ se_14_final
t_valor_15=beta_swammy_15 :/ se_15_final
t_valor_21=beta_swammy_21 :/ se_21_final
t_valor_22=beta_swammy_22 :/ se_22_final
t_valor_23=beta_swammy_23 :/ se_23_final
t_valor_24=beta_swammy_24 :/ se_24_final
t_valor_25=beta_swammy_25 :/ se_25_final


p_valor_1= ttail(df_1, abs(t_valor_1))*2
p_valor_2= ttail(df_2, abs(t_valor_2))*2
p_valor_11= ttail(df_11, abs(t_valor_11))*2
p_valor_12= ttail(df_12, abs(t_valor_12))*2
p_valor_13= ttail(df_13, abs(t_valor_13))*2
p_valor_14= ttail(df_14, abs(t_valor_14))*2
p_valor_15= ttail(df_15, abs(t_valor_15))*2

p_valor_21= ttail(df_21, abs(t_valor_21))*2
p_valor_22= ttail(df_22, abs(t_valor_22))*2
p_valor_23= ttail(df_23, abs(t_valor_23))*2
p_valor_24= ttail(df_24, abs(t_valor_24))*2
p_valor_25= ttail(df_25, abs(t_valor_25))*2
 
 
*Taking  p-valor to stata

p_valor_1=st_matrix("p_valor_1",p_valor_1)
p_valor_2=st_matrix("p_valor_2",p_valor_2)

p_valor_11=st_matrix("p_valor_11",p_valor_11)
p_valor_12=st_matrix("p_valor_12",p_valor_12)
p_valor_13=st_matrix("p_valor_13",p_valor_13)
p_valor_14=st_matrix("p_valor_14",p_valor_14)
p_valor_15=st_matrix("p_valor_15",p_valor_15)

p_valor_21=st_matrix("p_valor_21",p_valor_21)
p_valor_22=st_matrix("p_valor_22",p_valor_22)
p_valor_23=st_matrix("p_valor_23",p_valor_23)
p_valor_24=st_matrix("p_valor_24",p_valor_24)
p_valor_25=st_matrix("p_valor_25",p_valor_25)


end 

*Testing null hypothesis that beta from different between themselfs:
*H0: beta_ji=bjh para h!= i.

*Auxiliary calculations:
forvalues j=1/2{
mat define beta_cruz_`j'= inv((1/sigma_`j'_11)*aux_xx_`j'_1+(1/sigma_`j'_22)*aux_xx_`j'_2+(1/sigma_`j'_33)*aux_xx_`j'_3+(1/sigma_`j'_44)*aux_xx_`j'_4+(1/sigma_`j'_55)*aux_xx_`j'_5)*((1/sigma_`j'_11)*aux_xx_`j'_1* b_`j'_1+(1/sigma_`j'_22)*aux_xx_`j'_2* b_`j'_2+(1/sigma_`j'_33)*aux_xx_`j'_3* b_`j'_3+(1/sigma_`j'_44)*aux_xx_`j'_4* b_`j'_4+(1/sigma_`j'_55)*aux_xx_`j'_5* b_`j'_5)
}
forvalues j=1/2{
forvalues i=1/5{
mat define aux_t_`j'_`i'=(b_`j'_`i'-beta_cruz_`j')' * ((1/sigma_`j'_`i'`i')*aux_xx_`j'_`i') * (b_`j'_`i'-beta_cruz_`j')
}
}

forvalues j=1/2{
mat define T_`j'=(aux_t_`j'_1+aux_t_`j'_2+aux_t_`j'_3+aux_t_`j'_4+aux_t_`j'_5)
mat list T_`j'
}

*Squared Xi, with k(P-1) df. k: estimated parameters (8),  P panels (5), then 40 df.
*Values indicate that estimation for both products (j=1 diesel and j=2 gasoline) in every panel are statistically different.

*Exporting results:
 
*Auxiliary steps to simplify syntax:
forvalues j=1/2{
drop identificador_`j'
}
bysort anio mes canal_reclasificado Provincia  empresa_reclasificada: gen identificador_1= 1 if (identificador==1 & utilidad_media>0 & gasoleos==1 & anio!=2015  & canal_reclasificado!="EESS propias"  & p1>0)																													
bysort anio mes canal_reclasificado Provincia  empresa_reclasificada: gen identificador_2= 1 if (identificador==1 & utilidad_media>0 & naftas==1 & anio!=2015  & canal_reclasificado!="EESS propias"  & p2>0)																													


 forvalues j=1/2{
 mat define beta_sombrero_swammy_`j'= (beta_sombrero_swammy_`j')'
 mat define p_valor_`j'=(p_valor_`j')'
 }

   
 
 *Coefficients by region and product:
 
 *Auxiliary steps:
 forvalues j=1/2{
 forvalues i=1/5{
 mat define beta_swammy_`j'`i'= (beta_swammy_`j'`i')'
 mat define p_valor_`j'`i'= (p_valor_`j'`i')'
 }
 }
 
 *Matrices with coefficients by region and its  p-value:
  
 
*Diesel:
 mat define matriz_1=(beta_swammy_11, p_valor_11,beta_swammy_12, p_valor_12 , beta_swammy_13, p_valor_13 ,beta_swammy_14, p_valor_14 , beta_swammy_15, p_valor_15)
 matrix rownames matriz_1=	 "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const."				
*Gasoline:
 mat define matriz_2=(beta_swammy_21, p_valor_21,beta_swammy_22, p_valor_22 , beta_swammy_23, p_valor_23 ,beta_swammy_24, p_valor_24 , beta_swammy_25, p_valor_25)
 matrix rownames matriz_2=	"Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const."				

 forvalues j=1/2{
mat define matriz_`j'=(matriz_`j')'
*mat list matriz_`j'
}

 forvalues j=1/2{
 mat rownames matriz_`j' = "Cuyo" "p-valor" "Patagonia" "p-valor" "Northeast" "p-valor" "Northwest" "p-valor" "Pampeana" "p-valor"
}
 forvalues j=1/2{
 esttab matrix(matriz_`j', fmt(2)) , varwidth(7) modelwidth(7) abbrev substitute(p_valor  ) title("Random coefficient model results by region, and product (1= diesel, 2= gasoline).")  replace
 }
 
 *Auxiliary data manipulation:
 
 forvalues j=1/2{
 mat define std_error_`j'=se_`j'_final'
 forvalues i=1/5{
 mat define std_error_`j'`i'=se_`j'`i'_final'
 }
 }
 
 
 *Diesel
 mat define betas_diesel= (beta_sombrero_swammy_1\beta_swammy_11\beta_swammy_12\ beta_swammy_13\beta_swammy_14\ beta_swammy_15)
 mat define std_deviation_diesel= (std_error_1\std_error_11\std_error_12\std_error_13\std_error_14\std_error_15)
 mat define p_values_diesel= ( p_valor_1\ p_valor_11\p_valor_12\p_valor_13\p_valor_14\p_valor_15)
 matrix rownames betas_diesel=	 "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const."				
 matrix rownames std_deviation_diesel= "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const."				
 matrix rownames p_values_diesel=	 "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const."				
  
 mat define diesel_swamy=(betas_diesel, std_deviation_diesel, p_values_diesel)
 matrix rownames  diesel_swamy=	 "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Cuyo: Price" "Cuyo: Premium" "Cuyo: Reward card" "Cuyo: Wholesale competitors" "Cuyo: Flagged outlets" "Cuyo: Urbanization"  "Cuyo: Vehicle fleet" "Cuyo: Const." "Patagonia: Price" "Patagonia: Premium" "Patagonia: Reward card" "Patagonia: Wholesale competitors" "Patagonia: Flagged outlets" "Patagonia: Urbanization" "Patagonia: Vehicle fleet" "Patagonia: Const." "Northeast: Price" "Northeast: Premium" "Northeast: Reward card" "Northeast: Wholesale competitors" "Northeast: Flagged outlets" "Northeast: Urbanization" "Northeast: Vehicle fleet" "Northeast: Const." "Northwest: Price" "Northwest: Premium" "Northwest: Reward card" "Northwest: Wholesale competitors" "Northwest: Flagged outlets" "Northwest: Urbanization"  "Northwest: Vehicle fleet" "Northwest: Const." "Pampeana: Price" "Pampeana: Premium" "Pampeana: Reward card" "Pampeana: Wholesale competitors" "Pampeana: Flagged outlets" "Pampeana: Urbanization"  "Pampeana: Vehicle fleet" "Pampeana: Const."				
 matrix colnames diesel_swamy= Coef Std_Error p-value
 
 *Diesel markup by region and firm:
 esttab matrix(diesel_swamy, fmt(2)), title(Swamy random coefficient's model results with instrumented price.  General coefficients and coefficients by region. Product: Diesel.)  t(3) b(a3) p(2) 
 
*Gasoline:

 mat define betas_gasoline= (beta_sombrero_swammy_2\beta_swammy_21\beta_swammy_22\ beta_swammy_23\beta_swammy_24\ beta_swammy_25)
 mat define std_deviation_gasoline= (std_error_2\std_error_21\std_error_22\std_error_23\std_error_24\std_error_25)
 mat define p_values_gasoline= ( p_valor_2\ p_valor_21\p_valor_22\p_valor_23\p_valor_24\p_valor_25)
 matrix rownames betas_gasoline=	 "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const."				
 matrix rownames std_deviation_gasoline= "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const."				
 matrix rownames p_values_gasoline=	 "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const."				
  
 mat define  gasoline_swamy=(betas_gasoline, std_deviation_gasoline, p_values_gasoline)
 matrix rownames  gasoline_swamy=	 "Price" "Premium" "Reward card" "Wholesale competitors" "Flagged outlets" "Urbanization"  "Vehicle fleet" "Const." "Cuyo: Price" "Cuyo: Premium" "Cuyo: Reward card" "Cuyo: Wholesale competitors" "Cuyo: Flagged outlets" "Cuyo: Urbanization"  "Cuyo: Vehicle fleet" "Cuyo: Const." "Patagonia: Price" "Patagonia: Premium" "Patagonia: Reward card" "Patagonia: Wholesale competitors" "Patagonia: Flagged outlets" "Patagonia: Urbanization" "Patagonia: Vehicle fleet" "Patagonia: Const." "Northeast: Price" "Northeast: Premium" "Northeast: Reward card" "Northeast: Wholesale competitors" "Northeast: Flagged outlets" "Northeast: Urbanization" "Northeast: Vehicle fleet" "Northeast: Const." "Northwest: Price" "Northwest: Premium" "Northwest: Reward card" "Northwest: Wholesale competitors" "Northwest: Flagged outlets" "Northwest: Urbanization"  "Northwest: Vehicle fleet" "Northwest: Const." "Pampeana: Price" "Pampeana: Premium" "Pampeana: Reward card" "Pampeana: Wholesale competitors" "Pampeana: Flagged outlets" "Pampeana: Urbanization"  "Pampeana: Vehicle fleet" "Pampeana: Const."				
 matrix colnames gasoline_swamy= Coef Std_Error p-value
 
 *Gasoline markup by region and firm:
 esttab matrix(gasoline_swamy, fmt(2)), title(Swamy random coefficient's model results with instrumented price.  General coefficients and coefficients by region. Product: Gasoline.)  t(3) b(a3) p(2) 

 
 **********************************************************MARKUPS*************************************************************************************************************************************************																															

 /*Once coefficients by product have been estimated , markups are calculated for regular and premium products, because prices 
 as well as market shares, differ significantly between varieties.
 Products are indexed here (Markup section), as follows:
 1= regular diesel
 2= premium diesel
 3= regular gasoline
 4= premium gasoline
 
 */

 *Random coefficient model's markups estimation.
 
 *Auxiliary calculations: 
 																														
forvalues j=1/4{																															
drop px_promedio_bandera_p`j'
bysort empresa_reclasificada: egen px_promedio_bandera_p`j'= wtmean(precios_netos_actuales) if (producto_reclasificado==`j')	, weight(volumen_empresa_periodo)																														
}
*Average prices over the country by product and firm:
 estimates clear
 forvalues j=1/4{
 forvalues i=1/13{
 summarize px_promedio_bandera_p`j' if  empresa==`i' & producto_reclasificado==`j', meanonly  
 scalar define px_promedio_empresa`j'`i' = r(mean)
 }
}
*Average market share over the country by product and firm:
 forvalues j=1/4{
 forvalues i=1/13{
 summarize partic_mercado_empresa if  empresa==`i' & producto_reclasificado==`j' & identificador==1, meanonly  
 scalar define partic_mercado_empresa`i'_prod_`j' = r(mean)
 }
 }
 
 *Prices' coefficients:
 
 forvalues j=1/2{
 scalar coef_precio_`j'= beta_sombrero_swammy_1[1,1]
 scalar define inv_coef_precio_`j'=-(1/coef_precio_`j')
 }
 
  forvalues j=3/4{
 scalar coef_precio_`j'= beta_sombrero_swammy_2[1,1]
 scalar define inv_coef_precio_`j'=-(1/coef_precio_`j')
 }
 
 
 forvalues j=1/4{
 forvalues i=1/13{
 scalar define markup_prod_`j'_`i'= (partic_mercado_empresa`i'_prod_`j'/px_promedio_empresa`j'`i')*inv_coef_precio_`j'
 }
 }
 
 forvalues j=1/4{
 mat define markup_prod_`j'= (markup_prod_`j'_1 \markup_prod_`j'_2 \markup_prod_`j'_3 \markup_prod_`j'_4 \markup_prod_`j'_5 \markup_prod_`j'_6 \markup_prod_`j'_7 \markup_prod_`j'_8 \markup_prod_`j'_9 \ markup_prod_`j'_10 \ markup_prod_`j'_11 \markup_prod_`j'_12 \markup_prod_`j'_13)
 }
 
 *Average markups by product and firm for the whole country:
 
 mat define markups_predichos= (markup_prod_1, markup_prod_2, markup_prod_3, markup_prod_4)
 matrix colnames markups_predichos=   "Diesel" "Premium Diesel" "Gasoline" "Premium Gasoline"
 matrix rownames markups_predichos=	 Axion Dapsa Gulf Oil_Comb Others Pampa_En PDVSur Petrobras Puma Refinor Shell SCP YPF																							

 /*IMPORTANT: "Others" markup was not taken into accoount in the manuscript, because it does not represent a unique firm's sales, 	
but the sum of small companies sales, without significant individual market share .*/

 *General results:
 esttab matrix(markups_predichos, fmt(2)) , varwidth(7) modelwidth(7) abbrev title("Estimated markups by firm and product, general results for Argentina .")  replace



																													
*Prices by region and firm:
 
forvalues j=1/4{
forvalues i=1/5{																															
drop px_promedio_bandera_p`j'_r_`i'
bysort empresa_reclasificada: egen px_promedio_bandera_p`j'_r_`i'= wtmean(precios_netos_actuales) if (producto_reclasificado==`j' & REGIONES==`i')	, weight(volumen_empresa_periodo)																														
}
}
 
*Auxiliary calculations:

 forvalues j=1/4{
 forvalues e=1/13{
 forvalues i=1/5{
 summarize px_promedio_bandera_p`j'_r_`i' if  empresa==`e', meanonly  
 scalar define px_promedio_bandera_p`j'_r`i'_e`e' = r(mean)
 summarize partic_mercado_empresa if  empresa==`e' & producto_reclasificado==`j' & REGIONES==`i', meanonly  
 scalar define partic_mercado_prod_`j'_r`i'_e`e' = r(mean)
 scalar define aux_swamy_p`j'_r_`i'_e`e'=partic_mercado_prod_`j'_r`i'_e`e'/px_promedio_bandera_p`j'_r`i'_e`e'
 }
 }
 }
 

 *Coefficients by product and region:
 
 forvalues j=1/2{
 forvalues i=1/5{
 scalar coef_precio_`j'`i'= beta_swammy_1`i'[1,1]
 scalar define inv_coef_precio_`j'`i'=-(1/coef_precio_`j'`i')
 }
 }
 
  forvalues j=3/4{
 forvalues i=1/5{
 scalar coef_precio_`j'`i'= beta_swammy_2`i'[1,1]
 scalar define inv_coef_precio_`j'`i'=-(1/coef_precio_`j'`i')
 }
 }
 
 *Markups by region and firm:
 
 forvalues j=1/4{
 forvalues e=1/13{
 forvalues i=1/5{
 mat define markup_p`j'_e`e'_r`i'= aux_swamy_p`j'_r_`i'_e`e' * inv_coef_precio_`j'`i'
 mat list markup_p`j'_e`e'_r`i'
 }
 }
 }
 
 forvalues j=1/4{
 forvalues e=1/13{
 mat define markups_predichos_p`j'_e`e'=(markup_p`j'_e`e'_r1,markup_p`j'_e`e'_r2, markup_p`j'_e`e'_r3, markup_p`j'_e`e'_r4, markup_p`j'_e`e'_r5)
 }
 }
 
*Products from 1/4:
*Firm 5 is not taken into account because it represents "others" firms sales.
 forvalues j=1/4{
 mat define markups_predichos_p`j' =(markups_predichos_p`j'_e1\ markups_predichos_p`j'_e2\ markups_predichos_p`j'_e3\markups_predichos_p`j'_e4\markups_predichos_p`j'_e6\ markups_predichos_p`j'_e7\markups_predichos_p`j'_e8\markups_predichos_p`j'_e9\markups_predichos_p`j'_e10\ markups_predichos_p`j'_e11\markups_predichos_p`j'_e12\ markups_predichos_p`j'_e13)
 matrix colnames markups_predichos_p`j' = Cuyo Patagonia Northeast Northwest Pampeana
 matrix rownames markups_predichos_p`j'=  Axion Dapsa Gulf Oil_Comb Pampa_En PDVSur Petrobras Puma Refinor Shell SCP YPF																							
 }


 
 
