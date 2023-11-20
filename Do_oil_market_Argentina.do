
*Version: stata 14.

*Call database. Configure according to database location: 

*use /.../oil_market_database_2019.dta 
*use /.../base_mayorista_final.dta

do Swamy_estimation_IV_2019.do

*Herfindahl-Hirschman Index

tab producto_reclasificado if (anio!=2015  & identificador_unico==1), sum(HHI)  																															


*Graphics:

/*You may need to install the follownig packages to display the images with their corresponding format:

ssc install grstyle, replace
ssc install palettes, replace
ssc install colrspace, replace
*/

grstyle init
grstyle set imesh 
grstyle set compact
grstyle set color economist 
graph drop _all

*Diesel 
graph hbox partic_empresa_3 if producto_reclasificado==1 & canal_reclasificado!="EESS propias" & anio!=2015 & identificador_unico==1, over(empresa_reclasificada, relabel(1 "AXION" 2 "DAPSA" 3 "GULF" 4 "OIL COMBUSTIBLES" 5 "OTHERS" 6 "PAMPA ENERGIA" 7 "PDV SUR" 8 "PETROBRAS" 9 "PUMA" 10 "REFINOR" 11 "SHELL" 12 "YPF" ) label(labsize(small))) allcategories nooutside ytitle(Percentage (%)) ytitle(, size(small) justification(left)) title(Regular Diesel, size(medsmall) justification(left)) legend(off) name(pr1, replace) box(1, color(dkgreen)) note(" ", span) 
graph hbox partic_empresa_3 if producto_reclasificado==2 & canal_reclasificado!="EESS propias" & anio!=2015 & identificador_unico==1, over(empresa_reclasificada, relabel(1 "AXION" 2 "DAPSA" 3 "GULF" 4 "OIL COMBUSTIBLES" 5 "OTHERS" 6 "PAMPA ENERGIA" 7 "PDV SUR" 8 "PETROBRAS" 9 "PUMA" 10 "REFINOR" 11 "SHELL" 12 "YPF" ) label(labsize(small))) allcategories nooutside ytitle(Percentage (%)) ytitle(, size(small) justification(left)) title(Premium Diesel, size(medsmall) justification(left)) legend(off) name(pr2, replace) box(1, color(teal)) note(" ", span) 
graph combine pr1 pr2, ycommon 
*Gasoline
graph hbox partic_empresa_3 if producto_reclasificado==3 & canal_reclasificado!="EESS propias" & anio!=2015 & identificador_unico==1, over(empresa_reclasificada, relabel(1 "AXION" 2 "DAPSA" 3 "GULF" 4 "OIL COMBUSTIBLES" 5 "OTHERS" 6 "PAMPA ENERGIA" 7 "PDV SUR" 8 "PETROBRAS" 9 "PUMA" 10 "REFINOR" 11 "SHELL" 12 "YPF" ) label(labsize(small))) allcategories nooutside ytitle(Percentage (%)) ytitle(, size(small) justification(left)) title(Regular Gasoline, size(medsmall) justification(left)) legend(off) name(pr3, replace) box(1, color(blue)) note(" ", span) 
graph hbox partic_empresa_3 if producto_reclasificado==4 & canal_reclasificado!="EESS propias" & anio!=2015 & identificador_unico==1, over(empresa_reclasificada, relabel(1 "AXION" 2 "DAPSA" 3 "GULF" 4 "OIL COMBUSTIBLES" 5 "OTHERS" 6 "PAMPA ENERGIA" 7 "PDV SUR" 8 "PETROBRAS" 9 "PUMA" 10 "REFINOR" 11 "SHELL" 12"YPF" ) label(labsize(small))) allcategories nooutside ytitle(Percentage (%)) ytitle(, size(small) justification(left)) title(Premium Gasoline, size(medsmall) justification(left)) legend(off) name(pr4, replace) box(1, color(navy)) note(" ", span) 
graph combine pr3 pr4, ycommon 


/* Estimation results
You may need to install the follownig packages to display Swamy regression results in a proper way:
  
  ssc install estout, replace */

*Ordinary Least Squares, Instrumental Variables (Fixed Effects) and Random Coefficient estimations for diesel and gasoiline :

*Diesel
reg utilidad_media promedio_bandera premium tarjeta_puntos competencia_mayorista_provincia  eess_prop_tot_prov urbanizacion ln_parque_automotor i.REGIONES if (anio!=2015 & identificador==1 & utilidad_media>0 & canal_reclasificado!="EESS propias" & gasoleos==1 & max_comp_REGION>0 & p1>0), vce(robust)
ivregress 2sls utilidad_media (promedio_bandera =  max_comp_REGION)  premium tarjeta_puntos competencia_mayorista_provincia  eess_prop_tot_prov urbanizacion   ln_parque_automotor i.REGIONES if (anio!=2015 & identificador==1 & utilidad_media>0 & canal_reclasificado!="EESS propias" & gasoleos==1 ), vce(robust)
esttab matrix(diesel_swamy, fmt(2)), title(Swamy random coefficient's model results with instrumented price.  General coefficients and coefficients by region. Product: Diesel.)  t(3) b(a3) p(2) 
																																																											
*Gasoline
reg utilidad_media promedio_bandera premium tarjeta_puntos competencia_mayorista_provincia  eess_prop_tot_prov urbanizacion  ln_parque_automotor i.REGIONES if (anio!=2015 & identificador==1 & utilidad_media>0 & canal_reclasificado!="EESS propias" & naftas==1 & promedio_otros_canales>0 & p2>0), vce(robust)
ivregress 2sls utilidad_media (promedio_bandera =  promedio_otros_canales refineria) premium tarjeta_puntos competencia_mayorista_provincia  eess_prop_tot_prov urbanizacion ln_parque_automotor i.REGIONES if (anio!=2015 & identificador==1 & utilidad_media>0 & canal_reclasificado!="EESS propias" & naftas==1 & p2>0 ), vce(robust)
esttab matrix(gasoline_swamy, fmt(2)), title(Swamy random coefficient's model results with instrumented price.  General coefficients and coefficients by region. Product: Gasoline.)  t(3) b(a3) p(2) 
																																																											
																																																										

*Estimated markups by firm:	
*General results:
esttab matrix(markups_predichos, fmt(2)) , varwidth(7) modelwidth(7) abbrev title("Estimated markups by firm and product, general results for Argentina .")  replace

*Results using coefficients by region:
 /*IMPORTANT: "Others" markup was not taken into account in the manuscript, because it does not represent a unique firm's sales, 	
but the sum of small companies sales, without significant individual market share 
Product 1: diesel
Product 2: diesel premium
Product 3: gasoline
Product 4: premium gasoline.*/
 
forvalues j=1/4{
esttab matrix(markups_predichos_p`j', fmt(2)) , varwidth(7) modelwidth(7)   title("Estimated markups by firm, region and product, for product `j'.")  replace
}

