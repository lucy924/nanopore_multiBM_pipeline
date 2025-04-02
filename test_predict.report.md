# BCG Susceptibility Panel for NMIBC
*This test combines the results for relevant single nucleotide variants, methylation variants, and immune infiltrate prediction to reach a final score for prediction of BCG treatment success for NMIBC.*  
**Panel name** BCG on NMIBC<br>
**URL:** https://github.com/lucy924/nanopore_multiBM_pipeline
___
### Clinical Details
**Sample Name:** test4predict<br>
Bladder lesion, transurethral resection:<br>
Papillary urothelial carcinoma, high grade (grade 3/3 - WHO 1973)<br>
Noninvasive, pTa<br>
**Specimen:** TURBT biopsy  
___
## Overall Result and Interpretation
<span style="font-size:16px;">Tumour is **Somewhat Likely** to respond to treatment.</span><br>

*Guidance: Result is out of Highly Unlikely, Somewhat Unlikely, Somewhat Likely, Highly Likely*

**Comments:** The result is an aggregate score obtained by combining the scores of biomarkers that reach read depth threshold and fall within respective scoring parameters.

___
## Immune infiltrate results

| Immune infiltrate proportions (full) | Immune infiltrate proportions (non-malignant only) |
| --- | --- |
| ![Immune infiltrate proportions (full)](results/test6/immune_infiltrate/immune_infiltrate_barplot.png) | ![Immune infiltrate proportions (non-malignant only)](results/test6/immune_infiltrate/immune_infiltrate_barplot_noca.png) |

<br>  
<br>  

___
## Raw results
*Guidance: Total (normalised) Score is a scale of 0% – 100%. 0% will respond and 100% will not respond to BCG treatment. Most results will fall between 25% and 75%. We recommend treatment with BCG when the result is \<50% for patients who are otherwise healthy, and \<40% for all other patients.*


| **ID** | Biomarker name | Biomarker type | Result Options | Result | Score | 
|---:|---|---|---|---|---|
|**004**|NOS3|snv|T\|T/T\|A/T\|G/A\|A/A\|G/G\|G|G\|G|0.16|
|**005**|NOS3|snv|C\|C/C\|T/T\|T|T\|C|0.53|
|**009**|CHST11|mod|0.0-1.0|0.1202386872202884|2.63118846345102|
|**010**|KLF8|mod|0.0-1.0|0.2421676545300592|1.6798814563928874|
|**011**|GPR158|mod|0.0-1.0|0.2618603411513859|2.9289312366737743|
|**012**|C12orf42|mod|0.0-1.0|0.2392036753445635|1.781592649310873|
|**013**|WDR44|mod|0.0-1.0|0.6964944649446494|1.1329889298892988|
|**014**|FLT1|mod|0.0-1.0|0.0413255360623781|1.5478070175438603|
|**022**|PMF1|mod|methylated/unmethylated|0.7149390243902439|1.21|
|**025**|IL1B|snv|G\|G/G\|A/A\|A|0.7149390243902439|1.0|
|**026**|IL1RAP|snv|G\|G/G\|A/A\|A|0.7149390243902439|1.0|
|**027**|IL1RAP|snv|G\|G/G\|C/C\|C|0.7149390243902439|1.0|
|**028**|IL1RAP|snv|T\|T/T\|A/T\|C/A\|A/A\|C/C\|C|0.7149390243902439|1.0|
|**029**|IL1RAP|snv|G\|G/G\|A/A\|A|0.7149390243902439|1.0|
|**030**|IL18R1|snv|C\|C/C\|T/T\|T|T\|T|0.16|
|**031**|IL18R1|snv|G\|G/G\|C/C\|C|T\|T|1.0|
|**032**|IL18R1|snv|C\|C/C\|A/A\|A|T\|T|1.0|
|**033**|IL18R1|snv|C\|C/C\|T/T\|T|T\|T|0.16|
|**035**|IL18R1|snv|G\|G/G\|A/A\|A|A\|A|0.16|
|**036**|IL18R1|snv|G\|G/G\|A/A\|A|A\|A|1.0|
|**037**|IL18R1|snv|C\|C/C\|T/T\|T|T\|T|0.16|
|**038**|IL18R1|snv|G\|G/G\|A/A\|A|T\|T|1.0|
|**401**|LMR|immune_ratio|0.0-10.0|2.973359663832316|4.946719327664632|
|**402**|NLR|immune_ratio|0.0-10.0|0.3423641081742103|1.6576358918257903|
|**403**|monocyte_inf|immune_inf|0.0-1.0|nan|nan|
|**405**|Bcell_inf|immune_inf|0.0-1.0|0.9044892569439369|1.7155914055551502|
|**406**|CD4_inf|immune_inf|0.0-1.0|0.851659066167062|0.336663626466825|
|**408**|CD8_inf|immune_inf|0.0-1.0|0.0668717129484986|0.0798769722309009|
|**409**|Treg_inf|immune_inf|0.0-1.0|0.7722534591745389|1.5218027673396315|
|**Total (raw)**|nan|nan|nan|nan|33.500679744344644|
|**Total (normalised to 0.0 - 100.0)**|nan|nan|nan|nan|47.05501206028529|

<!-- </div>   -->

___
## Methodology
Nanopore sequencing was used to detect single nucleotide variants, methylation status and immune status (by methylation) of the tumour sample. 
