twoprime-seq protocol
=====================

Preparing RNA
-------------
RNA prepared from log-phase cultures of yeast grown at 30C in YEPD.

Strains:

1. tpt1 (10X-tRNA TRP1) (TPT1 URA3)
2. tpt1 (10X-tRNA TRP1) 

RNA prepared by extraction with acid phenol into TES at 65C for 20
minutes. Precipitated with NaOAc.

Quantitation: XXX

Removal of 2-prime groups
-------------------------

Buffer (10X):
200 mM Tris-HCl (pH 7.5)
50 mM MgCl2
1 mM DTT
4% Triton X-100
100 mM NAD+

- Treat 10 micrograms of RNA in 1X Tpt1 buffer with 1 uM Tpt1 for 1 hour at
  30C.

- P:C:I extract and resuspend in H2O

Alkaline hydrolysis
-------------------

- Denature the RNA in a 200 μl PCR tube at 90°C for 2 min  a
  thermocycler.

- Add 0.1 M NaHCO3/Na2CO3 pH 9.9 in amount as the RNA solution and
  incubate at 90°C for 5 min

- Immediately add 10 μg glycogen and transfer the sample to a tube
  with 2.5 x volume 96% cold ethanol, 0.1 x volume 3 M NaAc and 0.025
  x volume 1 M Tris pH 7.5

- Precipitate by incubating on dry-ice for 15 minutes and spin at
  max speed for 30min.

Removal of 3-prime phosphate groups
-----------------------------------

Products of alkaline hydrolysis have 3-prime phosphate termini that must
be removed prior to the TGIRT-seq step (see PMID XXX).

Buffer (10X):
1 M Tris-acetate (pH 6.0)
100 mM MgCl2
20 mM DTT

  =============== ==============
  RNA (10 ug)     10 ug        
  H2O             X uL
  10 X buffer     2 uL
  T4 PNK (10 U)   2 uL
                  ----
                  10 uL

Incubate for 1 hour at 22 C.

Gel purification
----------------

- Purify fragments on 10% TBE-urea gel with small MW marker. Excise region
  with 20-40 nt fragments.

- Recover by agitating slice in 400 uL H2O, 400 uL P:C:I overnight at 4C.

- Chloroform extract aq phase and preciptate.

TGIRT-seq
---------

