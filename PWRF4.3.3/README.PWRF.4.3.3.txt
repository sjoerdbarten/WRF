Polar WRF 4.3.3 README
March 1, 2022 (edited March 10, 2022)

Well, here it is. This is the polar-optimized supplement for WRF version 4.3.3.
It follows the same form as earlier Polar WRF supplements developed by the Polar 
Meteorology Group of the Byrd Polar and Climate Research Center at The Ohio State University. 
Additional support is available at http://polarmet.osu.edu/hines/PWRF/ .
Please see the README files for earlier version of Polar WRF 
The standard WRF code can be downloaded from NCAR at 
https://www2.mmm.ucar.edu/wrf/users/download/get_source.html .

First, we request a favor that all PWRF users re-register at http://polarmet.osu.edu/PWRF/registration.php
so we can update our list. Thank you!

Next there have been ongoing issues with compiling the polar-modified version of the code
- possibly related to the files on the dyn_em subdirectory.
We recommend compiling the initial condition program real.exe with the standard/classic WRF
to produce the initial and boundary files. Then run wrf.exe with Polar WRF to produce 
the standard model output.

A new feature is the Vignon option for the P3 microphysics scheme. 

It is a compiler option in phys/module_mp_p3.F.PWRF4.3.3 (needs to be enabled)

The text must be editted to 

#define VIGNON_ICE_PHYSICS

for the option to be active. We only recommend this for the Antarctic and the Southern Ocean 
where it has been tested (Hines et al. 2021; Vignon et al. 2021). 
It sets small INP concentration numbers based upn McCluskey et al. (2018) observations.
It is possible something similar may work for pristine/supercooled cases in the Arctic 
with small ice-forming nuclei concentrations. However, it has not been tested for that,
and it is likely a difference ice nuclei formula should be used in the Arctic.

Again, it is not a default and needs to be enabled in module_mp_p3.F.PWRF4.3.3 . 
We found this improves the simulation of supercooled liquid clouds
over Antarctica. Again, it is not tested for the Arctic. The change is described in Hines et al. (2021, 
http://polarmet.osu.edu/PMG_publications/hines_bromwich_jgr_2021.pdf) and 
Vignon et al. (2021, https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2020JD033490).

The VIGNON adjustment is now added to the Morrison microphysics as a compiler option that can be enabled.


Additionally, Lesheng Bai of our Polar Meteorology Group has found through testing that
the NoahMP routines for the LSM physics work better than the older Noah LSM. 
We recommend PWRF users now use NoahMP for the LSM.

Here are list of modifications for PWRF 4.3.3:

QVQVSI (ICE SATURAION RATIO)  number control to reduce the warm temperature bias over the Antarctic ice sheet 
(change QVQVSI from 108% to 120%). 

Modification of the Noah MP glacier albedo over Antarctica: Albedo limitation > 0.8 over Antarctic 
(above 1000m) to reduce the strong warm bias and ~0.5 (below 1000m) to remove cold bias along coastal in summer. 

Revise the parameters of the Noah MP soil (snow) temperature scheme over ice sheet to reduce the strong diurnal cycle 
(cold bias at night time and warm bias during daytime) over the Antarctic ice sheet during summer. 


Added note: Please do not use the version of module_sf_noahmp_groundwater.F that originially came with the PWRF 4.3.3 
update. Apparently, it is designed for WRF version 4.4 not version 4.3 .

