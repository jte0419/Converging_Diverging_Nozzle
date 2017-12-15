# Converging_Diverging_Nozzle

Contained in this respository is the code needed to compute everything you've ever wanted to know about converging-diverging (CD) nozzles.  Two files are needed in order to run the code, and they must be placed in the same directory/folder:

* [GUI_CD_Nozzle_v2.m](GUI_CD_Nozzle_v2.m)
* [GUI_CD_Nozzle_v2.fig](GUI_CD_Nozzle_v2.fig)

When you run the *.m* file, a GUI will open, with the default solution plotted, corresponding to a specific heat ratio of 1.4, an area ratio of 2, and a back pressure to stagnation pressure ratio of 0.55.  You can change any of the three values in the *Inputs* panel, and the solutions will automatically update.  I have updated the function used to compute the area ratio and Mach number, so the plotting is now much faster than it used to be, although it messes up the plotting near the throat in the converging portion of the nozzle.
