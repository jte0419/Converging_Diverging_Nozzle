# Converging_Diverging_Nozzle

Contained in this respository is the code needed to compute everything you've ever wanted to know about converging-diverging (CD) nozzles.  Two files are needed in order to run the code, and they must be placed in the same directory/folder:

* [GUI_CD_Nozzle_v2.m](GUI_CD_Nozzle_v2.m)
* [GUI_CD_Nozzle_v2.fig](GUI_CD_Nozzle_v2.fig)

When you run the *.m* file, a GUI will open, with the default solution plotted, corresponding to a specific heat ratio of 1.4, an area ratio of 2, and a back pressure to stagnation pressure ratio of 0.55.  You can change any of the three values in the *Inputs* panel, and the solutions will automatically update.  I have updated the function used to compute the area ratio and Mach number, so the plotting is now much faster than it used to be, although it messes up the plotting near the throat in the converging portion of the nozzle.

## Solution Variables Explanation

The fields in the *Solutions* panel might need a little explanation.  We will start from the top and move our way down.  The first text box shows whether the flow in the nozzle is choked or not.  The only time the flow will not be choked is when we have fully subsonic isentropic flow throughout the nozzle, corresponding to a back-to-reservoir pressure ratio greater than *Pb/Po|sub*.

The next text box shows the state that the nozzle is operating in.  To get a primer on the different nozzle states, check out the [first video](https://www.youtube.com/watch?v=p8e8A3sdVOg) in the *YouTube Videos* section below.

The next two entries are nozzle exit Mach numbers.  The first one (*Msub*) is the exit Mach number when the Nozzle is operating in the choked subsonic isentropic flow state.  The second one (*Msup*) is the exit Mach number when the nozzle is operating in the choked supersonic isentropic flow state.  These will not change when the *Pb/Po* input is changed.  Note that when there is a normal shock in the nozzle or at the nozzle exit plane, the exit Mach number will be different than the Mach numbers listed in these two boxes.  I don't currently dispaly it in the GUI, but I'll update that in the future.

The next three boxes are the back-to-reservoir pressure ratios corresponding to the following three nozzle flow states: choked isentropic subsonic, choked isentropic supersonic, and normal shock at the exit plane.  Note that the pressure ratio for the normal shock at the exit is higher than the supersonic isentropic pressure ratio.  These are the only three pressure ratios values that you need to be able to delineate between flow states (check out the [third video](https://www.youtube.com/watch?v=b5q022xNgp0) in the *YouTube Videos* section below).

The last box is the area ratio at the location of the normal shock in the nozzle.  Note that a value will only appear in this box when the back-to-reservoir pressure ratio is below *Pb/Po|Sub* and above *Pb/Po|NS*, because in between these values is when a normal shock will form in the diverging section of the nozzle.  Also note that this is the area ratio of the normal shock location, not necessarily the distance along the nozzle of the normal shock.  That depends on how the area of the nozzle changes along its length.

## YouTube Videos

To learn more about the computations included in this code, check out my YouTube videos.  More videos will be coming soon.

* [Explained: Converging-Diverging Nozzle](https://www.youtube.com/watch?v=p8e8A3sdVOg)
* [Explained: Maximum Thrust Nozzle Exit Pressure Condition](https://www.youtube.com/watch?v=F3ylH8onOlw)
* [Converging-Diverging Nozzle Pressure Delineations](https://www.youtube.com/watch?v=b5q022xNgp0)
