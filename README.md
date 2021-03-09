# FLASH_RADIOTHERAPY
Geant4 simulation for flash radiotherapy. It is part of my ongoing Master's Thesis at the University of Pisa.
![GEANT$](/imagesforredame/simulation_eith_dosemap_zoom_nice.PNG)
## What is being simulated ?
![Novac 11](/imagesforredame/novac11.jpg)
The aim is to simulate the dose distribution of a Novac 11 electron beam in FLASH conditions in a water phantom.
After having obtained the PDD and profiles, a detector will be added to see what effects its presence will make.

## Specs
![Optics](/imagesforredame/ottica.PNG)
With SSD=80cm the applicator is made of a titanium window followed by aluminum monitor chambers covered with pmma, and then the applicator is pmma.

## Energy spectrum
![Spectrum](/imagesforredame/spettronovac.PNG)
The application uses a custom spectrum for the beam.
## Simulated Geometry
![Geo](/imagesforredame/geant4simulgoodnoev.PNG)
The geometry is as described in the optics sheet. In the images one can see in blue the area that will be reserved for placing the detector, reshaped accordingly.

## Dose distribution
![Dose](/imagesforredame/simulation_eith_dosemap_zoom_front.PNG)
The application outputs a 3d dosemap txt file, however one can visualize the distribution at end of run to see if it makes sense.

## Validation
![Validation](/Flash_ex_novo/VALIDATION/pdddata.png)
There is also a python script to confront experimental PDD with simulation.
## Profile 
![r100](/Flash_ex_novo/VALIDATION/r100.png)
![r80](/Flash_ex_novo/VALIDATION/r80.png)
![r50](/Flash_ex_novo/VALIDATION/r50.png)
Because of low statistics due to computation limitations, there are some fluctuations in the constant part of the profiles.
