# Hot Carrier Solar Cell DB (Detailed Balance) solver

Upper limit of hot carrier solar cells efficiency is solved by solving both particle balance and energy balance models with consideration of carrier cooling rate and energy selective contact energy level. Hence, non-radiative recombinations are not included. Realization of hot carrier solar cells depends on finding absorber materials with slow carrier cooling rate.

<p align="center">
<img src="fig\SBEffEgw400p69.png" width=48%>
</p>

Hot carrier solar cell is a third generation solar cell aiming to overcome Shockley-Queisser limit by utilizing excess energy above band gap which is normally lost in standard solar cells.

Ross and Nozik [1], as well as Wurfel [2] suggested the hot carrier solar cell (see figure below). The structure introduces the energy selective contact (ESC) in which carriers in the absorber are extracted at elevated energy E<sub>esc</sub>, e.g. electrons are extracted above the conduction band edge. The extraction energy has an energy band, w<sub>esc</sub>. But they assumed it is a single level of energy. ESC needs to be carrier selective. In the figure, electrons are extracted in the right hand side contact at energy E<sub>esc</sub> and holes are extracted at energy -E<sub>esc</sub> in the left hand side contact.

Chemical potential, or quasi-Fermi level in the absorber is &mu;<sub>c</sub> and that in contactc is denoted as &mu;<sub>e</sub>. If the carriers are extracted before thermalization to the band edges but after electron electron scattering which is super fast process, the carrier distributions can be described as heated Fermi-Dirac distribution with its own carrier temperature T<sub>c</sub>, which is higher than the room temperature T<sub>rm</sub>. For hot carrier solar cells to be effective, it is desirable to have narrow band gap materials although Wurfel suggests Auger recombinations are imminent. But in Ross and Nozik model, non-radiative recombinations are not considered. Nevertheless, it will provide the upper limit, just like in the Shockley-Queisser limit. The figure bottom right shows the upper limit of hot carrier solar cell in Ross and Nozik model with various carrier temperatures. If the hot carriers are significantly utilized, it can go beyond the Shockley-Queisser limit.

<img src="fig\SchematicHCSC.png" width=48%><img src="fig\Eff_Eg_T_1000_sun.png" width=48%>

The carrier thermalization process is very rapid , typically in the order of less than pico seconds. It is quite challenging to find an absorber material with slow cooling.


**References:**

[1] Robert T. Ross and Arthur J. Nozik. Eficiency of hot-carrier solar energy converters. Journal of Applied
Physics, 53(5):3813{3818, 1982.

[2] Peter Wurfel. Solar energy conversion with hot electrons from impact ionisation. Solar Energy Materials
and Solar Cells, 46(1):43{52, 1997.

## Example



## Usage


## Prerequisites

Confirmed working with Python 3.7

## Installing

Download files into a folder.
