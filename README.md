# Hot Carrier Solar Cell Particle and Energy balance solver

Upper limit of hot carrier solar cells efficiency is solved by solving both particle balance and energy balance models with consideration of carrier cooling rate, energy selective contact energy level and optical phonon bottleneck effect. Realization of hot carrier solar cells depends on finding absorber materials with slow carrier cooling rate. As a part of research project, this code was made. Below left shows the upper limit of conversion efficiency with respect to band gap and carrier thermalization time for single barrier structure, where carriers are extracted via thermionic emission. If carrier cooling is very rapid, which is typically the case, the upper limit reaches the Shockley-Queisser limit. Below right shows the case with extraction of carrier with resonant tunneling. These results vary with extraction energy level. As you can imagine, if the extraction energy level is too high, not much carrier can be extracted although thouse carriers are highly energetic. If too low, turning back to Shockley-Queisser limit. There should be optimized energy levels for extraction. But it also depends on the thermalization time as can be seen in the figures.

<p align="center">
<img src="fig\SBEffEgw400p69.png" width=48%><img src="fig\RTDsw400p69.png" width=48%>
</p>

## Example

Import hcsc module

```python
>>> import hcsc
```

First run init() function in hcsc module. It will define 3 boolean global variables and set up fonts and errors.

```python
>>> hcsc.init()
```

Then create a hcsc object.

```python
>>> scell = hcsc.hcsc()
```
Material properties are defined in hcscAttribute class in hcscAttribute.py. Some properties are in bulk class in bulk.py and contacts class in contacts.py those are loaded into hcscAttribute class. Carrier extraction methods are also defined as functions.

Then create a Photon_in object.

```python
>>> ph_in = photon_in.Photon_in()
```
By default, it will assume AM15 spectrum. You might want to set a concentration factor, e.g. 10 (default 1)

```python
>>> ph_in.c = 10
```
Shine the light spectrum,

```python
>>> scell.shine(ph_in)
```
This will tell the solar cell how much photon particles and energy can be absorbed. Quantum efficiency is 1 for above band gap and 0 for below.

Almost ready. At this moment, you might want to display current properties of the cell.
```python
>>> scell.display_attributes()
Band gap (eV): 0.65
ESC(right) extraction energy (eV): 1.0
ESC(right) extraction energy width (eV): 0.001
Thermalization time (ps): 1.0
Electron phonon balance model (bool): False
Optical phonon decay time (ps): 8.0
Effective mass (electron mass): 0.1
```

To show hot carrier solar cell current (A/m<sup>2</sup>) at 0.2V bias,
```python
>>> print("Jout:{:.3f}(A/m^2)".format(scell.Jouthc(0.2/2*nu.eV)[0]))
Jout:40.030(A/m^2)
```
The parameter is the half the bias multiplied by q (elemntary charge). Jouthc returns 2 variables, current flux and boolean value telling solution is found or not.

To show hot carrier solar cell power (W/m<sup>2</sup>) at 0.2V bias,
```python
>>> print("Pout:{:.3f}(A/m^2)".format(scell.Pouthc(0.2/2*nu.eV)[0]))
Pout:8.006(W/m^2)
```

To show the max power under the given condition, use maxPouthc() method. First it will search through to find open circuit voltage Voc, then it will tell the max power.

```python
>>> print("MaxPout:{:.3f}(W/m^2)".format(scell.maxPouthc()[0]))
Solving Voc
Count:0 Tc:335.71 mue(eV):0.000 Bias(V):0.000 J(A/m^2):40.03 P(W/m^2):0.00 
Count:1 Tc:335.71 mue(eV):0.002 Bias(V):0.004 J(A/m^2):40.03 P(W/m^2):0.16 
.......
Count:249 Tc:335.97 mue(eV):0.498 Bias(V):0.996 J(A/m^2):7.74 P(W/m^2):7.70 
Count:250 Tc:335.99 mue(eV):0.500 Bias(V):1.000 J(A/m^2):5.14 P(W/m^2):5.14 
Count:251 Tc:336.01 mue(eV):0.502 Bias(V):1.004 J(A/m^2):2.34 P(W/m^2):2.35 
Count:252 Tc:336.03 mue(eV):0.504 Bias(V):1.008 J(A/m^2):-0.69 P(W/m^2):-0.70 
Stop iteration at mue(eV):5.040000e-01
MaxPout:32.423(W/m^2)
```

To show IV characteristics, you can easily make one now, but use test module, then use IV method

```python
>>> import test
>>> test.IV(scell)
```

for PV, 
```
>>> test.PV(scell)
```

Furthermore, it is necessary to find optimized extraction energy level at energy selective contact. To find this, use maxPouthc_opt_ESC method. It will scan through Power and show the list. The second paremter tells number of trials between 50% band gap energy and 300% band gap energy. It will also scan energy width (currently set at 1 meV and 5 meV)

```python
>>> test.maxPouthc_opt_ESC(scell, 50)
cnt1:0/17 cnt2:0/1
deltaEesc:0.65(eV) Ewesc:5.000e+00(meV)
Solving Voc
Count:0 Tc:300.12 mue(eV):0.000 Bias(V):0.000 J(A/m^2):4.34 P(W/m^2):0.00 
Count:1 Tc:300.12 mue(eV):0.002 Bias(V):0.004 J(A/m^2):4.34 P(W/m^2):0.02 
......
Count:137 Tc:300.12 mue(eV):0.274 Bias(V):0.548 J(A/m^2):4.31 P(W/m^2):2.36 
Count:138 Tc:300.09 mue(eV):0.276 Bias(V):0.552 J(A/m^2):-0.13 P(W/m^2):-0.07 
Stop iteration at mue(eV):2.760000e-01
less than before. skip optimization
Power (W/m^2) [[2.47181878e+03]
 [0.00000000e+00]
 [2.45677289e+03]
 [1.26399532e+03]
 [4.20439403e+01]
 [2.31927764e+00]
 [2.35939288e+00]
 [0.00000000e+00]
........
 [0.00000000e+00]
 [0.00000000e+00]]
Temperature (K) [[300.1147805 ]
 [  0.        ]
 [300.81875867]
 [300.0845472 ]
 [300.11986514]
 [300.11963268]
 [300.12493338]
 [  0.        ]
........
 [  0.        ]
 [  0.        ]]
maxP:2471.82 (W/m^2)
T:300.115(K)
dEesc:0.65(eV) Ewesc:5.000e+00(meV)
```

## Usage

Some class attributes and methods are omitted for readability.

```python
class hcsc(hcscAttribute):

    """ A class used to calculate hot carrier solar cell characteristics

    Methods
    -------
    shine(self, ph_in)
        do this first so that photons are absorbed in the cell
        set light source defined by photon_in class
        Calculate photon flux, energy flux, etc
    Jouthc(self, mue)
        Current output of the hot carrier
        return
    Vochc(self, fmaxP=False, dsp_msg=True)
        solve Voc
        mue should be in unit eV so that brentq gives accurate value
    Pouthc(self, mue)
        return poewr output at given mue
        note mue is measured from Eg/2
        symmetric two band model
        hence mue difference at contacts is 2*mue
    maxPouthc(self)
        solve max power at given Eg, ESC E and Ew
    """

class Photon_in(object):
    """ defining spectrum
    default spectrum is AM15 if n=2

    attributes
    ----------
    c : float
        concentratoin factor
    """

test.py:
def IV(scell):
    """Plot IV characteristics

    parameters
    ----------
    scell : class hcsc

    Returns
    -------
    none
    """	
	
def PV(scell):
    """Plot PV characteristics
    parameters
    ----------
    scell : class hcsc

    Returns
    -------
    none
    """

def maxPouthc_opt_ESC(scell, numEesc=50, g_skip=False):
    """solve maxpower at given Eg optimizing ESC E and Ew
    Energy level of ESC are between 50% of Eg to 300% of Eg

    Parameters
    ----------
    scell : class hcsc
    numEesc : int
        number of trials for ESC energy level
    g_skip : bool
        if g_skip is True, scan all possible values during optimization process possibly more accurate, but significantly increase process time

    Returns
    -------
    maxP : float
        W/m^2
    maxT : float
        K
    optEesc : float
        J
    optEwesc : float
        J
    """

class hcscAttribute:
    def __init__(self):
        """
        Attributes
        ----------
        Trm : float
            solar cell temperature
        tau_th : float
            carrier thermalization time
        tau_pp : float
            optical phonon decay time
        NM : float
            number of phonon modes
        Tph : float
            optical phonon temperature
        EP : bool
             if true, solve electron phonon balance equations
        d : float
            thickness of the cell if bulk structure is assumed
        Ephn : float
            dominant phonon energy (only used in polar optical (POP) scattering)
        hw0 : float
            dominant longitudinal optical (LO) phonon energy
        absb : class
            absorber class object defined in a separate file
        lcnt : class
            left electrical contact defined in a separte file
        rcnt : class
            right electrical contacts defined in a separete file
        lesc : class
            left energy selective contact defined in a separate file
        resc :class
            right energy selective contact defined in a separate file
        Jext : function
            current through energy selective contact
        Uext : function
            energy through energy selective contact
        Uth : function
            carrier energy thermalization rate
        _c1 : float
            constant variables
        _cJ : float
            constant variables
        _cS : float
            constant variables

        Methods
        -------
        display_attributes(self):
            print attributes
        JextTherm(muc, mue)
            current flux extraction via thermionic emission
        UextTherm(muc, mue)
            energy flux extraction via thermionic emission
        JextESC(self, muc, mue)
            extracted current through tunneling
            assume carrier selective contact
            electron current flows from absb to right contact (absb->resc->rcnt)
        UextESC(self, muc, mue)
            extracted energy flux
        JextRN(self, muc, mue)
            unused
        UextRN(self, muc, mue)
            unused
        Uthstandard(self, muc):
            standard thermalization
            consider dimension
            try 3D
            parameter is tau_th thermalization time

        """
```

## Background

Hot carrier solar cell is a third generation solar cell aiming to overcome Shockley-Queisser limit by utilizing excess energy above band gap which is normally lost in standard solar cells.

Ross and Nozik [1], as well as Wurfel [2] suggested the hot carrier solar cell (see figure below). The structure introduces the energy selective contact (ESC) in which carriers in the absorber are extracted at elevated energy E<sub>esc</sub>, e.g. electrons are extracted above the conduction band edge. The ESCs are located between the absorber and the regular contacts. The extraction energy has an energy band, w<sub>esc</sub>. In Ross and Nozik model, it is assumed to be a single energy level. ESC needs to be carrier selective. In the figure, electrons are extracted in the right hand side contact at energy E<sub>esc</sub> and holes are extracted at energy -E<sub>esc</sub> in the left hand side contact.

Chemical potential, or quasi-Fermi level in the absorber is &mu;<sub>c</sub> and that in contactc is denoted as &mu;<sub>e</sub>. If the carriers are extracted before thermalization to the band edges but after electron electron scattering which is super fast process, the carrier distributions can be described as heated Fermi-Dirac distribution with its own carrier temperature T<sub>c</sub>, which is higher than the room temperature T<sub>rm</sub>. For hot carrier solar cells to be effective, it is desirable to have narrow band gap materials although Wurfel suggests Auger recombinations are imminent. But in Ross and Nozik model, non-radiative recombinations are not considered. Nevertheless, it will provide the upper limit, just like in the Shockley-Queisser limit. The figure bottom right shows the upper limit of hot carrier solar cell in Ross and Nozik model with various carrier temperatures. If the hot carriers are significantly utilized, it can go beyond the Shockley-Queisser limit.

<img src="fig\SchematicHCSC.png" width=48%><img src="fig\Eff_Eg_T_1000_sun.png" width=48%>

The carrier thermalization process is very rapid , typically in the order of less than pico seconds. It is quite challenging to find an absorber material with slow cooling.

<p align="center">
<img src="fig\particle_balance.PNG" width=75%>
<img src="fig\energy_balance.PNG" width=75%>
</p>

**References:**

[1] Robert T. Ross and Arthur J. Nozik. Eficiency of hot-carrier solar energy converters. Journal of Applied
Physics, 53(5):3813{3818, 1982.

[2] Peter Wurfel. Solar energy conversion with hot electrons from impact ionisation. Solar Energy Materials
and Solar Cells, 46(1):43{52, 1997.


## Speeding up

3 pyx files and setup.py file are ready for Cython. If someone can speed up scipy.integrate.quad, it will be great. scipy.LowLevelCallable seems to be the way to go but I can't figure it out how to use it.

Other than that, solver can be improved by using differential method so that it can skip lots of computation along flat IV curve, for example.


## Prerequisites

Confirmed working with Python 3.7

Numpy, Scipy and MatPlotLib are required.

## Installing

Download files into a folder.
