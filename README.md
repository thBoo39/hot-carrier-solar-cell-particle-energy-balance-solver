# Hot Carrier Solar Cell Particle and Energy balance solver

Upper limit of hot carrier solar cells efficiency is solved by solving both particle balance and energy balance models with consideration of carrier cooling rate, energy selective contact energy level and optical phonon bottleneck effect. Realization of hot carrier solar cells depends on finding absorber materials with slow carrier cooling rate. As a part of research project, this code was made. Below left shows the upper limit of conversion efficiency with respect to band gap and carrier thermalization time for single barrier structure, where carriers are extracted via thermionic emission. If carrier cooling is very rapid, which is typically the case, the upper limit reaches the Shockley-Queisser limit. Below right shows the case with extraction of carrier with resonant tunneling. These results vary with extraction energy level. Optimization will be necessary. Slow carrier thermalization time is a very important factor realizing hot carrier solar cells. Those parameters can be explored in this project. Note that figures below are produced separately using this solver, compilation of results.

<p align="center">
<img src="fig\SBEffEgw400p69.png" width=48%><img src="fig\RTDsw400p69.png" width=48%>
</p>

## Example

Make sure you have ASTMG173.csv available by downloading from [NREL web site](https://rredc.nrel.gov/solar//spectra/am1.5/).

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
Band gap (eV): 0.354
ESC(right) extraction energy (eV): 0.37170000000000003
ESC(right) extraction energy width (eV): 0.001
Thermalization time (ps): 1.0
Electron phonon balance model (bool): False
Optical phonon decay time (ps): 8.0
Effective mass (electron mass): 0.023
```

To show hot carrier solar cell current (A/m<sup>2</sup>) at 0.2V bias,
```python
>>> print("Jout:{:.3f}(A/m^2)".format(scell.Jouthc(0.2/2*nu.eV)[0]))
Jout:5100.802(A/m^2)
```
The parameter is the half the bias multiplied by q (elemntary charge). Jouthc method returns 2 variables, current flux and boolean value telling solution is found or not. Note that the solar cell is under 10 suns concentration. [Short circuit current](https://www.pveducation.org/pvcdrom/solar-cell-operation/short-circuit-current) of typical solar cells is about 700 A/m<sup>2</sup> under 1 sun concentration.

To show hot carrier solar cell power output (W/m<sup>2</sup>) at 0.2V bias,
```python
>>> print("Pout:{:.3f}(A/m^2)".format(scell.Pouthc(0.2/2*nu.eV)[0]))
Pout:1020.160(W/m^2)
```

To show the max power under the given condition, use maxPouthc() method. First it will search through to find open circuit voltage Voc, then it will tell the max power.

```python
>>> print("MaxPout:{:.3f}(W/m^2)".format(scell.maxPouthc()[0]))
Solving Voc
Count:0 Tc:551.68 mue(eV):0.000 Bias(V):0.000 J(A/m^2):6786.75 P(W/m^2):0.00
Count:1 Tc:556.21 mue(eV):0.002 Bias(V):0.004 J(A/m^2):6786.75 P(W/m^2):27.15
.......
Count:59 Tc:300.62 mue(eV):0.118 Bias(V):0.236 J(A/m^2):738.11 P(W/m^2):174.19
Count:60 Tc:300.61 mue(eV):0.120 Bias(V):0.240 J(A/m^2):-195.59 P(W/m^2):-46.94
Stop iteration at mue(eV):1.200000e-01
maxP (W/m^2):1083.2
```

If you want to plot current voltage (IV) characteristics, you can try using test module (test.py). Use IV function. To plot power voltage (PV) characteristics, use PV function. Note this is under 10 suns concentration.

```python
>>> import test
>>> test.IV(scell)
>>> test.PV(scell)
```

<p align="center">
<img src="fig\IV_InAs.png" width=48%><img src="fig\PV_InAs.png" width=48%>
</p>

To find optimized extraction energy level at energy selective contact, use maxPouthc_opt_ESC function. It will scan through max powers and show the list of max powers and carrier temperatures. The second parameter tells number of trials between 50% band gap energy and 300% band gap energy. It will also scan energy width (currently set at 1 meV and 5 meV). If carrier extraction is via thermionic emission, it won't affect.

```python
>>> test.maxPouthc_opt_ESC(scell, 50)
cnt1:0/9 cnt2:0/1
deltaEesc:0.35(eV) Ewesc:5.000e+00(meV)
Solving Voc
Count:0 Tc:311.94 mue(eV):0.000 Bias(V):0.000 J(A/m^2):6786.15 P(W/m^2):0.00
Count:1 Tc:311.06 mue(eV):0.002 Bias(V):0.004 J(A/m^2):6786.05 P(W/m^2):27.14
......
Count:65 Tc:300.63 mue(eV):0.130 Bias(V):0.260 J(A/m^2):0.06 P(W/m^2):0.02
Count:66 Tc:300.63 mue(eV):0.132 Bias(V):0.264 J(A/m^2):-0.01 P(W/m^2):-0.00
Stop iteration at mue(eV):1.320000e-01
less than before. skip optimization
Max Power (W/m^2) [[1.10025840e+03]
 [1.11097004e+03]
 [1.06870710e+03]
 [7.05822553e+02]
 [6.85471318e+00]
 [1.26748304e-01]
 [0.00000000e+00]
 [0.00000000e+00]
 [0.00000000e+00]]
Carrier Temperature (K) [[ 303.30276712]
 [ 302.87446135]
 [ 301.47017415]
 [1756.70678502]
 [ 300.65141722]
 [ 300.62030558]
 [   0.        ]
 [   0.        ]
 [   0.        ]]
maxP:1110.97 (W/m^2)
T:302.874(K)
dEesc:0.55(eV) Ewesc:5.000e+00(meV)
```
Note that band gap energy is 0.354eV (InAs). The optimization found that the extraction energy is 0.55eV wide (or 0.098eV above the band edge) for thermionic emission (default).

## Classes and functions

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
    Vochc(self, fmaxP=False, dsp_msg=True)
        solve open circuit voltage Voc
    Pouthc(self, mue)
        return power output at given mue
    maxPouthc(self)
        solve max power at given Eg, ESC E and Ew
    """

class Photon_in(object):
    """ defining spectrum
    default spectrum is AM15 if n=2

    attributes
    ----------
    c : float
        concentration factor
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
            left electrical contact defined in a separate file
        rcnt : class
            right electrical contacts defined in a separate file
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

## File dependency

```
hcsc.py <-
    hcscAttribute.py
    recombination.py
    photon_in.py
    nu.py
    myutil.py

hcscAttribute.py <-
    bulk.py
    contacts.py
    extraction_via_ESC.py
    nu.py
```

## Background

Hot carrier solar cell is a third generation solar cell aiming to overcome Shockley-Queisser limit by utilizing excess energy above band gap which is normally lost in standard solar cells.

Ross and Nozik [1], as well as Wurfel [2] suggested the hot carrier solar cell (see figure below). The structure introduces the energy selective contact (ESC) in which carriers in the absorber are extracted at elevated energy E<sub>esc</sub>, e.g. electrons are extracted above the conduction band edge. The ESCs are located between the absorber and the regular contacts. The extraction energy has an energy band, w<sub>esc</sub>. In Ross and Nozik model, it is assumed to be a single energy level. ESC needs to be carrier selective. In the figure, electrons are extracted in the right hand side contact at energy E<sub>esc</sub> and holes are extracted at energy -E<sub>esc</sub> in the left hand side contact.

Chemical potential, or quasi-Fermi level in the absorber is &mu;<sub>c</sub> and that in contacts is denoted as &mu;<sub>e</sub>. If the carriers are extracted before thermalization to the band edges but after electron electron scattering which is super fast process, the carrier distributions can be described as heated Fermi-Dirac distribution with its own carrier temperature T<sub>c</sub>, which is higher than the room temperature T<sub>rm</sub>. For hot carrier solar cells to be effective, it is desirable to have narrow band gap materials (although Wurfel suggests Auger recombinations are imminent). In Ross and Nozik model, non-radiative recombinations are not considered. The figure bottom right shows the upper limit derived in Ross and Nozik model with various carrier temperatures (very artificial rersults, nevertheless). If the hot carriers are significantly utilized, it can go beyond the Shockley-Queisser limit.

<img src="fig\SchematicHCSC.png" width=48%><img src="fig\Eff_Eg_T_1000_sun.png" width=48%>

Here are the two equations to solve.

#### Particle balance equation

Photon flux from the sun will be absorbed with quantum efficiency of 1 above energy band gap, and 0 below the band gap. Because there are always radiative recombinations, some photon flux will escape described by generalized Planck law, assuming one electron hole pair will generate one photon. Without considering any other non-radiative path, net carriers should be extracted through energy selective contacts.

#### Energy balance equation

Energy extracted through energy selective contacts should balance to the amount of energy absorbed subtracted by radiative recombinations and carrier thermalization. This is what is described in the equation below.

<p align="center">
<img src="fig\particle_balance.PNG" width=70%>
<img src="fig\energy_balance.PNG" width=70%>
</p>

**References:**

[1] Robert T. Ross and Arthur J. Nozik. Efficiency of hot-carrier solar energy converters. Journal of Applied
Physics, 53(5):3813{3818, 1982.

[2] Peter Wurfel. Solar energy conversion with hot electrons from impact ionization. Solar Energy Materials
and Solar Cells, 46(1):43{52, 1997.

[3] Yasuhiko Takeda, Akihisa Ichiki, Yuya Kusano, Noriaki Sugimoto, and
Tomoyoshi Motohiro. Resonant tunneling diodes as energy-selective
contacts used in hot-carrier solar cells. Journal of Applied Physics,
118(12):124510, 2015.

[4] Steven Limpert, Stephen Goodnick, Christiana Honsberg, Gavin Conibeer,
and Stephen Bremner. A hot carrier solar cell device model using a coupled
electron phonon energy balance model. In 2013 IEEE 39th Photovoltaic
Specialists Conference (PVSC), pages 1054{1059. IEEE, 2013.

and more. Let me know if you need more references.

## Limitations

It will solve both equations decently. But sometimes solution is not so smooth as can be seen in IV and PV curves.

## Speeding up

3 pyx files and setup.py file are ready for Cython. It will cut the process time in half.

If anyone can figure out use of scipy.LowLevelCallable to speed up scipy.integrate.quad, it will be great.

Solver can be improved. For example, considering differential in each points, it can skip lots of computation along flat IV curve.


## Prerequisites

Confirmed working with Python 3.7 on Windows 10

It also require:

* Numpy
* Scipy
* MatPlotLib.
* ASTMG173.csv (available by downloading from [NREL web site](https://rredc.nrel.gov/solar//spectra/am1.5/). "Text Files in Compressed Format" is the one needs to be unpacked.)


## Installing

Download files into a folder.
