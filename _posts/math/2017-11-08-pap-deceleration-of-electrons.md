---
layout: post
title:  Deceleration of electrons for the PAP
tag: test
category: Math

excerpt: Checking
---
The following comes from J.-L. Pouchou and F. Pichoir,
“Quantitative analysis of homogeneous or stratified microvolumes applying
the model ‘PAP,’” in Electron probe quantitation, Springer, 1991, pp. 35–36.


Bethe Law (used for ZAF)(good for high energy ~ > 30keV):

Get bethe paper to confirm eq - from wiki:

$$ \frac{dE}{ds} = -\frac{4\pi z^2}{m_e\nu^2}
\left(\frac{N_A Z \rho}{A}\right)
\left(\frac{e^2}{4\pi \epsilon_0}\right)^2
\cdot ln\left(\frac{\sqrt{2.718/2}E}{J_i}\right)$$

Combinining the constants and adjusting units to get keV/cm:

$$ \frac{dE}{ds} \mathrm{\left[\frac{keV}{cm}\right]}
= -C\left(\frac{Z \rho\mathrm{[g/cm^3]}}
{E\mathrm{[keV]}A\mathrm{[g/mol]}}\right)
\cdot ln\left(\frac{1.166E\mathrm{[keV]}}{J_i\mathrm{[keV]}}\right)$$

Where:

$$C\mathrm{\left[\frac{keV^2\cdot cm^{2}}{mol}\right]}
= \frac{4\pi N_A\mathrm{[mol^{-1}]}(e\mathrm{[C]})^4}
{32\pi^2 (\epsilon_0\mathrm{[C^2V^{-1}m^{-1}]})^2}
\cdot
\left(\frac{(100\mathrm{[cm/m])^2}}
{(e\mathrm{[V/eV]})^2(1,000\mathrm{[eV/keV])^2}}\right)$$
$$=78454\mathrm{\left[\frac{keV^2\cdot cm^{2}}{mol}\right]}$$

<!--((1.602e-19)**2*(6.023e23)/(8*np.pi*(8.854e-12)**2))/100-->

Which turns into:

$$ \frac{dE}{d \rho s}
\mathrm{\left[\frac{keV \cdot cm^2}{g}\right]}
= - \frac{C}{E}
\cdot \sum_i \frac{C_iZ_i}{A_i\mathrm{[g/mol]}} \cdot
ln\left(\frac{1.166E\mathrm{[keV]}}
{J_i\mathrm{[keV]}}\right)
\tag{4}\label{bethe}$$

$$J_i$$ is the mean ionization potential of each constituent

At low energies eq.\ref{bethe} produces lower results then experimental
penetrations.

At very low energies it is also desirable to have the log term be able to
change sign.

So we represent the energy term as $$1/f(V)$$ where $$V=E/J)$$:

<!--PAP Method-->
$$ \frac{dE}{d \rho s} \mathrm{\left[\frac{keV \cdot cm^2}{g}\right]}
= - \frac{M}{J}\cdot
\frac{1}{f(V)}
\tag{5}\label{eq5}$$

Where:

$$M\mathrm{\left[\frac{1}{g}\right]}
=\sum_i\frac{C_iZ_i}{A_i}\tag{5.1}$$

and

$$ln(J) = \sum_i C_i\frac{Z_i}{A_i}ln(J_i)/M
\tag{6}$$

The expression for $$J_i$$ was obtained from empirical fist of experiment
see refs:

$$J_i \text{[keV]}= 10^{-3}\cdot Z_i \cdot
\left[10.04+8.25exp\left(\frac{-Z_i}{11.22}\right)\right]
\tag{7}$$


$$f(V)=\sum_{k=1}^3 D_k \cdot V^{P_k}\tag{8}$$

with:

|$$D_1=6.6\times 10^{-6}$$ |$$P_1=0.78$$|
|$$D_2=1.12\times 10^{-5}(1.35-0.45J^2)$$|$$P_2=0.1$$
|$$D_3=2.2\times 10^{-6}/J$$|$$P_3=-(0.5-0.25J)$$

