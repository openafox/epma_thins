---
layout: post
title: Questions to answer to help finish port
tag: help
category: Info

excerpt: Some questions I wrote down at one point
---



# In PAP

1. [1]p35e7 need to check refs for full explanation of where this comes from {PAP_41}
  ~~~ python
  javg += (el_i.c1 * el_i.z / el_i.mass * np.log(el_i.z*(10.04 + 8.25*np.exp(-1.0*el_i.z/11.22))/1000.))
  ~~~
2. Total Trajectory [g/cm^2]  [1]p35e9
  ~~~ python
  dr0 += 1/mavg*(j_mip**p2*d_k[k])*(el_p.e0**p1 - el_p.xray**p1)/p1
  ~~~
  updated? with E0^p1-Eq^p1 from just Eo^p1 need to find ref {PAP_62}
3. Ionization crosssection [1]p36e10-11
  mparam from Hutchins (get Ref) {PAP_67}
4. [1]p54e39 PAP Weighting laws. with extra factors (/5, /2, /3) where do they come from?  {PAP_89}
  See Scanning 1990 PAP p 218 as well…   N is integral..??



