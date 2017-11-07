---
layout: post
title:  "Welcome to Jekyll!"
date:   2017-11-07 13:26:30 -0800
categories: jekyll update
---

##checking the math in PAP eq39

Treating each electron individually and neglecting absorption effects we can
estimate the thickness from eq 8.

We must also be assuming an element is
only contained in 1 layer:

$$p(\rho z)=N\cdot (\rho z-L)^2\cdot (\rho z-R)^2\label{eq39}\tag{39}$$

$$= N\cdot (\rho z^2-2\rho zL+L^2)(\rho z^2-2\rho zR+R^2)$$

$$=N\cdot (\rho z^4-2\rho z^3R-\rho z^2R^2-2\rho z^3L+
4\rho z^2RL-2\rho zR^2L-\rho z^2L^2-2\rho zRL^2+L^2R^2)$$

$$=N\cdot (\rho z^4-2\rho z^3(R+L)+\rho z^2(R^2+4RL+L^2)-
           2\rho z(R^2L+RL^2)+L^2R^2)$$

with N such that:

$$\int^R_0 p(\rho z)d\rho z=1$$

reference to eq39: \ref{eq39},
