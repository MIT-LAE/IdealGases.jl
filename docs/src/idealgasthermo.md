# Thermodynamic properties of ideal gases

## Calculating ``c_p, h, s``
The ``c_p, h, s`` of the gas is calculated using the following `NASA-9` polynomials:

$$\begin{aligned}
\frac{\hat{c}_p^\circ(T)}{\overline{R}} &=
&& a_1T^{-2}
&+ &a_2 T^{-1} 
&+ &a_3 
&+ &a_4 T 
&+ &a_5T^2 
&+ &a_6 T^3 
&+ &a_7 T^4 &&&&\\
   
\frac{\hat{h}^\circ (T)}{\overline{R} T} &= 
&-&a_1T^{-2}
&+ &a_2\frac{\ln(T)}{T}
&+ &a_3 
&+ &a_4\frac{T}{2} 
&+ &a_5\frac{T^2}{3} 
&+ &a_6\frac{T^3}{4} 
&+ &a_7\frac{T^4}{5} 
&+ &\frac{b_1}{T}&&\\

\frac{\hat{s}^\circ(T)}{\overline{R}} &= 
&-&a_1\frac{T^{-2}}{2}
&- &a_2T^{-1} 
&+ &a_3 \ln(T) 
&+ &a_4 T 
&+ &a_5 \frac{T^2}{2} 
&+ &a_6 \frac{T^3}{3} 
&+ &a_7 \frac{T^4}{4} &&
&+ &b_2 
 \end{aligned}$$

where the $^\circ$ and $\hat{\ }$ denote standard conditions 
($T_\mathrm{std} = 298.15\, \mathrm{K},\; P_\mathrm{std}=101.325\, \mathrm{kPa}$). These are dropped in the subsequent 
documentation for convinence.

```@docs
IdealGases.Cp
IdealGases.h
IdealGases.𝜙
```
## Thermodynamic derivatives
Additionally it is also useful to calculate the following derivatives,

$$\begin{aligned}
\frac{d h}{d T} & = c_p\\
\frac{ds}{dT} = \frac{d\phi}{dT} & = \frac{c_p}{T}\\
\frac{d c_p}{d T} & = R_\mathrm{univ}\left(-2a_1T^{-3} - a_2 T^{-2} 
+ a_4 + 2a_5T + 3a_6T^2 + 4a_7T^3\right)\\
\end{aligned}$$

where the last derivative $\displaystyle {\frac{dc_p}{dT}}$ is obtained by differentiating 
the polynomial representation of $c_p$.