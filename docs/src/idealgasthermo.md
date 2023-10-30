# Thermodynamic properties of ideal gases

## Calculating ``c_p, h, s``
The ``c_p, h, s`` of a particular species is calculated using the following `NASA-9` polynomials:

$$\begin{aligned}
\cphatR &= %\frac{\hat{c}_p^\circ(T)}{\overline{R}} &=
&& a_1T^{-2}
&+ &a_2 T^{-1} 
&+ &a_3 
&+ &a_4 T 
&+ &a_5T^2 
&+ &a_6 T^3 
&+ &a_7 T^4 &&&&\\
   
\hhatRT &= 
&-&a_1T^{-2}
&+ &a_2\frac{\ln(T)}{T}
&+ &a_3 
&+ &a_4\frac{T}{2} 
&+ &a_5\frac{T^2}{3} 
&+ &a_6\frac{T^3}{4} 
&+ &a_7\frac{T^4}{5} 
&+ &\frac{b_1}{T}&&\\

\frac{\hat{\phi}^\circ(T)}{\Ru} &= 
&-&a_1\frac{T^{-2}}{2}
&- &a_2T^{-1} 
&+ &a_3 \ln(T) 
&+ &a_4 T 
&+ &a_5 \frac{T^2}{2} 
&+ &a_6 \frac{T^3}{3} 
&+ &a_7 \frac{T^4}{4} &&
&+ &b_2 ,
 \end{aligned}$$

where the $^\circ$ and $\hat{\ }$ denote standard pressure ($P_\mathrm{std}=101.325\, \mathrm{kPa}$) and molar basis respectively and ``\Ru = 8.3144\, \mathrm{J\, mol^{-1}\,  K^{-1}}`` is the universal gas constant. 

``\hat{\phi}^\circ``is the entropy complement function (we use ``\phi`` instead of ``s`` to emphasize that the entropy complement is only a function of temperature and not, in general, equal to the entropy) and is equal to ``\hat{s}^\circ`` only at standard pressure. At any pressure other than ``\Pstd`` the entropy ``\hat{s}`` is given by,

```math
\frac{\hat{s}(T)}{\Ru} = \phihatR - \ln{\frac{P}{\Pstd}}.
```


```@docs
IdealGases.Cp
IdealGases.h
IdealGases.𝜙
```

## [Representing mixtures](@id mixthermo)

The molar specific heat, enthalpy, and entropy of a mixture of ideal gases at standard pressure, ``\Pstd``, is given by,

```math
\begin{aligned}
\cpbarR &= \sum_{i=1}^n \Xi\left.\cphatR\right|_i\\
\hbarRT &= \sum_{i=1}^n \Xi\left.\hhatRT\right|_i\\
\phibarR &= \sum_{i=1}^n \Xi\left.\phihatR\right|_i + \Xi \ln{\Xi}.
\\

\end{aligned}
```
The term ``\Xi \ln{\Xi}`` in the entropy complement function represents the [*entropy of mixing*](https://en.wikipedia.org/wiki/Entropy_of_mixing) when multiple species are present in the gas mixture.

When ``P \neq \Pstd``, ``\hat{c}_p``, ``\hat{h}``, and ``\phi`` are as above but ``s`` has a pressure term.

```math
\begin{aligned}
\overline{\frac{\hat{c}_p (T)}{\Ru}} = \cpbarR &= \sum_{i=1}^n \Xi\left.\cphatR\right|_i\\
\overline{\frac{\hat{h} (T)}{\Ru T}} = \hbarRT &= \sum_{i=1}^n \Xi\left.\hhatRT\right|_i\\
\overline{\frac{\hat{\phi}(T)}{\Ru}} = \phibarR &= \sum_{i=1}^n \Xi\left.\phihatR\right|_i + \Xi \ln{\Xi} \\
\overline{\frac{\hat{s}(T)}{\Ru}} &= \sum_{i=1}^n \Xi\left.\phihatR\right|_i + \Xi \ln{\Xi} - \Xi \ln{\frac{P}{\Pstd}}.
\\

\end{aligned}
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