# Uncertain Differential Equations

## Abstract

A visual demonstration of chaotic behaviour in the simple model of a
Lorenz attractor system. Implemented in Python using NumPy and VTK.
Two strategies of evolution are employed. First, ordinary evolution of a
cloud of points packed around some point. Second, resampling of the cloud
of points after time dt. At times dt ≤ 0.01 unexpected behaviour becomes
apparent.

## Introduction

Chaos is usually described as the inability to predict the future behaviour of a deter-
ministic system after some time because of uncertain initial conditions or uncertainties
in the parameters of the underlying differential equation. Thus, it seems that the out-
come of a series of identical experiments is random. In contrast, random behaviour is
characterized by the property that absolutely the same experiment can give different
outcomes with some probability distribution. However, chaotic behaviour can seem
and can share some of the properties of intrinsically random systems.
Examples of intrinsically random systems are the quantum laws, where the best one
can hope to know for an observable quantity is some probability distribution. Other
systems which exhibit random behaviour usually are influenced by great many factors
without any clear relations between them. Examples of this are stock prices and many
processes in biology.
2 The Lorenz Attractor
The system of three ordinary differential equations called the Lorenz attractor were
first published in 1963 by the meteorologist and mathematician from MIT Edward
Lorenz. They were obtained as an approximation to the convective motion governed
by partial differential equations in fluid mechanics of a fluid cell which is warmed from
below and cooled from above. The set of nonlinear ODE can be written as:
```
dx
--  = σ(y − x)
dt

dy
-- = x(τ − z) − y
dt
dz
-- = xy − βz
dt
```

The values of the three parameters σ, τ, β determine the qualitative behaviour of the
system. To get chaotic effects we need τ as big as 28.0. For the purposes of my
simulation I will use σ = 10.0, τ = 28.0, β = 8/3.

## Simulation

In order to trace the sensitivity of the solution on perturbations in the initial data we
will adopt two approaches.

### Without Sampling

We pick up N points with normal distribution N (μ_0 , σ^2 I), i.e. the points are uniformly
distributed around μ 0 with variance σ 2 . Then we integrate the system for each initial
point up to some prescribed time T .

### With Sampling

The second procedure is as follows:
1. pick up the points as before – N points with normal distribution N (μ_0 , σ^2 I).
2. let the system evolve for some prescribed time dt.
3. assuming the cloud of points has normal distribution, we estimate the mean μ
and the covariance matrix Σ.
4. sample a new normal distribution N (μ, σ^2 I).
5. if the elapsed time is less than T , go to step 2.
