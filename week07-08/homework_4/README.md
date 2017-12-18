Instructions

Note: &lt;X&gt; denotes the expectation of the random variable X which we denoted E(X) so far.

In this homework, you will grasp the intimate relation between sampling and integration and sense the full power of the Markov-chain Monte Carlo method, especially in high dimensions. Our goal will be to compute the volume of the sphere in 200 dimensions. This volume is smaller than

0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001...

Computing this volume is a sampling problem, and it literally corresponds to finding a needle in a haystack (or rather, a specific atom in all the observable universe). Nevertheless, we will succeed.

Two remarks before we start:

Often, one speaks of the "hypersphere in 200 dimensions", and we will also need the "hypercylinder". In this homework, we simply leave out the "hyper", and speak of "spheres" in d dimensions and "cylinders" in d dimensions. You will discover, among others, that the disk is a two-dimensional sphere, and that the square is a two-dimensional cylinder (see below). Of course there are precise mathematical definitions. In particular, we consider unit spheres and unit cylinders.
In one dimension, normally, one speaks of a "length", in two dimensions, we have an "area", and in three dimensions, there are "volumes". In this homework, we speak of "volume" for all these measures. A "volume", in our definition, has no units (like cm^3) and, although it might sound strange to some of you, we may compare volumes across different dimensions.

Some mathematical definitions, and some discussion:

Unit sphere in d dimensions

Definition: d-dimensional object of all points (x_0, ..., x_{d-1}) such that (x_0^2 + .... + x_{d-1}^2) < 1.

Special cases, and general case for unit spheres

In d=1, the unit sphere is the segment \[-1,1\] and its volume (or rather the length of the segment) is equal to 2.
In d=2, the unit sphere is equal to the disk of radius 1, and its volume (or rather its area) is equal to pi
In d=3, the unit sphere has volume 4/3 pi, as you can look up in your old Math books from .
In general d, we call the volume of the unit sphere V_sph(d).

Unit cylinder in d dimensions

Definition: d-dimensional object of all points (x_0, ..., x_{d-1}) such that (x_0^2 + .... + x_{d-2}^2) < 1 and -1 < x_{d-1} < 1.

The unit cylinder in d dimensions is put together from a unit sphere in d-1 dimensions, and a segment of length 2 along one dimension.

Special cases, and general case for unit cylinders

In d=1, the unit cylinder does not exist.
In d=2, the unit cylinder is equal to the square with x and y between -1 and 1. Its volume is equal to 4.
In d=3, the unit cylinder is put together from a disk (area pi) and a segment of length 2, its volume is equal to 2 pi. Again you can look this up in your old math books, or on wikipedia.
In general d, we call the volume of the unit cylinder V_cyl(d).

A crucial table containing exact volumes and relations between volumes in different dimensions in file sphy_size.png.

## A

In this section, we interpret the table and the mathematical definitions of the introduction. We strongly suggest that you print out the definitions and the table or copy them onto paper, to have a vade mecum (cheat-sheet) for this entire homework session.

## A1

Explain, using the mathematical definition given in the introduction, why the segment \[-1, 1\] represents a one-dimensional sphere, and why the square from -1 to 1 in x and from -1 to 1 in y is a unit cylinder in two dimensions (two sentences).
Confirm and explain that the unit cylinder in two dimensions has twice the volume of the unit sphere in one dimension (two sentences).
Explain (prove) why the above relations

V_cyl(2) = 2 V_sph(1)

V_cyl(3) = 2 V_sph(2)

are valid in arbitrary dimensions, that is, prove that

V_cyl(d + 1) = 2 V_sh(d)

In the above table, the quantity Q(d) equals twice the sphere volume in d dimensions divided by the unit cylinder volume in d dimensions. Prove that Q(d) equals the ratio of the sphere volume in dimension (d) and the sphere volume in dimension (d-1), that is, prove that

Q(d) = V_sph(d) / V_sph(d-1)

NB: In the following, we will compute Q(d) by Markov-chain Monte Carlo sampling.

## A2

In this section, you will write a simple Markov-chain Monte Carlo program to compute the quantity Q(d=3). To this aim,

download the program markov_pi.py (from week2 directory of the repository), and modify it such that it samples points (x, y) inside the unit disk, rather than inside the square. Start the simulation at the origin, that is, the point (0.0, 0.0), rather than at (1.0, 1.0).
At each step, sample a new value of z as random.uniform(-1.0, 1.0), and count as a "hit" if x^2 + y^2 + z^2 < 1.0.

When your program is ready, print it.

Compute the average value &lt;Q_3&gt; = 2 N_hits / N_trials, and check that it is close to 4/3, that is the ratio of the sphere volume for d=3 to the sphere volume for d=2 (see table). Run your program for as long as possible, to be sure that the results come out OK. Compare your result for &lt;Q_3&gt; with the value 4/3. It should be extremely close.
 
## A3

Write a new program, a modification of the program you wrote in Section A2 (itself a modification of markov_pi.py), to compute the volume of the unit sphere in four dimensions. To do so:

Generalize the program of Section A2 (by simply changing a few lines) so that it samples points (x,y,z) inside the three-dimensional unit sphere.
At each step, sample a alpha-value as random.uniform(-1.0, 1.0), and count as a "hit" if (x^2 + y^2 + z^2 + alpha^2) < 1.
Compute the average value of Q_4, as you did for Q_3 in section A2.

Then answer the following questions:

- Print your program, and verify that the simulation yields the correct estimate of Q_4.
- The volume of the four-dimensional unit sphere is equal to pi^2/2. What does this imply for the value of Q_4? (answer without any programming, using the table)
- Check that the volume V_sph(d=4) of the four-dimensional unit sphere can be given by three equivalent formulas

V_sph(d=4) = pi^2 / 2

V_sph(d=4) = V_sph(3) \* Q_4, approximately equal to: V_sph(3) \* &lt;Q_4&gt;

V_sph(d=4) = V_sph(2) \* Q_3 \* Q_4, approximately equal to: V_sph(2) \* &lt;Q_3&gt; \* &lt;Q_4&gt;

## B

In this section, you will write a Markov-chain program, generalization of markov_pi.py, that samples reasonably well in any dimension d.

## B1

Generalize the program markov_pi.py to implement an efficient Markov-chain Monte Carlo algorithm to sample uniformly distributed points inside the d-dimensional unit sphere, using the following hints:

Represent the configuration point as a list: x = \[x_0, x_1,...., x_k, ..., x_{d - 1}\]

Start the simulation at the origin \[0.0, 0.0, 0.0, ... 0.0\] (you can program this as: x = \[0\] * d), as in Sections A2 and A3.

Instead of modifying all components of x at a time, as we did in markov_pi.py, modify only one component at each iteration i (with i=0, 1, 2,...., n_trials). This can be done as in:

```python
k = random.randint(0, d - 1)
x_old_k = x[k]
x_new_k = x_old_k + random.uniform(-delta, delta)
```

Then you should accept the move if the new radius is <1, and reject otherwise (remaining in the same configuration).

Use an optimized way of computing the new radius (useful when d is large):

```python
new_radius_square = old_radius_square + x_new_k ** 2 - x_old_k ** 2
```

Do not forget to initialize correctly the variable old_radius_square, and to update it as needed.

Once the code is ready:

Print your program, which works for a general d.


## B1 (continued)

Test your program for d=4, by plotting the normed histogram of r = sqrt(x\[0] \*\* 2 + x\[1] \*\* 2 + ... + x\[3] \*\* 2). This can be done by the command pylab.hist(..., normed=True). Also plot the analytic curve P(r)=4 r^3 (0 < r< 1) on the same graph, and compare it to your histogram. Print the graphics file with this comparison.

## B1 (continued)

Test your program for d=20, again by plotting the histogram of r and plotting on the same graph also the analytic function P(r)=20 r^19 (0<r<1). Print the graphics file with this comparison.

## B2

Now,

- introduce the calculation of &lt;Q&gt; into the program you wrote in Section B1, on the line of what was done in Sections A2 and A3.
- Use your program to estimate the values of Q(4)=V_sph(4)/V_sph(3) and of Q(200)=V_sph(200)/V_sph(199)
- Compare your result with the analytical formula that you can obtain from the program

```python
import math

def V_sph(dim):
    return math.pi ** (dim / 2.0) / math.gamma(dim / 2.0 + 1.0)

for d in range(1, 20):
    print d, V_sph(d)
```

Once your code is ready,print your program.

## B2 (continued)

- Give your result for V_sph(4)/V_sph(3), as precisely as possible. Compare it with the number obtained in Section A3.
- Give your result for V_sph(200)/V_sph(199).

## C

In this section, you compute the actual volume V_sph(200) of the 200-dimensional sphere, from a series of calculations. Finally you do some error estimations.

## C1

- Modify your program of section B so that it loops over dimension d=1,2,3,4,5,6,7,..., d_max, and for each dimension d, compute the average value &lt;Q(d+1)&gt;. Each iteration should take a constant number of iterations, independent of d.
- Use the values of Q(d+1) which you are computing iteratively, to evaluate V_sph(4) once more (start from the value V_sph(1)=2). Write down your result, and compare it with pi^2 / 2.
- Use this program to compute, in one run, V_sph(200) = 2 * &lt;Q_2&gt; * &lt;Q_3&gt; * &lt;Q_4&gt; * .... * &lt;Q(200)&gt; (NB: the initial factor of 2 corresponds to the volume of the d=1 sphere: V_sph(1)=2). To do that, before switching from d to (d+1) dimensions, use the value of Q(d+1) to obtain V_sph(d+1) as a product of V_sph(d) by Q(d+1). At the end, write down your result for V_sph(200)
- print this program, which allows you to compute V_sph for any dimensionality.

## C1 (continued)

- Use this program to plot all the values of V_sph(d) as a function of d for d=1,....,200 (NB: use a logarithmic scale for the y-axis with pylab.yscale('log'), or similar). On the same plot, show the analytical curve that was given in Section B2. Print the graphics file containing the plot of the analytic formula and of the Monte Carlo result.

## C2

In every Monte Carlo calculation, we should control the error. To save time, you will do this for d=20, rather than for d=200.

- Modify the program from Section C1, by adding an external cycle which loops over values of n_trials = 1, 10, 100, 1000, ...
- For each choice of n_trials, perform n_runs independent runs (take for instance n_runs=10) of the calculation of V_sph(20). For each run, always start from the origin as initial condition.
- Compute the averages &lt;V_sph(20)&gt; and &lt;V_sph(20)^2&gt; for the n_runs runs, and estimate the error on the average as sqrt(&lt;V_sph(20)^2&gt; - &lt;V_sph(20)&gt;^2)/sqrt(n_runs)
- Print your modified program.
- Produce the following table:

```
n_trials | <V_sph(20)> |  V_sph(20) (exact) | error | difference
10 ...
100 ...
1000 ...
```

 where the last column corresponds to the difference between the Monte Carlo and exact results for V_sph(20).

Comment: are the error bars consistent with the difference from the exact result, for very short Markov-chain calculations?
