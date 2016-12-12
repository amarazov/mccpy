## A

In this section, you study the simulated-annealing method in simple cases. First, in Section A1, you observe the tricky consequences of automatic step-size control in a general Markov-Chain Monte Carlo program. In Section A2, you discover the essence of simulated annealing in a simple two-dimensional optimization problem, and in Section A3, you consider a simple non-generic example where simulated annealing breaks down, at least sometimes.

## A1

Consider the following normalized probability distribution prob(x), given by the sum of two Gaussians:

```python
def prob(x):
    s1 = math.exp(-(x + 1.2) ** 2 / 0.72)
    s2 = math.exp(-(x - 1.5) ** 2 / 0.08)
    return (s1 + 2.0 * s2) / math.sqrt(2.0 * math.pi)
```

Write a Markov-chain program that draws samples {x_0, x_1, x_2, x_3, ...} from prob(x) using the Metropolis algorithm, and compute the average of the variable x in the simulation:

(x_0 + x_1 + x_2 + ... + x_{nsteps-1}) / nsteps

which gives you an estimate of the true average value of x:

\int (x prob(x)) dx = -0.12

For your convenience, here is this program:

```python
import random, math

def prob(x):
    s1 = math.exp(-(x + 1.2) ** 2 / 0.72)
    s2 = math.exp(-(x - 1.5) ** 2 / 0.08)
    return (s1 + 2.0 * s2) / math.sqrt(2.0 * math.pi)

delta = 10.0
nsteps = 10000
acc_tot = 0
x = 0.0
x_av = 0.0
for step in xrange(nsteps):
    xnew = x + random.uniform(-delta, delta)
    if random.uniform(0.0, 1.0) < prob(xnew) / prob(x):
        x = xnew
        acc_tot += 1
    x_av += x

print 'global acceptance ratio:', acc_tot / float(nsteps)
print '<x> =', x_av / float(nsteps)
```

Download (copy/paste) this program. Run it a few times, modifying delta and increasing nsteps. Once you have found a good choice of delta, run your programs at least four times, and communicate your four estimates of x_av (also communicate the values of delta and nsteps that you are using, and the typical acceptance rate across the runs). Comment: are the MC estimates similar to the exact result (x_av=-0.12)? (NB: the full error analysis is not required, four independent results are sufficient).



## A1 (continued)

- Modify this program, so that it includes an additional variable acc_tmp, which is used to perform an automatic step-size control. This variable is increased by one unit every time that a move is accepted (in the same way as acc_tot). Then, every 100 steps, you should check if acc_tmp is inside the interval [40,60]: if acc_tmp is larger than 60, you should increase the step size (from the current delta to delta*1.1), and if acc_tmp is smaller than 40 you should decrease the step size (from the current delta to delta/1.1). After this check, don't forget to reset the acc_tmp variable to 0.
- Run you program with some "extreme" initial values of delta (for instance delta=0.001 or delta=100.0), and verify that the global acceptance rate is close to 50% (NB: the global acceptance rate is computed using acc_tot, not acc_tmp), which is a direct consequence of the step-size control. This is only a check, nothing is required for this point.
- Run your program at least four times, and communicate the values you find for x_av. Comment: are the MC estimates similar to the exact result (x_av= -0.12)? Can you explain why?

(Hint: a comparison of the histogram for the MC samples with the exact probability distribution is not required, but it can be useful to look at it to explain the results).

## A2

Consider the function V(x, y), defined on the usual heliport (x between -1 and 1, y between -1 and 1):

```python
def V(x, y):
    pot  = -4.0 * x ** 2 - x ** 3 + 4.0 * x ** 4
    pot += -4.0 * y ** 2 - y ** 3 + 4.0 * y ** 4
    return pot
```

Here we show a color plot of this function: plot_A2_contour.png

from which it is clear that V(x,y) has one global minimum (with x=y=0.807044513157..) and three local minima (located next to the other three corners of the square). Our goal is to find the global minimum with a Monte Carlo method, without remaining trapped in one of the other three local minima. This problem is quite difficult to address in higher dimensions, where there may be many local minima, and where you cannot draw a map to orient yourself.Download (copy/paste) the following program:

```python
import math, random

def V(x, y):
    pot  = -4.0 * x ** 2 - x ** 3 + 4.0 * x ** 4
    pot += -4.0 * y ** 2 - y ** 3 + 4.0 * y ** 4
    return pot

xmin, ymin = 0.807044513157, 0.807044513157
gamma = 0.4
n_runs = 10
n_success = 0
for run in range(n_runs):
    T = 4.0
    x, y = 0.0, 0.0
    delta = 0.1
    step = 0
    acc = 0
    while T > 0.00001:
        step += 1
        if step == 100:
            T *= (1.0 - gamma)
            if acc < 30:
               delta /= 1.2
            elif acc > 70:
               delta *= 1.2
            step = 0
            acc = 0
        xnew = x + random.uniform(-delta, delta)
        ynew = y + random.uniform(-delta, delta)
        if abs(xnew) < 1.0 and abs(ynew) < 1.0 and \
           random.uniform(0.0, 1.0) < math.exp(- (V(xnew, ynew) - V(x, y)) / T):
            x = xnew
            y = ynew
            acc += 1
    if math.sqrt((x - xmin) ** 2 + (y - ymin) ** 2) < 0.1:
        n_success += 1
print gamma, '\t', n_success / float(n_runs)
```

and study it carefully. You should recognize its different elements:

- n_runs independent minimization runs are performed, and for each one there is a check of whether the global minimum was found (this is encoded in the n_success counter),
- the program includes an annealing schedule, in which the temperature is reduced every 100 iterations: T -> T * (1 - gamma),
- each run is stopped when a certain final temperature is reached,
- an automatic step-size control is implemented, to keep the acceptance ratio in the range [0.3,0.7] (see "if acc < 30" and following lines).

These are common features for a simulated-annealing program, and you will recognize them also in sections B and C.

- Run this program for at least six values of gamma between 0.0025 and 0.6. For each choice of gamma, compute the approximate success rate, i.e. the probability with which the correct solution is found (if you have enough CPU time available, you can increase n_runs for more precise results).
- Supply a table of your results (the values of gamma and the corresponding success rates).
- Comment: is it true that for a finite annealing rate gamma, the correct solution is always found? What is the tendency of the success rate for gamma going to zero?

## A3

We now imagine a "Gedankenexperiment" (thought experiment) of a hard disk in the container shown below. An actual implementation is neither required nor expected. The container has two boxes connected through a hard bottleneck. All its contours, and in particular the bottleneck, are made of hard walls.
see container.png

We imagine solving the problem of finding the radius of the largest disk that can be placed into this container. To do so, we perform simulated annealing for one hard disk in this container, with displacements (x,y) -> (x + delx, y + dely), where delx, dely (that can be positive or negative) are much smaller than the dimensions of the container.

After each 100 iterations, the disk radius is increased

```
 sigma -> sigma * (1 + gamma)  with gamma > 0 (simulated annealing step)
```

if this is possible (if there is an overlap, no action is taken). Our goal is to find the largest disk fitting into the container.

- Convince yourself that simulated annealing for this case does not find with probability 1 the optimal solution (disk to the right), even for gamma -> 0 in the limit t -> infinity. Write down a short explanation of why this is so.

## B

Here you perform simulated annealing for the disk-packing problem. For simplicity, please download (cut-and-paste) the below program.

```python

import random, math

def unit_sphere():
    x = [random.gauss(0.0, 1.0) for i in range(3)]
    norm =  math.sqrt(sum(xk ** 2 for xk in x))
    return [xk / norm for xk in x]

def minimum_distance(positions, N):
    dists = [math.sqrt(sum((positions[k][j] - positions[l][j]) ** 2 \
             for j in range(3))) for l in range(N) for k in range(l)]
    return min(dists)

def resize_disks(positions, r, N, gamma):
    Upsilon = minimum_distance(positions, N) / 2.0
    r = r + gamma * (Upsilon - r)
    return r

N = 13
gamma  = 0.5
min_density = 0.78
for run in range(10):
    print 'run', run
    sigma  = 0.25
    r = 0.0
    positions = [unit_sphere() for j in range(N)]
    n_acc = 0
    step = 0
    while sigma > 1.e-8:
        step += 1
        if step % 500000 == 0:
            eta = N / 2.0 * (1.0 - math.sqrt(1.0 - r ** 2))
            print r, eta, sigma, acc_rate
        k = random.randint(0, N - 1)
        newpos = [positions[k][j] + random.gauss(0, sigma) for j in range(3)]
        norm = math.sqrt(sum(xk ** 2 for xk in newpos))
        newpos = [xk / norm for xk in newpos]
        new_min_dist = min([math.sqrt(sum((positions[l][j] - newpos[j]) ** 2 \
                       for j in range(3))) for l in range(k) + range(k + 1, N)])
        if new_min_dist > 2.0 * r:
            positions = positions[:k] + [newpos] + positions[k + 1:]
            n_acc += 1
        if step % 100 == 0:
            acc_rate = n_acc / float(100)
            n_acc = 0
            if acc_rate < 0.2:
                sigma *= 0.5
            elif acc_rate > 0.8 and sigma < 0.5:
                sigma *= 2.0
            r = resize_disks(positions, r, N, gamma)
            R = 1.0 / (1.0 / r - 1.0)
            eta = 1.0 * N / 2.0 * (1.0 - math.sqrt(1.0 - r ** 2))
    print 'final density: %f (gamma = %f)' % (eta, gamma)
    if eta > min_density:
        f = open('N_' + str(N) + '_final_'+ str(eta) + '.txt', 'w')
        for a in positions:
           f.write(str(a[0]) + ' ' + str(a[1]) + ' ' + str(a[2]) + '\n')
        f.close()
```


Study this program carefully, and compare it with the program in Section A2.

- Several independent runs are performed, and each time one solution is found.
- The annealing schedule consists in increasing the disk radius (which corresponds to decreasing the temperature). The parameter gamma constitutes an annealing rate (gamma=1 gives a fast annealing, while a small gamma makes the annealing slower).
- Automatic step-size control is used (by changing sigma) to keep the acceptance ratio in a reasonable interval.
- The relation between outer-sphere radius R and disk radius r reads R = 1 / (1 / r -1).
- The parameter eta = N * (1 - sqrt(1 - r^2)) / 2 is the surface area of the spherical caps formed by the disks.
- For each run, the final configuration is written on a file. This is done only for configurations that are above a given density min_density (to avoid storing configurations which are clearly not optimal).

## B1

- Run the simulated annealing program for N=13 disks on a sphere. Experiment with its various parameters, especially with the annealing rate gamma, and the control parameter sigma, that sets the step width of the Markov chain.
- Use somewhat smaller values of gamma to recover the optimal solution for N=13, that has density eta=0.79139, and a sub-optimal solution with density eta=0.78639. Communicate the parameter gamma you found useful. Upload both configurations (optimal and sub-optimal), i.e. the x,y,z positions of all spheres, which are written on a file by the given program. If you want to produce a 3D plot of the configuration (not required), see example_pylab_visualization.py and direct_sphere_disks_movie.py, from Week 9 programs (notice that the second one requires the mayavi library, see http://docs.enthought.com/mayavi/mayavi).
- Answer the following question: For small annealing rates, do you always recover the optimal solution (see discussion in section A2) or do you keep finding different solutions even for small annealing rates (case discussed in section A3)?

NB: For N=13, it has been proven in 1951 that R has to be lower than 1 (Schütte and van der Waerden), but it took until 2012 to show that the eta=0.79139 solution is the optimal one (Musin and Tarasov).

## B2

- Use the program of Section B1 to find the optimal solution for N = 15, which should have density eta = 0.80731. Upload the configuration (again, as a list of positions).
- In fact, for N=15, prove that there are TWO optimal solutions at the same density. To do so, write a program that checks whether the 5-connected disks touch each other, and use this criterion to distinguish the two solutions.

## B3

- Use the program of Section B1 to find the optimal solution for N=19, which should have density eta = 0.81096. Upload the configuration (again, as a list of positions).
- Optional (zero points to be gained): Show that there are infinitely many optimal solutions, as one of the 19 disks is not jammed. Write a program checking that one of the 19 disks is zero-connected: this "rattler" is not in contact with any other disk.

NB: Configurations of Sections B2 and B3 were found by D. A. Kottwitz (1991) using a very complicated optimization program.\

## C

Here you consider the travelling salesman problem (TSP), one of the classic problems in combinatorial optimization: given N cities, with distances d(i,j) = d(j,i) between them, find the shortest closed tour visiting all of them. The TSP is a prominent example of the class of NP complete problems, for which it is very difficult to find a solution, but this solution is easy to check. We will apply simulated annealing to the TSP in two dimensions. Simulated annealing provides a great "first approach" to this problem.

## C1

Before simulated annealing, let us try a direct-sampling approach. For this, download (cut-and-paste) the below program.

```python
import random, math, pylab

def dist(x, y):
    return math.sqrt((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2)

def tour_length(cities, N):
    return sum(dist(cities[k + 1], cities[k]) for k in range(N - 1)) + dist(cities[0], cities[N - 1])

N = 10
random.seed(54321)
cities = [(random.uniform(0.0, 1.0), random.uniform(0.0, 1.0)) for i in range(N)]
random.seed()
energy_min = float('inf')
for sample in xrange(1000000):
    random.shuffle(cities)
    energy =  tour_length(cities, N)
    if energy < energy_min:
        print sample, energy
        energy_min = energy
        new_cities = cities[:]
cities = new_cities[:]
for i in range(1,N):
    pylab.plot([cities[i][0], cities[i - 1][0]], [cities[i][1], cities[i - 1][1]], 'bo-')
pylab.plot([cities[0][0], cities[N - 1][0]], [cities[0][1], cities[N - 1][1]], 'bo-')
pylab.title(str(energy_min))
pylab.axis('scaled')
pylab.axis([0.0, 1.0, 0.0, 1.0])
pylab.savefig('plot_tsp_direct_sampling_N' + str(N) + '_energy' + str(energy_min) + '.png')
pylab.show()
```


Study this program. Notice that any tour is given by a random permutation of the cities. Notice also that the random number generator is initialized by a seed in order to generate an 'instance' (positions of the N cities), but is then randomized. Subsequent runs of the program produce the same instance, but then the program uses different random numbers at each run.

- Run the program at least four times for N=10, for 1 Million iterations each (you may change the 54321 seed, but keep it constant across the different runs), and upload one of the final tours you find, as a graphics file.

C1 (continued)

- Communicate the minimal lengths that you obtain. Comment: by looking at the optimal lengths in the four runs, and by looking at the figure, do you get the impression that you have found the global minimum of the problem?

C1 (continued)

- Repeat the same, with N=20: Run the program several times for 1 Million iterations, and upload one graphics file.
- communicate the different lengths, and comment.

C2

Now let us do simulated annealing, using as an energy the length of the tour. Download (cut-and-paste) the below program, that you are free to run and to modify.

```python
import random, math, pylab

def dist(x, y):
    return math.sqrt((x[0] - y[0]) ** 2 + (x[1] - y[1]) ** 2)

def tour_length(cities, N):
    return sum(dist(cities[k + 1], cities[k]) for k in range(N - 1)) + dist(cities[0], cities[N - 1])

N = 10
random.seed(54321)
cities = [(random.uniform(0.0, 1.0), random.uniform(0.0, 1.0)) for i in range(N)]
random.seed()
random.shuffle(cities)
beta = 1.0
n_accept = 0
best_energy = float('inf')
energy =  tour_length(cities, N)
for step in xrange(1000000):
    if n_accept == 100:
        beta *=  1.005
        n_accept = 0
    p = random.uniform(0.0, 1.0)
    if p  < 0.2:
        i = random.randint(0, N / 2)
        cities = cities[i:] + cities[:i]
        i = random.randint(0, N / 2)
        a = cities[:i]
        a.reverse()
        new_cities =  a + cities[i:]
    elif p < 0.6:
        new_cities = cities[:]
        i = random.randint(1, N - 1)
        a = new_cities.pop(i)
        j = random.randint(1, N - 2)
        new_cities.insert(j, a)
    else:
        new_cities = cities[:]
        i = random.randint(1, N - 1)
        j = random.randint(1, N - 1)
        new_cities[i] = cities[j]
        new_cities[j] = cities[i]
    new_energy =  tour_length(new_cities, N)
    if random.uniform(0.0, 1.0) < math.exp(- beta * (new_energy - energy)):
        n_accept += 1
        energy = new_energy
        cities = new_cities[:]
        if energy < best_energy:
           best_energy = energy
           best_tour = cities[:]
    if step % 100000 == 0:
        print energy, step, 1.0 / beta

cities = best_tour[:]
for i in range(1, N):
    pylab.plot([cities[i][0], cities[i - 1][0]], [cities[i][1], cities[i - 1][1]], 'bo-')
pylab.plot([cities[0][0], cities[N - 1][0]], [cities[0][1], cities[N - 1][1]], 'bo-')
pylab.title(str(best_energy))
pylab.axis('scaled')
pylab.axis([0.0, 1.0, 0.0, 1.0])
pylab.savefig('plot_tsp_simulated_annealing_N' + str(N) + '_energy' + str(best_energy) + '.png')
pylab.show()
```

Study this program carefully, and try to make connections with the program in section A2.

- The annealing schedule consists in increasing beta (i.e. decreasing the temperature) every 100 accepted moves.
- Three different moves are implemented, and proposed with probabilities 20%, 40%, 40%. Make sure you understand these three moves (the explanation is not required). Do not hesitate to experiment with different moves, and see whether they help in finding good solutions also for large N.
- The acceptance probability corresponds to the usual Metropolis condition, which makes use of the current value of beta.
- Run the program for the same instance as in C1, for N=10. Upload the graphics file of the solution you found.
- Do you find the same solution as in C1? Or even do you find a better solution?
- Run the program for the same instance as in C1, for N=20. Upload the graphics file of the solution you found.
- Do you find the same solution as in C1? Or even do you find a better solution?
- Comment: Do you think that the program finds the optimal solution or can you see, by visual inspection, that the solution found is sub-optimal? 



