## A

We first study the direct-sampling algorithm direct_pi.py and modify it somewhat.

## A1

Consider the program "direct_pi_multirun.py". Modify it so that it computes the root mean square (rms) deviation:

rms deviation `=\sqrt(1/n_runs \sum_{i=1}^{n_runs} (\pi^est_i−\pi)^2)`

, where \pi^est_i is the estimation of pi in run i-1 of direct_pi_multirun.py, while \pi = 3.1415926... is the mathematical constant.

NB: In the rms deviation, one squares the difference so that positive and negative deviations add up, rather than compensate each other. At the end, one takes a square root in order to "undo" the squaring. Take n_runs = 500, and plot the rms deviation according to the above formula as a function of n_trials for n_trials = 2^4, 2^5,...,2^12 (use log-log-scaling of axes).

For convenience (week 1 bonus), find below the code that does all this for you.

```python
import random, math, pylab

def direct_pi(N):
    n_hits = 0
    for i in range(N):
        x, y = random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0)
        if x ** 2 + y ** 2 < 1.0:
            n_hits += 1
    return n_hits

n_runs = 500
n_trials_list = []
sigmasqs = []
for poweroftwo in range(4, 13):
    n_trials = 2 ** poweroftwo
    sigmasq = 0.0
    for run in range(n_runs):
        pi_est = 4.0 * direct_pi(n_trials) / float(n_trials)
        sigmasq += (pi_est - math.pi) ** 2
    sigmasqs.append(math.sqrt(sigmasq/(n_runs)))
    n_trials_list.append(n_trials)

pylab.plot(n_trials_list, sigmasqs, 'o')
pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('number of trials')
pylab.ylabel('root mean square deviation')
pylab.title('Direct sampling of pi: root mean square deviation vs. n_trials')
pylab.savefig('direct_sampling_rms_deviation.png')
pylab.show()
```
Simply cut and paste this program into a file, then run it.

* Save the plot as a graphics file (normally, in png format). 

## A2

Modify the program from Section A1 by adding a single line in order to plot, in addition to the data of Section A1, the function 1.642/sqrt(N_trials). This line should look somewhat like:

Run the modified program.The function should nicely fit the experimental data.

* Save the new plot as a graphics file (normally, in png format).

## A2 (continued)

* Explain where in the program you added the above line.
* Explain in a few words why one can say that the "error" of the direct_pi calculation goes like 1.642 / sqrt(N_trials) (maybe explain what that means at n_trials = 100). Does the "error" of the direct-sampling algorithm go to zero as N_trials goes to infinity?

## B

We now study the Markov-chain sampling algorithm markov_pi.py. It also computes pi, but its behavior depends on the step size delta.

## B1

Copy the program of Section A1 into a new file, replace in it the function direct_pi through the function markov_pi (from markov_pi_multirun.py). Again compute the rms deviation, exactly as you did for the direct sampling algorithm. Also produce graphics output for different values of the step size delta. Choose delta = 0.062, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0. Take n_runs = 500, and plot the rms deviation according to the above formula as a function of n_trials for n_trials = 2^4, 2^5,...,2^12 (use a log-log-scaling of axes).

For convenience (week 1 bonus), find below the code that does all this for you.

```python
import random, math, pylab

def markov_pi(N, delta):
    x, y = 1.0, 1.0
    n_hits = 0
    for i in range(N):
        del_x, del_y = random.uniform(-delta, delta), random.uniform(-delta, delta)
        if abs(x + del_x) < 1.0 and abs(y + del_y) < 1.0:
            x, y = x + del_x, y + del_y
        if x**2 + y**2 < 1.0: n_hits += 1
    return n_hits

n_runs = 500
for delta in [0.062, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0]:
    n_trials_list = []
    sigmas = []
    for poweroftwo in range(4, 13):
        n_trials = 2 ** poweroftwo
        sigma = 0.0
        for run in range(n_runs):
            pi_est = 4.0 * markov_pi(n_trials, delta) / float(n_trials)
            sigma += (pi_est - math.pi) ** 2
        sigmas.append(math.sqrt(sigma/(n_runs)))
        n_trials_list.append(n_trials)
    pylab.plot(n_trials_list, sigmas, 'o', ms = 8, label = '$\delta = $' + str(delta))

pylab.xscale('log')
pylab.yscale('log')
pylab.xlabel('number of trials')
pylab.ylabel('root mean square deviation')
pylab.plot([10,10000],[1.642 / math.sqrt(10.0), 1.642 / math.sqrt(10000.0)], label = 'direct')
pylab.title('Markov-chain sampling of pi: root mean square deviation vs. n_trials')
pylab.legend(loc='upper right')
pylab.savefig('markov_sampling_rms_deviation.png')
pylab.show()
```

Simply cut and paste this program into a file, then run it.

* Upload the new plot as a graphics file (normally, in png format).

## B1 (continued):

* Which of the values of delta gives the most precise results?
* Explain why VERY small values of delta and VERY large values of delta yield a less precise result than intermediate values.
* Explain in a few words why the error is larger than for the direct sampling algorithm, even for the optimal value of delta.

## B2

The **"1/2 thumb rule"** predicts that the best performance of a Markov Chain Monte Carlo algorithm is for an acceptance ratio of approximately 1/2: half the moves are accepted, and half of them are rejected. Download the program markov_pi.py from the coursera website, and modify it so that it computes the acceptance rate. You may do a single run for n_trials = 2^12. Simply run this program for delta = 0.062, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0.

* Write down the acceptance rates you find for the different values of delta from, each time, a single run of the modified markov_pi.py. The acceptance rates should be between 0 and 1. Simply make a table as follows in your answering box.

```
delta | acceptance rate
0.062 | ??
0.125 | ??
0.25  | ??
0.5   |
1.0   |
2.0   |
4.0   |
```

* Compare with the results of Section **B1**, do you confirm that an acceptance rate of approximately 1/2 gives the best results?

## C

In section A, we saw that the error of direct_pi behaves as 1.642 / sqrt(N_trials)while, in Section B, we noticed that the error of markov_pi follows the law: const / sqrt(N_trials)for large N_trials. The constant is larger (sometimes much larger) than 1.642 and it depends on the stepsize delta.In this Section we understand what the value 1.642 means and how we can compute it with and without knowing the mathematical value of pi. We then find a good way to compute the error in markov_pi.py from a single run and without knowing the mathematical value of pi . This is the bunching method.

## C1

Modify the program direct_pi.py (from the Coursera website), so that it computes the "variance" of the relevant observable (Obs = 0.0 if the pebble is outside the circle, Obs = 4.0 if the pebble is inside the circle). To do so, simply compute the variance, namely the mean value of (Obs - math.pi) **2 and the standard deviation, namely the square root of the variance.

For convenience (week 1 bonus), this modified program is given right below. Copy and paste the following program into a file and check whether it indeed yields a value close to 1.642 for the standard deviation:

```python
import random, math
n_trials = 400000
n_hits = 0
var = 0.0
for iter in range(n_trials):
    x, y = random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0)
    Obs = 0.0
    if x**2 + y**2 < 1.0:
        n_hits += 1
        Obs = 4.0
    var += (Obs - math.pi)**2
print 4.0 * n_hits / float(n_trials), math.sqrt(var / n_trials)
```

* Write down your value for the standard deviation that you found from a single run of this modified version of direct_pi.py.

## C2

The program of Section C1 used the exact value of pi in order to compute the error. In a real calculation, we don't know this result beforehand. One therefore replaces pi through its best estimate, the mean value &lt;Obs&gt;. Using


&lt;(Obs−pi)^2&gt;∼&lt;Obs^2&gt;−pi^2∼&lt;Obs^2&gt;−&lt;Obs&gt;^2,

we see that we have to estimate &lt;Obs^2&gt; - &lt;Obs&gt;^2.

Modify the program from Section C1 so that it not only estimates the mean value of Obs (Obs can be 0 or 4), but also of Obs^2 (which can be 16 or 0). Check that the square root of &lt;Obs^2&gt; - &lt;Obs&gt;^2 is approximately equal to 1.642.

* Upload your modified program, where a few added lines compute &lt;Obs&gt; and &lt;Obs^2&gt;, &lt;Obs^2&gt; - &lt;Obs&gt;^2, and its square root.

## C2 (continued):

* Write down your value for the standard deviation, obtained from a single run of the modified direct_pi.py. It should be close to 1.642.

NB: This procedure for computing the standard deviation can be used quite generally for direct-sampling Monte Carlo algorithms. The error is quite generally given by the ratio of the standard deviation and the square root of the number of samples.

## C3

The previous Section C2 gives a general solution for the direct-sampling, and we now have to do the same for Markov chains. As we saw in Section B1, the error is again as const / sqrt(N_trials), but the constant does not equal 1.664, it may be much larger. Nevertheless, as the distribution of pebbles is the same in the direct-sampling and in Markov-chain sampling, the variance and the standard deviation are the same. The quantity sigma / sqrt(N_trials) strongly UNDERESTIMATES the error, as we saw in Section B1. This is due to the correlation between subsequent samples.To estimate the error of a Markov-chain simulation that has obtained observables Obs_0, Obs_1,...,Obs_(n_trials -1), we should rather BUNCH the data, as follows:

(Obs_0 + Obs_1)/2, (Obs_2 + Obs_3)/2, ..., (Obs_(n_trials -2) + Obs_(n_trials -1))/2

This gives a new chain of data, with new observables (mean value of two subsequent old ones), a new number of trials that is half the old one, an unchanged mean value, but a new value of the variance, and a new estimate of the error, using the naive formula for direct sampling.This famous bunching procedure can then be iterated, making pairs, then pairs of pairs, then pairs of pairs of pairs. This very successful algorithm is implemented for a long sequence of data produced by markov_pi.py ( x_i = 4 if it corresponds to a hit, and x_i = 0 otherwise) in the below program, that is again provided for convenience. Copy and past this program into a file and run it.

```python
import random, pylab, math

def markov_pi_all_data(N, delta):
    x, y = 1.0, 1.0
    data = []
    for i in range(N):
        del_x, del_y = random.uniform(-delta, delta), random.uniform(-delta, delta)
        if abs(x + del_x) < 1.0 and abs(y + del_y) < 1.0:
            x, y = x + del_x, y + del_y
        if x ** 2 + y ** 2 < 1.0:
            data.append(4.0)
        else:
            data.append(0.0)
    return data

poweroftwo = 20
n_trials = 2 ** poweroftwo
delta = 0.1
data = markov_pi_all_data(n_trials, delta)
errors  = []
bunches = []
for i in range(poweroftwo):
    new_data = []
    mean = 0.0
    mean_sq = 0.0
    N = len(data)
    while data != []:
        x = data.pop()
        y = data.pop()
        mean += x + y
        mean_sq += x ** 2 + y ** 2
        new_data.append((x + y) / 2.0 )
    errors.append(math.sqrt(mean_sq / N - (mean / N) ** 2) / math.sqrt(N))
    bunches.append(i)
    data = new_data[:]
    print mean / float(N), 'mean value, estimate of pi'
pylab.plot(bunches, errors, 'o')
pylab.xlabel('iteration')
pylab.ylabel('apparent error')
pylab.title('Bunching: naive error vs iteration number')
pylab.savefig('apparent_error_bunching.png')
pylab.show()
```

The program produces a plot of naive (apparent) error against the iteration of the bunching procedure.The observed error is found to increase with the iterations, then exhibits a plateau (to see this you may have to run the program several times). This plateau is an excellent estimation of the true error of a single Markov chain output. Run this program, then

* Upload the plot you obtained (usually in png format).

## C3 (continued):

* Explain in a few words why the apparent error initially increases with the iterations, then saturates (more or less) to a plateau.
* The program outputs the mean value (which does not change with iteration), the Monte Carlo evaluation pi^est of pi. Compare the absolute value of (pi^est - pi) with the plateau error. Are they similar?

NB: For your own interest, you can rerun this program for different values of N_trial and delta, to understand the robustness of the error evaluation. It is also very interesting to modify in the bunching algorithm the function markov_pi_all_data into a function direct_pi_all_data. You will then see a very large plateau, from the first iteration on.
