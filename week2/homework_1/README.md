A

We first study the direct-sampling algorithm direct_pi.py and modify it somewhat.

A1

Consider the program "direct_pi_multirun.py". Modify it so that it computes the root mean square (rms) deviation:

rms deviation :math:`=\sqrt(1/n_runs \sum_{i=1}^{n_runs} (\pi^est_i−\pi)^2)`

, where \pi^est_i is the estimation of pi in run i of direct_pi_multirun.py, while pi = 3.1415926... is the mathematical constant.

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

    Save the plot as a graphics file (normally, in png format). 

A2

Modify the program from Section A1 by adding a single line in order to plot, in addition to the data of Section A1, the function 1.642/sqrt(N_trials). This line should look somewhat like:

Run the modified program.The function should nicely fit the experimental data.

    Save the new plot as a graphics file (normally, in png format).

A2 (continued)

    Explain where in the program you added the above line.
    Explain in a few words why one can say that the "error" of the direct_pi calculation goes like 1.642 / sqrt(N_trials) (maybe explain what that means at n_trials = 100). Does the "error" of the direct-sampling algorithm go to zero as N_trials goes to infinity?