import random
 
def direct_pi(N):
    n_hits = 0
    for i in range(N):
        x, y = random.uniform(-1.0, 1.0), random.uniform(-1.0, 1.0)
        if x ** 2 + y ** 2 < 1.0:
            n_hits += 1
    return n_hits
 
n_runs = 1000
n_trials = 4000
for run in range(n_runs):
    print 4.0 * direct_pi(n_trials) / float(n_trials)
