## A Equiprobability

In this lecture, we insisted on the equiprobability principle, which governs the statistical physics of hard disks: All legal configurations are to be equally probable.

## A1

Check that equiprobability is satisfied in direct_disks_box.py. Check it not for all legal configurations, but for the three configurations (a, b, c) shown in the picture boxes.png in this dir. The probability to hit any of these configurations exactly is of course zero, so we must put little boxes around them, as shown in red:

Using small boxes [x - del_xy, x + del_xy], etc, modify the program direct_disks_box_multirun.py to show that the probability to sample configurations a, b, and c are the same (within the numerical precision), with
```python
a = ((0.30, 0.30), (0.30, 0.70), (0.70, 0.30), (0.70,0.70))
b = ((0.20, 0.20), (0.20, 0.80), (0.75, 0.25), (0.75,0.75))
c = ((0.30, 0.20), (0.30, 0.80), (0.70, 0.20), (0.70,0.70))
```

For your convenience, the program is already provided:
```python
import random, math
def direct_disks_box(N, sigma):
    condition = False
    while condition == False:
        L = [(random.uniform(sigma, 1.0 - sigma), random.uniform(sigma, 1.0 - sigma))]
        for k in range(1, N):
            a = (random.uniform(sigma, 1.0 - sigma), random.uniform(sigma, 1.0 - sigma))
            min_dist = min(math.sqrt((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) for b in L)
            if min_dist < 2.0 * sigma:
                condition = False
                break
            else:
                L.append(a)
                condition = True
    return L


sigma = 0.15
del_xy = 0.05
n_runs = 1000000
conf_a = ((0.30, 0.30), (0.30, 0.70), (0.70, 0.30), (0.70,0.70))
conf_b = ((0.20, 0.20), (0.20, 0.80), (0.75, 0.25), (0.75,0.75))
conf_c = ((0.30, 0.20), (0.30, 0.80), (0.70, 0.20), (0.70,0.70))
configurations = [conf_a, conf_b, conf_c]
hits = {conf_a: 0, conf_b: 0, conf_c: 0}
for run in range(n_runs):
    x_vec = direct_disks_box(4, sigma)
    for conf in configurations:
        condition_hit = True
        for b in conf:
            condition_b = min(max(abs(a[0] - b[0]), abs(a[1] - b[1])) for a in x_vec) < del_xy
            condition_hit *= condition_b
        if condition_hit:
            hits[conf] += 1

for conf in configurations:
    print(conf, hits[conf])
```

Cut-and-paste this program into a file, then run it. Familiarize yourself with it (add a few print statements to see how it works). Then answer a few questions:

1. There is a logical variable "condition_hit" (it is either True or False). Under which condition is it "True"? Explain in your own words: what are the meanings of the curious "min(max(...) ..) " statement, and of "if condition_hit"?

Hint: Modify the program by adding print statements in order to find out.

2. Run the above program as follows:

    3 times for n_runs = 10^4 (write down the number of hits for configuration a, b, c each time).
    3 times for n_runs = 10^5 (write down the number of hits for configuration a, b, c each time).
    3 times for n_runs = 10^6 (same as above)

Do you get indications for equiprobability?

3. Run your program again for sigma = 0.15, and for del_xy = 0.10. Can you explain what happens?

4. Explain in one sentence, and in your own words: How does the algorithm direct_disks_box.py implement equiprobability?

NB: Note that the combination of the four little red boxes makes a single red box in 8-dimensional configuration space. We are measuring the statistical weight of three arbitrary points (corresponding to legal configurations) in this space.

## A2 (continued):

Print the program for the previous section A2.

## A3

Consider the program event_disks_box.py. Download this program from the Coursera website. Study this program. It outputs the positions of event configurations.

Prove, without any numerical simulations, that no event-configuration will count for the configurations a, b, c of Section A1. Show that this means that the event positions are not equally probable (Note that this is rather natural).

Modify the program event_disks_box.py, so that it determines the configurations at constant time intervals, t=0, Delta t, 2 Delta t, 3 Delta t,.. It is for these configurations that equal probability should be valid.

For convenience, this modified program is provided below. Study it (and compare it to the original event_disks_box.py) in order to understand the constant-time-interval mechanism.

```python
import math, pylab

def wall_time(pos_a, vel_a, sigma):
    if vel_a > 0.0:
        del_t = (1.0 - sigma - pos_a) / vel_a
    elif vel_a < 0.0:
        del_t = (pos_a - sigma) / abs(vel_a)
    else:
        del_t = float('inf')
    return del_t

def pair_time(pos_a, vel_a, pos_b, vel_b, sigma):
    del_x = [pos_b[0] - pos_a[0], pos_b[1] - pos_a[1]]
    del_x_sq = del_x[0] ** 2 + del_x[1] ** 2
    del_v = [vel_b[0] - vel_a[0], vel_b[1] - vel_a[1]]
    del_v_sq = del_v[0] ** 2 + del_v[1] ** 2
    scal = del_v[0] * del_x[0] + del_v[1] * del_x[1]
    Upsilon = scal ** 2 - del_v_sq * ( del_x_sq - 4.0 * sigma **2)
    if Upsilon > 0.0 and scal < 0.0:
        del_t = - (scal + math.sqrt(Upsilon)) / del_v_sq
    else:
        del_t = float('inf')
    return del_t

conf_a = ((0.30, 0.30), (0.30, 0.70), (0.70, 0.30), (0.70,0.70))
conf_b = ((0.20, 0.20), (0.20, 0.80), (0.75, 0.25), (0.75,0.75))
conf_c = ((0.30, 0.20), (0.30, 0.80), (0.70, 0.20), (0.70,0.70))
configurations = [conf_a, conf_b, conf_c]
hits = {conf_a: 0, conf_b: 0, conf_c: 0}
del_xy = 0.10
pos = [[0.25, 0.25], [0.75, 0.25], [0.25, 0.75], [0.75, 0.75]]
vel = [[0.21, 0.12], [0.71, 0.18], [-0.23, -0.79], [0.78, 0.1177]]
singles = [(0, 0), (0, 1), (1, 0), (1, 1), (2, 0), (2, 1), (3, 0), (3, 1)]
pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
sigma = 0.10
t = 0.0
n_events = 5000000
for event in range(n_events):
    if event % 100000 == 0:
        print(event)
    wall_times = [wall_time(pos[k][l], vel[k][l], sigma) for k, l  in singles]
    pair_times = [pair_time(pos[k], vel[k], pos[l], vel[l], sigma) for k, l in pairs]
    next_event = min(wall_times + pair_times)
    t_previous = t
    for inter_times in range(int(t + 1), int(t + next_event + 1)):
        del_t = inter_times - t_previous
        for k, l in singles:
            pos[k][l] += vel[k][l] * del_t
        t_previous = inter_times
        for conf in configurations:
            condition_hit = True
            for b in conf:
                condition_b = min(max(abs(a[0] - b[0]), abs(a[1] - b[1])) for a in pos) < del_xy
                condition_hit *= condition_b
            if condition_hit:
                hits[conf] += 1
    t += next_event
    del_t = t - t_previous
    for k, l in singles:
        pos[k][l] += vel[k][l] * del_t
    if min(wall_times) < min(pair_times):
        collision_disk, direction = singles[wall_times.index(next_event)]
        vel[collision_disk][direction] *= -1.0
    else:
        a, b = pairs[pair_times.index(next_event)]
        del_x = [pos[b][0] - pos[a][0], pos[b][1] - pos[a][1]]
        abs_x = math.sqrt(del_x[0] ** 2 + del_x[1] ** 2)
        e_perp = [c / abs_x for c in del_x]
        del_v = [vel[b][0] - vel[a][0], vel[b][1] - vel[a][1]]
        scal = del_v[0] * e_perp[0] + del_v[1] * e_perp[1]
        for k in range(2):
            vel[a][k] += e_perp[k] * scal
            vel[b][k] -= e_perp[k] * scal

for conf in configurations:
    print(conf, hits[conf])
```

Note that in this program, we use del_xy = 0.1, and sigma = 0.1, in order to improve statistics.

Then do the following:

- First answer the question about event-configurations from the beginning of the section.
- Run he above program for 5,000,000 events (as written in the code), write down how many times you were close to a, b, and c. (choose smaller number of events, if you cannot run that long). Do you get an indication in favor of equiprobability?
- What is the total running time at the end of the simulation, and at what times are the configurations analyzed (is it at time 0,1,2,3,... or at times 0,2,4,6,... or at times 0, 0.1, 0.2, 0.3....?). (Hint: put print statements to find out).

NB: It's well worth your time to run and modify this program quite a bit, just for your own instruction.


## B

In section B, we again consider four disks in the same square box of edge length 1 as before. We set the density equal to eta = 0.18, which corresponds to a disk radius sigma = 0.1197. This density is different from what it was in Section A.

Instead of the configurations themselves, and their probability distribution, we now consider an observable, in fact a particularly simple one, the position x: the x-coordinate of the center of a disk. We will compute its probability distribution, as the normed histogram of x-positions. This histogram is the same for all disks, so we can collect data for one disk or for all of them.

## B1

Compute the histogram of the x-positions for the direct-sampling Monte Carlo algorithm by modifying the algorithm direct_disks_box_multirun.py.

For convenience, simply cut and paste the below program into a file, and produce a normed histogram of the x-positions with 100 bins as output.

```python
import random, pylab

def direct_disks_box(N, sigma):
    overlap = True
    while overlap == True:
        L = [(random.uniform(sigma, 1.0 - sigma), random.uniform(sigma, 1.0 - sigma))]
        for k in range(1, N):
            a = (random.uniform(sigma, 1.0 - sigma), random.uniform(sigma, 1.0 - sigma))
            min_dist_sq = min(((a[0] - b[0]) ** 2 + (a[1] - b[1]) ** 2) for b in L)
            if min_dist_sq < 4.0 * sigma ** 2:
                overlap = True
                break
            else:
                overlap = False
                L.append(a)
    return L


N = 4
sigma = 0.1197
n_runs = 1000000
histo_data = []
for run in range(n_runs):
    pos = direct_disks_box(N, sigma)
    for k in range(N):
        histo_data.append(pos[k][0])
pylab.hist(histo_data, bins=100, normed=True)
pylab.xlabel('x')
pylab.ylabel('frequency')
pylab.title('Direct sampling: x coordinate histogram (density eta=0.18)')
pylab.grid()
pylab.savefig('direct_disks_histo.png')
pylab.show()
```

Print the histogram (check that it has an appropriate title explaining the program that was used, and also that axes are correctly labeled).

## B1 (continued)

Briefly explain (in two sentences) the paradox that gives a non-constant histogram for a manifestly equiprobable probability distribution, as we checked in section A1. Only give a hint or two... We will revisit this question during week 3 in great detail.

## B2

Use what was provided in Section B1 to produce an analogous histogram for markov_disks_box.py. Run it again for really long times, at least n_steps=2,000,000.

Print your modified program.

Attention: This part is quite subtle. You MUST record the x-positions at the right place (remember the piles in the pebble game in Lecture 1!). Compare to your result of Section B1, and don't hesitate to vary delta and to increase n_steps for high-quality comparisons with B1. Also, make sure you check consistency of the present program with the program in Section A2.

## B2 (continued)

Print the histogram. Make sure that there is an appropriate title saying which program was used. Make sure that the axes are correctly labeled
 
## B3

We now consider the event-driven Molecular dynamics algorithm. Modify your program (the modified event_disks_box.py) from Section A3 (replace the condition_hit check by the histogram update, add the figure output, etc). Make sure that you use the same radius sigma as in Sections B1 and B2.

Compute the histogram for the x-positions with 100 bins analogously to what you did in Sections B1 and B2, use as long a running time as reasonable. Make sure you compute the x-position at regular time intervals, and not at collision times. Make sure you include all the four disk positions into the histogram at all times.

Print the histogram (Provide a title and label the axes).

NB: Don't hesitate to spend extra time on the comparisons in B1, B2, B3, in order to think over what you did. Note also that, using statistical tools, we get really precise results.

## B3 (continued)

Print your modified event-driven molecular dynamics program.

## B3 (continued)

Comment on your findings: Do you find that direct-sampling, Markov-chain sampling and molecular dynamics give the same histogram for the x-position of the disks?