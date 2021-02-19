# Maplestory Star Force Calculator

## Game context

This application calculates the complete distribution of costs expected when
performing [star force
upgrades](https://strategywiki.org/wiki/MapleStory/Spell_Trace_and_Star_Force#Star_Force_Enhancement),
which enable you to essentially pay in-game gold (mesos) in exchange for a
*chance* at upgrading equipment stats. The linked article includes the exact
calculations for the costs and probabilities of success for a *single* attempt.
However, the actual cost required to achieve a specific goal (for example,
starting at 10 stars and going to 17 stars) is much less clear. In fact, to my
knowledge, there is no closed form solution to exactly calculate the complete
distribution. (It looks like [there used to
be](https://www.ocf.berkeley.edu/~ted/maplestory/starforce/methodology/), but
this information is now outdated; the costs now vary depending on the number of
stars, and you can realistically go up to 22 stars.)

## Existing work

There are existing tools like [this
one](https://blushiemagic.github.io/Maplestory-Starforce-Calculator/) and [this
one](https://brendonmay.github.io/starforceCalculator/) which make it easier to
estimate the cost of star forcing. However, these calculate only the expected
(*mean*) cost, which is sensitive to outliers and provides only a vague picture
of how much one can expect to spend, or simulate a given number of random
trials, which is sensitive to actual variance.

The former approach is extremely fast and easy to use, since calculating the
mean directly takes a small, constant number of arithmetic operations. The
latter approach is reasonably fast, straightforward to understand, and allows an
arbitrary trade-off of performance for precision. Given enough random trials,
the observed data is guaranteed to converge.

## Approach

I was personally interested in finding an approach that didn't have such
trade-offs, i.e. (1) calculated a complete distribution, and (2) was fast enough
to effectively generate a complete table mapping `(start, target)` pairs (and
some additional in-game factors) to those distributions in just a few seconds.
Of course, it is possible to simply run offline simulations overnight, and
save/publish that data. But I was interested more in the general problem than
the results - the existing methods are already quite sufficient for practical
usage.

This project takes the basic approach of computing sums of random variables. For
example, the 10->12 cost distribution is simply the sum of the distributions for
10->11 and 11->12. This is a little more complex that it seems. (Denote the cost
for a given (start, target) pair as `C(i,j)`, e.g. `C(10,11)`.)

First, there is a computational dependence: computing `C(11,12)` requires
computing `C(10,11)`, because one of the outcomes of an attempt at 11->12 is
*losing* a star. In fact this is not just a sum of two distributions, but a
compound/mixture distribution where the components are partial sums
corresponding to the number of attempts required.

Second, computing the sums of two distributions is not actually computationally
easy. The sum of two *independent* random variables is not just a scalar
product, i.e. `C(10,11) + C(10,11) != 2*C(10,11)`, but the convolution. And a
convolution, since it scales quadratically, is actually somewhat computationally
expensive once the distributions get large enough, which they do. The number of
possible costs of a sequence of attempts grows exponentially, as there is no
nice relation between the costs at different stars.

### Algorithm

I won't describe the complete details in this README, but there are two main
tricks:

To solve the issue of dependence, I compute `C` in order of increasing targets,
and decreasing starts. That is, the following sequence: `C(10,11), C(11,12),
C(10,12), C(12,13), C(11,13), C(10,13),...`. This also solves the issues of
"booms", i.e. the outcome of destroying an item which also resets it to 12
stars.

To make computing sums of distributions tractable, I apply exponential binning.
That is, we can approximate the actual cost of any particular event, with a
fixed number of costs which increase at a slow but exponential rate,
guaranteeing an upper bound on the relative error. In this case I chose a rate
of 1%, i.e. the bins grow as `1.01^n`, which also immediately guarantees the
rounding error from binning does not exceed 1% (in fact, ~0.5%).

### Results

I can compute all the distributions of costs for a single configuration and all
(start, target) pairs with `start >= 10` and `target <= 22` in ~0.5s.

### Limitations

This approach is still numerical; it approximates the true distribution,
rounding the actual costs of a given outcome, and truncating a distribution at
the last one-millionth of its tail. As such it gives results for the mean which
can be 1-2% off. Since the highest probability events never suffer from
truncation, only the rounding error, that difference is mostly likely a result
of artifacts in the tail such as missing outliers.
