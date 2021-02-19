# Maplestory Star Force Calculator

## Game context

This application calculates the complete distribution of costs expected when performing
[star force upgrades](https://strategywiki.org/wiki/MapleStory/Spell_Trace_and_Star_Force#Star_Force_Enhancement),
which enable you to essentially pay in-game gold (mesos) in exchange for a *chance* at upgrading equipment stats.
The linked article includes the exact calculations for the costs and probabilities of success for a *single* attempt.
However, the actual cost required to achieve a specific goal (for example, starting at 10 stars and going to 17 stars)
is much less clear. In fact, to my knowledge, there is no closed form solution to exactly calculate the complete
distribution. (It looks like [there used to be](https://www.ocf.berkeley.edu/~ted/maplestory/starforce/methodology/),
but this information is now outdated; the costs now vary depending on the number of stars, and you can realistically
go up to 22 stars.)

## Existing work

There are existing tools like [this one](https://blushiemagic.github.io/Maplestory-Starforce-Calculator/) and
[this one](https://brendonmay.github.io/starforceCalculator/) which make it easier to estimate the cost of star forcing.
However, these calculate only the expected (*mean*) cost, which is sensitive to outliers and provides only a vague picture
of how much one can expect to spend, or simulate a given number of random trials, which is sensitive to actual variance.

The former approach is extremely fast and easy to use, since calculating the mean directly takes a small, constant number of
arithmetic operations. The latter approach is reasonably fast, straightforward to understand, and allows an arbitrary
trade-off of performance for precision. Given enough random trials, the observed data is guaranteed to converge.

## Approach

I was personally interested in finding an approach that didn't have such trade-offs, i.e. (1) calculated a complete
distribution, and (2) was fast enough to effectively generate a complete table mapping `(start, target)` pairs (and some
additional in-game factors) to those distributions in just a few seconds. Of course, it is possible to simply run
offline simulations overnight, and save/publish that data. But I was interested more in the general problem than the results -
the existing methods are already quite sufficient for practical usage.

This project takes the basic approach of computing sums of random variables. For example, the 10->12 cost distribution is
simply the sum of the distributions for 10->11 and 11->12. This is a little more complex that it seems. (We denote the cost
for a given (start, target) pair as `C(i,j)`, e.g. `C(10,11)`.)

First, there is a computational dependence: computing `C(11,12)` requires computing `C(10,11)`, because one of the outcomes of
an attempt at 11->12 is *losing* a star. In fact this is not just a sum of two distributions, but a compound/mixture distribution
where the components are partial sums corresponding to the number of attempts required. For example,

Second, c
