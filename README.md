bayesian_bits
=============

P(A|B) ∝P(B|A)P(A)

Personal collection of Bayesian statistics snippets and functions

Visualizing Bayesian Updating
----------

As in [Corey Chivers](http://madere.biol.mcgill.ca/cchivers/)'s [post](http://bayesianbiologist.com/2011/09/10/visualizing-bayesian-updating/):

One of the most straightforward examples of how we use Bayes to update our beliefs as we acquire more information can be seen with a simple *Bernoulli process*. That is, a process which has only two  possible outcomes.

Probably the most commonly thought of example is that of a coin toss. The outcome of tossing a coin can only be either heads, or tails (barring the case that the coin lands perfectly on edge), but there are many other real world examples of Bernoulli processes. In manufacturing, a widget may come off of the production line either working, or faulty.  We may wish to know the probability that a given widget will be faulty.  We can solve this using Bayesian updating.

I’ve put together this little piece of R code to help visualize how our beliefs about the probability of success (heads, functioning widget, etc) are updated as we observe more and more outcomes.

The result is a plot of posterior (which become the new prior) distributions as we make more and more observations from a Bernoulli process.

![bayesian_binomial_updating] (https://raw.github.com/boboppie/bayesian_bits/master/plots/bayesian_binomial_updating.png)

With each new observation, the posterior distribution is updated according to Bayes rule. You can change *p* to see how belief changes for low, or high probability outcomes, and *N* for to see how belief about *p* asymptotes to the true value after many observations.