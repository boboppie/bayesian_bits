bayesian_bits
=============

P(A|B) ∝P(B|A)P(A)

Personal collection of Bayesian statistics snippets and functions

Visualizing Bayesian Updating
----------

As in [Corey Chivers](http://madere.biol.mcgill.ca/cchivers/)'s [post](http://bayesianbiologist.com/2011/09/10/visualizing-bayesian-updating/):

One of the most straightforward examples of how we use Bayes to update our beliefs as we acquire more information can be seen with a simple *Bernoulli process*. That is, a process which has only two  possible outcomes.

Probably the most commonly thought of example is that of a coin toss. The outcome of tossing a coin can only be either heads, or tails (barring the case that the coin lands perfectly on edge), but there are many other real world examples of Bernoulli processes. In manufacturing, a widget may come off of the production line either working, or faulty.  We may wish to know the probability that a given widget will be faulty.  We can solve this using Bayesian updating.

I’ve put together this little piece of R code ([**bayesian_binomial_updating.R**](https://raw.github.com/boboppie/bayesian_bits/master/bayesian_binomial_updating.R)) to help visualize how our beliefs about the probability of success (heads, functioning widget, etc) are updated as we observe more and more outcomes.

The result is a plot of posterior (which become the new prior) distributions as we make more and more observations from a Bernoulli process.

![bayesian_binomial_updating] (https://raw.github.com/boboppie/bayesian_bits/master/plots/bayesian_binomial_updating.png)

With each new observation, the posterior distribution is updated according to Bayes rule. You can change *p* to see how belief changes for low, or high probability outcomes, and *N* for to see how belief about *p* asymptotes to the true value after many observations.

In the comment by mamluk:

It looks like it might be helpful showing people how Bayesian updating works as a process. It’s often nice to show how different priors can affect the final result, so I added some code to make the Beta parameters variables that can be specified when calling the function. ([**bayesian_binomial_updating_v2.R**](https://raw.github.com/boboppie/bayesian_bits/master/bayesian_binomial_updating_v2.R))

![bayesian_binomial_updating_v2] (https://raw.github.com/boboppie/bayesian_bits/master/plots/bayesian_binomial_updating_v2.png)

An update on visualizing Bayesian updating
----------
Corey Chivers's [following post](http://bayesianbiologist.com/2012/08/17/an-update-on-visualizing-bayesian-updating/):

A while ago I wrote a post with some R code to visualize the updating of a beta distribution as the outcome of Bernoulli trials are observed. The code provided a single plot of this process, with all the curves overlayed on top of one another. Then John Myles White (co-author of Machine Learning for Hackers) piped up on twitter and said that he’d like to see it as an animation. Challenge accepted – and with an additional twist. ([**bayesian_updating_video.R**](https://raw.github.com/boboppie/bayesian_bits/master/bayesian_updating_video.R) or [gist](https://gist.github.com/3373348))

The video shows how two observers who approach the problem with different beliefs (priors) converge toward the same conclusion about the value of the unknown parameter after making enough observations. Watch it on [Youtube](http://www.youtube.com/watch?v=rUoJvogN7qQ).

Screenshot:
![ bayesian_updating_video_screenshot](https://raw.github.com/boboppie/bayesian_bits/master/plots/bayesian_updating_video_screenshot.png)