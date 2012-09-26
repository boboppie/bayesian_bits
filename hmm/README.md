Hidden Markov Model
==========

[Python Hidden Markov Model](http://www.cs.colostate.edu/~hamiltom/code.html)
----------
[hmm.py](https://raw.github.com/boboppie/bayesian_bits/master/hmm/hmm.py) provides an implementation by [Hamiltom](http://www.cs.colostate.edu/~hamiltom/) of a hidden Markov model (HMM) for discrete observation variables. The module was written with two requirements

1. Easy-to-follow code for teaching
2. Ability to handle large datasets

For those learning about HMMs, I recommend following along with the Rabiner tutorial [1] as a reference for variable names and techniques implemented in this module. Finally, this software requires the NumPy and matplotlib modules.

Installation

There is no packaging of this software. Simply place hmm.py into your Python path.

Download

[HMMPy](http://www.cs.colostate.edu/~hamiltom/_downloads/HMMPy.tar3.gz) Contains hmm.py module and Doxygen documentation (HMMPY_doc.pdf).

References

[1]  LR Rabiner. A tutorial on hidden Markov models and selected applications in speech recognition. Proceedings of the IEEE, 1989.


[Laplacian Smoother](https://bitbucket.org/les2/aiclass/src)
----------
[Smooth.py](https://raw.github.com/boboppie/bayesian_bits/master/hmm/Smooth.py): Implementation of the Laplacian smoother discussed in the [Stanford AI class](https://www.ai-class.com/) by Lloyd Smith.

Textbook example:
----------

I do three things on a day: shopping (s), cooking(c), watching tv (w).
If it is sunny, I do those things according to probability S:C:W = 40:30:30
If it rains, the probabilities are S:C:W = 10:40:50

Independent of what I do, the weather follows a Markov chain of two-states Rain(R), Sunny(S), if it rains on a day, the next day it has 50% raining, 50% sunny if it is sunny on a day, for the next day, raining:sunny = 30:70

For M = 10000 say,

Probability questions: Q0: On average, how many days in M days will be raining? On average, how many days in M days will I be cooking? Q1: If I tell you I'm cooking now, what's the probability that it's raining outside?

Monte-Carlo: Q2: Suppose Day 0 is raining. Generate a random series of weather of the next M days (Does that match your result in Q0?) Q3: Based on the series from Q2, generate a series of my activity (Does that match your result in Q0?)

Now the HMM questions: You have no idea about the underlying weather, and I give you a length M series of my recent activity up to today Q4: What's the probability that it's raining today Q5: What's the most likely sequence of weather in those M days (eg rain,rain,sunny,sunny,sunny,...) Q6: What's the most likely activity I will do tomorrow

Solution:
----------
Q0: Stationary distribution 
* P(R_i) = P(R_i-1), p(S_i) = 1 - P(R_i)
* P(c) = P(c|R)P(R) + P(c|S)P(S)

Q1: 
* P(R|c) = P(c|R)P(R)/P(c)

Q2:

Q3:

Q4:

Q5: http://en.wikipedia.org/wiki/Viterbi_algorithm

----------
* https://www.ai-class.com/course/video/videolecture/138 (https://github.com/boboppie/stanford-ai-class/blob/master/11_HMMs_and_filters.md)
* http://mintgene.wordpress.com/2012/01/28/hmm-implementation-of-viterbi-algorithm-durbin-1998-part-1/
* http://mintgene.wordpress.com/2012/01/29/hmm-implementation-of-viterbi-algorithm-durbin-1998-part-2/
* http://mintgene.wordpress.com/2012/02/09/discrimination-between-cpg-islands-and-random-sequences-using-markov-chains/
