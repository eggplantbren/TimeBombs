Time Bombs
==========

Inferring time delays from a single light curve, assuming a flaring mixture model. We are trying to understand the light curves of lensed (and unlensed) blazars observed with Fermi, by forward modeling the light curves in their entirety, and doing Bayesian model comparison. So far our model just has a single time delay and magnification ratio (as if there was no microlensing), but we plan to explore the more complex scenarios suggested by [Barnacka et al]() and [Neronov et al](). 

This is research in progress, with the Fermi data analysis being led by UC Davis grad student Nick Rumbaugh. If you make use of any of the ideas, algorithms, models or code in this repository, please cite it as *(Rumbaugh et al, in preparation\footnote{https://github.com/eggplantbren/TimeBombs})*.

### People, Licence, Credits etc

Currently working on this project:

* Brendon Brewer (Auckland)
* Phil Marshall (KIPAC)
* Nick Rumbaugh (UC Davis)
* Jim Chiang (SLAC)

Feel free to send us email, tweets, issues, pull requests etc. The code is still under development and not yet written up, so beware. It'd be lovely to hear from you though! 

(c) 2014,2015 Brendon J. Brewer, Philip J. Marshall, Nicholas A. Rumbaugh

LICENCE: GNU GPL v3.

### Getting Started

Dependencies include:

* DNest3 (http://github.com/eggplantbren/DNest3)
* RJObject (http://github.com/eggplantbren/RJObject)

Watch this space for installation notes.

Compile with:

    make

Then run with:

    ./main -t 8
The option causes the program to run with 8 threads.

To postprocess results and create the file posterior_sample.txt:

    python2 showresults.py



