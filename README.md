# epma_thins
Electron Probe MicroAnalysis - thin film analysis, based on [gmrfilm](https://github.com/openafox/gmrfilm) (some 'port' some not)

A more OOP style code base has been created and calculation algorithms are being traced and checked against the literature (citations are being added as comments).

The hope is that once complete other models (algorithms) can be easily implemented and a GUI wrapper could also be developed.

Currently the base OOP builds an analysis_sample with film_layers (stored from substrate to top) including one to many atomic_elements.

Much of the code is working to test the correctness of the algorithms but the main loop and minimization has not yet been fully implemented.

Some math checks and more info are in the [gh-pages](http://openafox.com/epma_thins) site for this project.

Currently bethe.py will reproduce plots seen in ref [1] but that is where things stop.

Any help making this a functional tool for epma thinfilm analysis would be greatly appreciated. Please see the contributing @ [epma_thins project site](http://openafox.com/epma_thins/contributing).
