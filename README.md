We report our development and subsequent analysis of a Monte-Carlo Markov-Chain (MCMC) simulation for a 2-dimensional
lattice following the Ising model, capable of studying the thermodynamic properties of the system. We reproduced plots for the
magnetisation, susceptibility, magnetic susceptibility, and specific heat per spin for a system with external magnetic field, $𝐻 = 0$,
and for a system with external magnetic field, $𝐻 = 1$. For the system absent of any external fields, we observe a phase transition
between the magnetic and paramagnetic phase, with a critical temperature comparable to literature values. For the case of an
external magnetic field, the system does not appear to undergo a similar phase transition, as indicated by our results, due to the
impact of the introduced asymmetry. Overall, our simulation appears to be capable of accurately measure the thermodynamical
properties of an Ising model, and reproduce the theoretically expected figures in both presence and absence of an external
magnetic field.



\textbf{main_file.py} 	Main program containig the simulation and relevant parameters. It allows the user to measure
		          thermodynamic properties of a lattice following the Ising Model, for different numbers of
		          particles and values of the coupling constant and external magnetic field. For more information regarding
		          its proper use we refer to its docstring description.
              
              
\textbf{functions.py}	File containing all the relevant functions used by the simulation. We note that the
		          ising_mc_simu function performs the simulation on its own and can be used idependently
