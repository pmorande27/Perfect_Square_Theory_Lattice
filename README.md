# Perfect Square Theory in the Lattice

Repository used to investigate a scalar 4D perfect square field theory in the Lattice.
The Euclidean action is given by:
$$S_E = \frac{1}{2}\int dx^4 \left(\Delta \phi(x) + \lambda(\partial \phi(x))^2\right)^2$$

For the simulation, we rescale the field by $\chi = \lambda \phi$, which gives:
$$S_E = \frac{1}{2\lambda^2}\int dx^4 \left(\Delta \chi(x) + (\partial \chi(x))^2\right)^2$$
And we introduce the inverse coupling $\beta = \frac{1}{\lambda^2}$, such that $\beta$ acts as an inverse temperature:
$$S_E = \frac{\beta}{2}\int dx^4 \left(\Delta \chi(x) + (\partial \chi(x))^2\right)^2$$
Given that this is an asymptotically free theory, we expect that the continuum limit is achieved as $\beta \to \infty$.
