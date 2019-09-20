# SpatioTemporalSystems.jl
Simulations of spatio temporal dynamical systems

## API
For consistency and readability, this is the API that we're aiming for:

1. Spatio temporal systems are each represented with a function (in CamelCase) representing the name of the system.
1. All functions return `t, u, prob`, time vector, solution vector and the `DEProb` instance (for further integration with **DynamicalSystems.jl**).
2. Arguments to the functions are parameters of the system. Keyword arguments are exclusively used for solving the system (through DiffEq).
3. Initial condition is given explicitly to all functions as the first argument. For systems that have multiple couple fields, the initial condition is simpy a Tuple of all the fields.
4. Details about the parameters, algorithms and systems are given in the documentation strings of the functions.

### To-dos
* Is it possible to also allow some kind of "seed" as the initial condition? Quite tricky to do, as information about the dimensionality of space, as well as how many coupled fields are there, needs to be provided.
