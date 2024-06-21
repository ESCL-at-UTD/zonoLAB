# **Welcome to zonoLAB!** 

# Introduction

The zonoLAB toolbox is a MATLAB-based software package for generating, manipulating, and analyzing _hybrid zonotopes_, a set representation that is equivalent to the union of $2^N$ constrained zonotopes (i.e. convex polytopes). Hybrid zonotopes are represented in _Hybrid Constrained Generator-representation_ (HCG-rep), with the shorthand notation $\mathcal{Z}_h = \langle G^c,G^b,c,A^c,A^b,b\rangle\subset\mathbb{R}^n$, such that 

$$\mathcal{Z}_h = \lbrace z = \left[G^c \ G^b\right] \begin{bmatrix}\xi^c \atop\xi^b \end{bmatrix}  + c \ | \ \begin{bmatrix}\xi^c \atop\xi^b \end{bmatrix} \in \mathcal{B}\_{\infty}^{n\_{g}} \times \lbrace -1,1\rbrace^{n\_{b}}, \ \left[A^c \ A^b\right] \begin{bmatrix}\xi^c \atop\xi^b \end{bmatrix}  = b \rbrace, $$

where:   
* $n$ denotes the dimension of the space,
* $n\_g$ denotes the number of continuous generators forming the columns of $G^c \in \mathbb{R}^{n \times n\_g}$, which are multiplied by the $n\_g$ factors $\xi^c \in \mathbb{R}^{n\_g}$,
* $n\_b$ denotes the number of binary generators forming the columns of $G^b \in \mathbb{R}^{n \times n\_b}$, which are multiplied by the $n\_b$ factors $\xi^b \in \mathbb{R}^{n\_b}$,
* $n\_c$ denotes the number of linear equality constraints,
* $c\in\mathbb{R}^n$ is the center,
* $\mathcal{B}\_{\infty}^{n\_{g}}$ is the $n\_g$-dimensional unit hypercube such that $\mathcal{B}\_{\infty}^{n\_{g}} = \lbrace x \in \mathbb{R}^{n} \ | \ \|| x \||\_\infty \leq 1 \rbrace$,
* $\lbrace -1,1\rbrace^{n\_{b}}$ is the set of all $n\_{b}$-dimensional binary vectors (i.e. vertices of $\mathcal{B}\_{\infty}^{n\_{b}}$),
* $A\_c \in \mathbb{R}^{n\_c \times n\_g}$ is the continuous factor constraint coefficient matrix,
* $A\_b \in \mathbb{R}^{n\_c \times n\_b}$ is the binary factor constraint coefficient  matrix, and 
* $b \in \mathbb{R}^{n\_c}$ is the constraint offset vector.

In this toolbox, hybrid zonotopes are represented using the class [hybZono](https://github.com/ESCL-at-UTD/zonoLab/tree/main/%40hybZono). Constrained zonotopes and zonotopes can be viewed as specific cases of hybrid zonotopes and are also represented using the classes [conZono](https://github.com/ESCL-at-UTD/zonoLab/tree/main/%40conZono) and [zono](https://github.com/ESCL-at-UTD/zonoLab/tree/main/%40zono), respectively. All three classes inherit properties from the abstract superclass [abstractZono](https://github.com/ESCL-at-UTD/zonoLab/tree/main/%40abstractZono), which contains the properties and methods commonly used by [hybZono](https://github.com/ESCL-at-UTD/zonoLab/tree/main/%40hybZono), [conZono](https://github.com/ESCL-at-UTD/zonoLab/tree/main/%40conZono), and [zono](https://github.com/ESCL-at-UTD/zonoLab/tree/main/%40zono).


# Installation
The zonoLab toolbox does not require installation beyond [adding the folder and subfolders to the path in MATLAB](https://www.mathworks.com/help/matlab/matlab_env/add-remove-or-reorder-folders-on-the-search-path.html). Running [zonoLab_Install](https://github.com/ESCL-at-UTD/zonoLAB/blob/main/zonoLab_Install.m) will automatically add all necessary folders to your MATLAB path. Additionally, it is highly recommended that Gurobi is installed and added to your MATLAB path. Information on installing Gurobi, including details on free license for academic use, can be found on the [Getting Started with Gurobi Optimizer page](https://support.gurobi.com/hc/en-us/articles/14799677517585-Getting-Started-with-Gurobi-Optimizer).

# Getting Started

The following code can be used to create and plot your first hybrid zonotope. 

```matlab
Gc = [3 -3 1; 2 1 -2];          % Continuous generator matrix
Gb = [6 -6 2; 4 2 -4];          % Binary generator matrix
c = [0;0];                      % Center
Ac = [1 1 1];                   % Continuous constraint matrix
Ab = [1 1 1];                   % Binary constraint matrix
b = 1;                          % Constraint offset vector
Zh = hybZono(Gc,Gb,c,Ac,Ab,b);  % Creates a hybrid zonotope

figure;
plot(Zh,'b',0.1)                % Plots hybrid zonotope in transparent blue
```

<img src="https://github.com/ESCL-at-UTD/zonoLAB/blob/main/dev/figures/GetStartedHybZono.png" width="400">


To learn more about generating hybrid zonotopes, plotting, and set-based operations, please visit the Wiki [Tutorial](https://github.com/ESCL-at-UTD/zonoLAB/wiki/Tutorial), [Functions](https://github.com/ESCL-at-UTD/zonoLAB/wiki/Functions), and [Examples](https://github.com/ESCL-at-UTD/zonoLAB/wiki/Examples) pages.

# License

The zonoLab toolbox is distributed under the [GNU General Public License (GPL)](https://github.com/ESCL-at-UTD/zonoLab/blob/main/LICENSE). The paper provided above must be referenced when zonoLab is used in published work. This toolbox is distributed without any warranty and, therefore, the user is responsible for assessing the correctness of the software.

# Cite

Please cite the following publication if you publish work based on the zonoLAB toolbox:

* Justin P. Koeln, Trevor J. Bird, Jacob Siefert, Justin Ruths, Herschel C. Pangborn, Neera Jain, [zonoLAB: A MATLAB toolbox for set-based control systems analysis using hybrid zonotopes](https://arxiv.org/abs/2310.15426), arXiv.2310.15426v2, 2024.

For general background on hybrid zonotopes, refer to the following publication:

* Trevor J. Bird, Herschel C. Pangborn, Neera Jain, Justin P. Koeln, [Hybrid zonotopes: A new set representation for reachability analysis of mixed logical dynamical systems](https://www.sciencedirect.com/science/article/pii/S0005109823002674), Automatica, Volume 154, 2023. ([preprint with full proofs](https://arxiv.org/pdf/2106.14831.pdf))

# Additional References

More information on the hybrid zonotopes and applications can be found in the following publications:

* Jacob A. Siefert, Trevor J. Bird, Andrew F. Thompson, Jonah J. Glunt, Justin P. Koeln, Neera Jain, Herschel C. Pangborn, [Reachability Analysis Using Hybrid Zonotopes and Functional Decomposition](https://arxiv.org/abs/2304.06827), arXiv:2304.06827v2, 2024.

* Joshua Ortiz, Alyssa Velluci, Justin Koeln, Justin Ruths, [Hybrid Zonotopes Exactly Represent ReLU Neural Networks](https://ieeexplore.ieee.org/abstract/document/10383944), IEEE Conference on Decision and Control, 2023.

* Jacob A. Siefert, Andrew F. Thompson, Jonah J. Glunt, Herschel C. Pangborn, [Set-Valued State Estimation for Nonlinear Systems Using Hybrid Zonotopes](https://ieeexplore.ieee.org/abstract/document/10383789), IEEE Conference on Decision and Control, 2023.

* Jacob A. Siefert, Trevor J. Bird, Justin P. Koeln, Neera Jain, Herschel C. Pangborn, [Successor Sets of Discrete-time Nonlinear Systems Using Hybrid Zonotopes](https://ieeexplore.ieee.org/abstract/document/10156300), American Control Conference, 2023.

* Jacob A. Siefert, Trevor J. Bird, Justin P. Koeln, Neera Jain, Herschel C. Pangborn, [Robust Successor and Precursor Sets of Hybrid Systems using Hybrid Zonotopes](https://ieeexplore.ieee.org/abstract/document/9815835), IEEE Control Systems Letters, 2022.

* Trevor J. Bird, Neera Jain, Herschel C. Pangborn, Justin P. Koeln, [Set-Based Reachability and the Explicit Solution of Linear MPC using Hybrid Zonotopes](https://ieeexplore.ieee.org/abstract/document/9867853), American Control Conference, 2022.

* Trevor J. Bird, Neera Jain, [Unions and Complements of Hybrid Zonotopes](https://ieeexplore.ieee.org/abstract/document/9638970), IEEE Control Systems Letters, 2021.

Dissertations:

* Jacob A. Siefert, [Reachability Analysis of Nonlinear and Hybrid Systems Using Hybrid Zonotopes and Graphs of Functions](https://etda.libraries.psu.edu/catalog/30361jas7031), Mechanical Engineering (PhD), Penn State University, 2023.

* Trevor J. Bird, [Hybrid Zonotopes: A Mixed-Integer Set Representation for the Analysis of Hybrid Systems](https://hammer.purdue.edu/articles/thesis/Hybrid_Zonotopes_A_Mixed-Integer_Set_Representation_for_the_Analysis_of_Hybrid_Systems/21225332), Mechanical Engineering (PhD), Purdue University, 2022.

