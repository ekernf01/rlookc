
## Appendix: leave-one-out knockoff construction (LOOKC)

In order to generate knockoffs, what quantities must be considered mathematically or computed explicitly? Gaussian knockoffs are constructed such that the centered, scaled data $X$ and the knockoffs $X$ have joint covariance 
$$G = \begin{bmatrix}\Sigma & \Sigma - S \\ \Sigma -S & \Sigma \end{bmatrix}$$. 
$S$ is a diagonal matrix that can be specified by the user. There are some constraints on $S$, and existing software can determine a good option based on $\Sigma$. Since $X$ is known but $\tilde X $ must be generated, the sample is drawn from $Pr(\tilde X| X)$. This distribution can be derived with standard techniques and is given in the model-X knockoffs paper. The exact formulas are given in the bulleted list below. 

Computing this involves matrix operations of order $O(ND^2)$ and $O(D^3)$ where $D$ is the number of variables and $N$ the number of observations. The challenge of leave-one-out construction is to avoid increasing the computation time by a factor of $D$. This turns out to be possible, but first, some preliminaries:

- Let $M$ and $C$ be the desired mean and covariance of the knockoffs with no variables omitted. That means 
    - $M = X - X\Sigma^{-1}S$
    - $C = 2S - S\Sigma^{-1} S$. 
- Let $L$ be such that $R^TR=C$. Any square root will suffice but the code uses a Cholesky factor.
- Without loss of generality, assume we wish to omit the final column of $X$, and call this variable $k$. 
- Let $S_{-k}$, $\Sigma_{-k}$, $X_{-k}$, $M_{-k}$, $C_{-k}$, and $R_{-k}$ denote the obvious matrices but with variable $k$ omitted. For $S$, $\Sigma$, and $C$, omitting a variable means omitting the column and the row. For $X$, $M$, and $R$, only the column is omitted. An important consequence is that $C_{-k} = R_{-k}^TR_{-k}$.
- Since $G$ is twice as tall and wide, let $G_{-k}$ denote $G$ but omitting variables $k$ and $k+D$. Both rows and columns are omitted. 
- Let $\tilde M$ and $\tilde C$ be the desired mean and covariance of the $k$th knockoff. That means 
    - $\tilde M = X_{-k} - X_{-k}(\Sigma_{-k})^{-1}S_{-k}$
    - $\tilde C = 2S_{-k} - S_{-k} (\Sigma_{-k})^{-1} S_{-k}$. 


To obtain valid knockoffs, one requirement is that $S_{-k}$ be chosen such that $G_{-k}$ remains positive definite. This is already satisfied because for a positive definite matrix, any principal submatrix is also positive definite. This is an example of a useful initial guess towards an efficient algorithm: pre-compute knockoffs with all variables and omit portions of them at each iteration. Beyond $S$, this method is not quite correct on its own, which is why we have defined $M_{-k}$ separately from $\tilde M$ and $C_{-k}$ separately from $\tilde C$, which we did not bother to do with $S$. These initial guesses $M_{-k}$ and $C_{-k}$ can be corrected quite efficiently, which we will now show by comparing $M_{-k}$ to $\tilde M$ and $C_{-k}$ to $\tilde C$. 

Before that comparison, there is one preliminary to discuss. Partition
 
$$
\Sigma = 
\begin{bmatrix} 
A & c^T \\
c & d
\end{bmatrix}
$$
and
$$
\Sigma^{-1} = 
\begin{bmatrix} 
E & g^T \\
g & h
\end{bmatrix} 
$$. In general, $A^{-1} \neq E$, but this can be resolved with a standard rank-one update: 
$$A^{-1} = E - gg^T/h$$
. This is useful in updating both the mean and the covariance. 

#### Mean

For the mean, notice 
$$M = [M_{-k}| M_k] = [X_{-k}|x_k] - [X_{-k}|x_k] \begin{bmatrix} 
E & g^T \\
g & h
\end{bmatrix} S$$, 
so
$$M_{-k} = X_{-k} - (X_{-k}E + x_kg)S_{-k}$$. By contrast, we need
$$\tilde M = X_{-k} - X_{-k}(\Sigma_{-k})^{-1}S_{-k} = X_{-k} - X_{-k}A^{-1}S_{-k}$$. The necessary correction is only of rank two:

$$
\begin{align}
\tilde M - M_{-k} 
&= -X_{-k}A^{-1}S_{-k} + (X_{-k}E + x_kg)S_{-k} \\
&= -X_{-k}ES_{-k} + X_{-k}gh^{-1}g^TS_{-k} + X_{-k}ES_{-k} + x_kgS_{-k} \\
&= X_{-k}gh^{-1}g^TS_{-k} + x_k g S_{-k}
\end{align}
$$.

#### Covariance

For the covariance, 

$$
\begin{align}
\tilde C 
&= 2S_{-k} - S_{-k} (\Sigma_{-k})^{-1} S_{-k} \\
&= 2S_{-k} - S_{-k} A^{-1} S_{-k} \\
&= 2S_{-k} - S_{-k} E S_{-k} + S_{-k} gh^{-1}g S_{-k} \\
&= C_{-k} + S_{-k} gh^{-1}g S_{-k} \\
\end{align}
$$ .

Thus, the covariance of the precomputed knockoffs can be fixed by adding an independently drawn random vector $S_{-k} gh^{-1/2}z_n$ where $z_n \sim N(0,1)$. This must be done once per row of $X_{-k}$.