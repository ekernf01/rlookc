
## Appendix: leave-one-out knockoff construction (LOOKC)

Gaussian knockoffs are constructed such that the centered, scaled data $X$ and the knockoffs $\tilde X$ have joint covariance 
$$G = \begin{bmatrix}\Sigma & \Sigma - S \\ \Sigma -S & \Sigma \end{bmatrix}$$. 
Since the mean is 0 and the distribution is Gaussian, this covariance matrix completely specifies the distribution. Here, $\Sigma$ is the covariance of $X$ or an estimate thereof. $S$ is a diagonal matrix that can be specified by the user. There are some constraints on $S$, and the choice can affect the method's power. Existing software can determine a good option for $S$ based on $\Sigma$. Since $X$ is known but $\tilde X $ must be generated, the sample is drawn not from $Pr(\tilde X, X)$ but from $Pr(\tilde X| X)$. This distribution can be derived with standard techniques and is given in the model-X knockoffs paper. The exact formulas in terms of $X$, $S$, and $\Sigma$ are reproduced below as needed. 

Generating knockoffs involves matrix operations of order $O(ND^2)$ and $O(D^3)$ where $D$ is the number of variables and $N$ the number of observations. In general, knockoffs depend on $X$ but not on $Y$, so whenever a new variable is treated as $Y$, the construction would need to be repeated with that variable left out. If done naively, this would add a painful factor of $D$ to the runtime (where $D$ is the number of variables). Fortunately, it is possible to avoid this and generate all leave-one-out knockoffs (LOOKs) within a constant factor of the original $O(ND^2 + D^3)$ computational requirement. The details are explained below, but first, some preliminaries:

- Let $M$ and $C$ be the desired mean and covariance of knockoffs with no variables omitted. That means 
    - $M = X - X\Sigma^{-1}S$
    - $C = 2S - S\Sigma^{-1} S$. 
- Without loss of generality, assume we wish to omit the final column of $X$, and call this variable $k$. 
- Let $S_{-k}$, $\Sigma_{-k}$, $X_{-k}$, $\tilde X_{-k}$, $M_{-k}$, and $C_{-k}$ denote the obvious matrices but with variable $k$ omitted. For $S$, $\Sigma$, and $C$, omitting a variable means omitting the column and the row. For $X$ and $M$, only the column is omitted. This implies
    $$E[\tilde X_{-k} | X] = M_{-k}$$
    and
    $$Cov[\tilde X_{-k} | X] = C_{-k}$$.   
- Let $G_{-k}$ denote $G$ but omitting variables $k$ and $k+D$. Both rows and columns are omitted. $G_{-k}$ is never computed explicitly, but it is important mathematically because it specifies the joint distribution $Pr(\tilde X_{-k}, X_{-k})$. To obtain valid knockoffs, one requirement is that $G_{-k}$ must remain positive definite. This is satisfied because for a positive definite matrix, any principal submatrix is also positive definite. Thus, no changes are needed to the entries of $S$ or $\Sigma$.
- Given that $G_{-k}$ specifies the correct joint distribution $Pr(X_{-k}, \tilde X_{-k})$, let $\tilde M$ and $\tilde C$ be the desired mean and covariance of the distribution we actually need to sample from: $Pr(\tilde X_{-k}| X_{-k})$. That means 
    - $\tilde M = E(X_{-k} | X_{-k}) = X_{-k} - X_{-k}(\Sigma_{-k})^{-1}S_{-k}$
    - $\tilde C = Cov(X_{-k} | X_{-k}) = 2S_{-k} - S_{-k} (\Sigma_{-k})^{-1} S_{-k}$. 


A reasonable guess would be to pre-compute knockoffs with all variables and omit portions of them at each iteration. This method is not quite correct on its own since $Pr(\tilde X_{-k}| X_{-k}) \neq Pr(\tilde X_{-k}| X)$. This is why we have defined $M_{-k} = E(\tilde X_{-k}| X_{-k})$ separately from $\tilde M = E(\tilde X_{-k}| X_{-k})$ and similar for $\tilde C$ and $C_{-k}$. These initial guesses, while not exact, are very close. They can be corrected efficiently, which we will now show by comparing $M_{-k}$ to $\tilde M$ and $C_{-k}$ to $\tilde C$. 

Before that comparison, there is one preliminary to discuss. Partition $\Sigma$ and $\Sigma^{-1}$ as
 
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
. This is useful in fixing both the mean and the covariance. 

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
$$
. 

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

Thus, the covariance of the precomputed knockoffs can be fixed by adding an independently drawn random vector $S_{-k} gh^{-1/2}z_n$ where $z_n \sim N(0,1)$. This must be $N$ times, once per observation in $X$.