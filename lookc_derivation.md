
#### Appendix: leave-one-out knockoff construction (LOOKC)

Plain vanilla (a.k.a. "second-order") knockoffs are constructed such that the centered, scaled data $X$ and the knockoffs $X$ have joint distribution 
$$G = \begin{bmatrix}\Sigma & \Sigma - S \\ \Sigma -S & \Sigma \end{bmatrix}$$. 
$S$ is a diagonal matrix that can be specified by the user. There are some constraints on $S$, and existing software can determine a good option based on $\Sigma$. Since $X$ is known but $\tilde X $ must be generated, the sample is drawn from $Pr(\tilde X| X)$, which is Gaussian with mean $X - X\Sigma^{-1} S$ and covariance $C = 2S - S \Sigma^{-1} S$. This involves matrix operations of order $O(ND^2)$ and $O(D^3)$ where $D$ is the number of variables and $N$ the number of observations.

The challenge of leave-one-out construction is to avoid increasing the computation time by a factor of $D$. This is possible via some rank-one updates as follows. Let $S_k$, $\Sigma_k$, and $X_k$ denote the obvious matrices but with variable $k$ omitted. Likewise let $G_k$ denote $G_k$ omitting row and column $k$ and $k+D$, and let $C_k = 2S_k - S_k \Sigma_k^{-1} S_k$. How do we determine these matrices cheaply?

- $S_k$: Mathematically, the main requirement is that $S_k$ be chosen such that $G_k$ remain positive definite, but this is already satisfied because for a positive definite matrix, and principal submatrix is also positive definite. For multipication by $\Sigma_k^{-1}$, a fast solution is possible via a rank-1 update to a pre-computed Cholesky decomposition of $\Sigma$. 
- $C_k$: working on it


#### Cholesky update to get $L_kL_k^T = \Sigma_k$ from $LL^T = \Sigma$

hmmmmm

#### References

Krause, O., & Igel, C. (2015, January). A more efficient rank-one covariance matrix update for evolution strategies. In Proceedings of the 2015 ACM Conference on Foundations of Genetic Algorithms XIII (pp. 129-136).