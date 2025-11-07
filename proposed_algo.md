
## Problem Definition

We are given:

-   A  **directed graph**  $G = (V, E)$ with $|V| = N$.
    
-   Each node $v \in V$ has a  **color (group)**  $c(v) \in {1, 2, \dots, m}$, with $m \le N$.
    
-   For each node $v$, define its  **color in-degree vector**:
    

$$
deg_{in}(v) = [d_1(v), d_2(v), \dots, d_m(v)]  
$$

where $d_k(v)$ = number of incoming edges to $v$  **from nodes of color $k$**.

The  **coloring is valid**  if for any two nodes $u, v$ with the same color $c(u) = c(v)$, we have:

$$deg_{in}(u) = deg_{in}(v) $$

##  Goal

Add a minimal set of edges to make this property true — i.e., to make every group have identical color-in-degree vectors across its members.  
(Since the source nodes are only constrained by color, not by identity, we can choose  _any_  nodes of a given color to provide those missing edges.)

----------

## Algorithm

### Step 1. Compute color–color in-degree matrix per node

For each node $v$:

1.  Initialize an array $deg_{in}[v][1..m] = 0$.
    
2.  For each incoming edge $(u \to v)$:
    
    -   Increment $deg_{in}[v][c(u)]$ by one.
        

This gives every node’s color-in-degree vector.

**Time complexity:**  $O(|E|)$.

----------

### Step 2. Compute target in-degree vector per color group

For each color $g = 1 \dots m$:

1.  Consider all nodes $v$ with $c(v) = g$.
    
2.  For each color $k = 1 \dots m$:
    
    -   Compute $t_{g,k} = \max_{v: c(v)=g} deg_{in}[v][k]$.
        

This defines the  **target color-in-degree vector**  for group $g$:  
$$ 
T_g = [t_{g,1}, t_{g,2}, \dots, t_{g,m}]  
$$

Justification:  
To make all nodes in color $g$ equal, we must raise all their $d_k(v)$ up to the  **maximum**  observed value, since we can only add edges (not remove).

**Time complexity:**  $O(Nm)$.

----------

### Step 3. Compute missing edges per node

For each node $v$ of color $g$:

For each color $k \in [1,m]$:

-   Let $\delta_{v,k} = t_{g,k} - deg_{in}[v][k]$
    
-   If $\delta_{v,k} > 0$, then we must add $\delta_{v,k}$ new edges  **from nodes of color $k$**  to $v$.
    

We can pick arbitrary nodes of color $k$ (e.g. in round-robin fashion or randomly) to add the required number of edges.

**Time complexity:**  $O(Nm)$ to compute deficits.  
Edge insertion generation depends on the number of added edges, say $M_{add} = \sum_{v,k} \delta_{v,k}$. So adding edges is $O(M_{add})$.

----------

### Step 4. Add missing edges

Implementation detail:

-   Maintain for each color $k$ a list of nodes $L_k$.
    
-   When we need to add an edge from color $k$ to some target $v$, choose any node from $L_k$ (possibly allowing self-loop if $c(v)=k$.
    
-   Insert new edge(s).
    

**Time complexity:**  $O(M_{add})$.

----------

## Total Time Complexity

$$
O(E + Nm + M_{add})  
$$

In the worst case, $M_{add} \le N m \cdot \max_{v,k} t_{g,k}$, but practically $M_{add} \le O(Nm)$ for dense graphs.

----------

## Possible Optimizations

-   If the graph is sparse, use adjacency lists.
    
-   If $m \ll N$, the algorithm is near-linear.
    
-   You can parallelize Step 1 and Step 3 easily since each node’s computation is independent.