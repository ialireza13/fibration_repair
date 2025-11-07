
## Problem Definition

We are given:

-   A  **directed graph**  ( G = (V, E) ) with ( |V| = N ).
    
-   Each node ( v \in V ) has a  **color (group)**  ( c(v) \in {1, 2, \dots, m} ), with ( m \le N ).
    
-   For each node ( v ), define its  **color in-degree vector**:
    

[  
\text{deg_in}(v) = [d_1(v), d_2(v), \dots, d_m(v)]  
]

where ( d_k(v) ) = number of incoming edges to ( v )  **from nodes of color ( k )**.

The  **coloring is valid**  if for any two nodes ( u, v ) with the same color ( c(u) = c(v) ), we have:

$$deg_{in}(u) = deg_{in}(v) $$

##  Goal

Add a minimal set of edges to make this property true — i.e., to make every group have identical color-in-degree vectors across its members.  
(Since the source nodes are only constrained by color, not by identity, we can choose  _any_  nodes of a given color to provide those missing edges.)

----------

## 🧠 Algorithm

### Step 1. Compute color–color in-degree matrix per node

For each node ( v ):

1.  Initialize an array ( \text{deg_in}[v][1..m] = 0 ).
    
2.  For each incoming edge ( (u \to v) ):
    
    -   Increment ( \text{deg_in}[v][c(u)]++ ).
        

This gives every node’s color-in-degree vector.

**Time complexity:**  ( O(|E|) ).

----------

### Step 2. Compute target in-degree vector per color group

For each color ( g = 1..m ):

1.  Consider all nodes ( v ) with ( c(v) = g ).
    
2.  For each color ( k = 1..m ):
    
    -   Compute ( t_{g,k} = \max_{v: c(v)=g} \text{deg_in}[v][k] ).
        

This defines the  **target color-in-degree vector**  for group ( g ):  
[  
T_g = [t_{g,1}, t_{g,2}, ..., t_{g,m}]  
]

Justification:  
To make all nodes in color ( g ) equal, we must raise all their ( d_k(v) ) up to the  **maximum**  observed value, since we can only add edges (not remove).

**Time complexity:**  ( O(Nm) ).

----------

### Step 3. Compute missing edges per node

For each node ( v ) of color ( g ):

For each color ( k \in [1,m] ):

-   Let ( \delta_{v,k} = t_{g,k} - \text{deg_in}[v][k] )
    
-   If ( \delta_{v,k} > 0 ), then we must add ( \delta_{v,k} ) new edges  **from nodes of color ( k )**  to ( v ).
    

We can pick arbitrary nodes of color ( k ) (e.g. in round-robin fashion or randomly) to add the required number of edges.

**Time complexity:**  ( O(Nm) ) to compute deficits.  
Edge insertion generation depends on the number of added edges, say ( M_{add} = \sum_{v,k} \delta_{v,k} ). So adding edges is ( O(M_{add}) ).

----------

### Step 4. Add missing edges

Implementation detail:

-   Maintain for each color ( k ) a list of nodes ( L_k ).
    
-   When we need to add an edge from color ( k ) to some target ( v ), choose any node from ( L_k ) (possibly allowing self-loop if ( c(v)=k )).
    
-   Insert new edge(s).
    

**Time complexity:**  ( O(M_{add}) ).

----------

## ⏱️ Total Time Complexity

Step

Complexity

Step 1 (counting)

( O(E) )

Step 2 (target vectors)

( O(Nm) )

Step 3 (compute deficits)

( O(Nm) )

Step 4 (add edges)

( O(M_{add}) )

Overall:

[  
O(E + Nm + M_{add})  
]

In the worst case, ( M_{add} \le N m \cdot \max_{v,k} t_{g,k} ), but practically ( M_{add} \le O(Nm) ) for dense graphs.

----------

## ✅ Correctness Justification

-   After repair, for every color ( g ), every node ( v ) with ( c(v)=g ) has:  
    [  
    \text{deg_in}[v][k] = t_{g,k} = \max_{v': c(v')=g} \text{deg_in}[v'][k]  
    ]
    
-   Thus all nodes in color ( g ) share the same color-in-degree vector.
    
-   No unnecessary edge removals are done (only additions).
    
-   The number of added edges is  **minimal**, since each node only receives exactly enough edges to reach the maximum per color.
    

----------

## ⚙️ Example

Let’s revisit your example:

Node

Color

In from c1

In from c2

In from c3

i

1

2

3

0

j

1

1

2

1

Then for color group 1 (nodes i, j):

-   Target vector ( T_1 = [\max(2,1), \max(3,2), \max(0,1)] = [2,3,1] )
    

Missing edges:

-   For i: +1 edge from color 3
    
-   For j: +1 from color 1, +1 from color 2
    

Exactly as desired.

----------

## 🧮 Possible Optimizations

-   If the graph is sparse, use adjacency lists.
    
-   If ( m \ll N ), the algorithm is near-linear.
    
-   You can parallelize Step 1 and Step 3 easily since each node’s computation is independent.
    

----------

## 🧾 Summary

**Algorithm:**

1.  Compute color-in-degree per node.
    
2.  For each color group, determine max per incoming color.
    
3.  For each node, compute deficit and add edges accordingly.
    

**Time Complexity:**  ( O(E + Nm + M_{add}) )  
**Space Complexity:**  ( O(Nm) )  
**Guarantees:**

-   Produces minimal edge additions.
    
-   Balances color-in-degree vectors per color group.
    
-   Works with self-loops.
    

----------

Would you like me to give you a  **Python implementation**  of this algorithm next (that works with adjacency lists and color arrays)?