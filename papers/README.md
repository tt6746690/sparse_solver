https://sympiler.github.io/


### basics...


+ tiling (loop nest optimization)
    + loop transformations for locality/parallelism


+ compressed sparse column (CSC) format
    + https://www.zhihu.com/search?type=content&q=sparse%20matrix%20

```c
/* --- primary CSparse routines and data structures ------------------------- */
typedef struct cs_sparse    /* matrix in compressed-column or triplet form */
{
    csi nzmax ;     /* maximum number of entries */
    csi m ;         /* number of rows */
    csi n ;         /* number of columns */
    csi *p ;        /* column pointers (size n+1) or col indices (size nzmax) */
    csi *i ;        /* row indices, size nzmax */
    double *x ;     /* numerical values, size nzmax */
    csi nz ;        /* # of entries in triplet matrix, -1 for compressed-col */
} cs ;

// Matrix:
// 1，2，0，0；
// 7，0，0，4；
// 0，0，0，1；
// 0，9，8，1。

// p：0，2，4，5，8；
// i：0，1，0，3，3，1，2，3；
// x：1，7，2，9，8，4，1，1。
```

+ pivot element
    + in Gauss-Jordan elimination, the diagonal of  a matrix is called a pivot element, which has to be a nonzero element. 
+ partial pivoting
    + interchanging rows to place "good" (i.e. large, nonzero) elements in the diagonal position prior to forward/backward substitution

+ elimination tree
    + https://graal.ens-lyon.fr/~bucar/CR07/lecture-etree.pdf

+ polyhedral model 
    + http://icps.u-strasbg.fr/~bastoul/research/papers/Bas04-PACT.pdf
    + http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.73.7141&rep=rep1&type=pdf

+ https://www3.nd.edu/~zxu2/acms60212-40212-S12/Lec-05.pdf
    + task granularity
        + fine grained
        + coarse grained means small number of takss
        + granularity has a finite bound
    + degree of concurrency
        + number of tasks that can be executed in parallel
    + critical path (weighted)
        + longest directed path from start node to finish node
    + average degree of concurrency
        + total work / critical path length

### sparse matrix algorithms

+ [1986_sparse_partial_pivoting_in_time_proportional_to_arithmetic_operations](1986_sparse_partial_pivoting_in_time_proportional_to_arithmetic_operations.pdf)
    + idea
        + determine nonzero structure
        + use this info to avoid unnecessary arithmetics



### domain-specific compiler for sparse methods


+ [2017_sympiler_transforming_sparse_matrix_codes_by_decoupling_symbolic_analysis](2017_sympiler_transforming_sparse_matrix_codes_by_decoupling_symbolic_analysis.pdf)
    + idea
        + sympiler: a sparsity-aware code generator for sparse matrix algos that leverages symbolic information to generate fast code for a specific matrix structure
    + motivation 
        + sparse matrix operation 
            + indirect array accesses
            + hard to apply tiling/vectorization
        + traditional method
            + specialized __library__ with hand-tuned implementation (CHOLMOD) for different architectures
                + bad: must be manually ported to new architecture
        + better approach
            + compiler optimization
                + good: handles portability for you
                + bad: still has indirect acceses, complex dependence structure
        + compiler loop optimization 
            + polyheral model 
                + _5,10,41,54,65,67_
                    + problem with non-affine loop bounds and/or array subscript
                + _63,66,68-70_
                    + extended to operate on kernels with static index arrays (inspector-executor)
                    + limited to transforming sparse kernels with static index arrays
            + sympiler
                + symbolic analysis at compile time
                + remove dynamic index arrays -> only static index array is left
                + constraint:
                    + matrix nonzero remains constant
        + sparse matrix method (LU, cholesky)
            + view computation as a graph, then apply method-dependent graph algorithm gives information about dependencies, which is used to more efficiently compute numerical method
            + 14: libraries utilize symbolic information but this is __coupled__ with numerical computation. -> makes compiler optimization difficult
    + sympiler
        + generates high-perf sparse matrix code
            + decouple symbolic analysis from numerical computation
        + from symbolic info, applies inspector-guided transformations, e.g. variable-sized blocking, resulting in similar perf with hand-tuned library
        + also generates code for a specific nonzero structure (inspector-guided, low lvl transformations)
        + improves perf by applying optimization for single-core with _vectorization_ and _increased data locality_ 
        + extension: shared/distributed memory system
    + motivating example
        + `Lx=b` where `L,b` are sparse
        + library
            + skip values where `b` is sparse
            + `O(|b| + n + f)` where `|b|` is number of nonzero elements in `b` and `f` is the number of floating point operations
        + decoupled code
            + dependency graph DG_L = (V,E) computes nonzero pattern of `x`
            + _34_: only needed to compute on nodes in reachset (by dfs) starting at nonzeros of `b` 
                + eigen actually does not use this (https://eigen.tuxfamily.org/dox/TriangularSolver_8h_source.html)
            + `O(|b|+ f)`
        + sympiler
            + use dependency graph idea but generates reach-set for a given `L,b` at compile time 
            + properties
                + iterates over reached column
                + peels iteration where number of nonzeros in a column is greater in some threshold
                + peeled loop -> vectorization
            + central idea
                + decouple symbolic analysis from numerical computation
                + use reachset to guide transformation
            + perf
                + SuiteSparse ...
                + avg `1.3x` compared to library
    + static sparsity pattern
        + central assumption:
            + nonzero pattern is fixed during compile time
        + centrla idea:
            + values might change, but nonzero pattern change rarely ...
            + sparsity pattern in different domains does not change in many applications...
        + examples
            + power system modeling, circuit simulation 
                + interconnection rarely change, unless a circuit breaks ...
            + graphics
                + sparse structure originates from physical discretization and sparsity pattern is same, except where there are deformations ...
    + contributions
        + compile-time symbolic inspectors
        + inspector-guided transformations, leveraging compile-time information to transform sparse matrix code
        + implementation of inspector/inspector-guided transformations for 
            + sparse triangular solve
            + sparse cholesky factorization
        + better performance ...
    + overview
        + method-specific symbolic inspector
        + domain-specific optimizations while lowering the code
        + lowered code annotated with additional low-level transformations (unrolling)
        + annotated code apply low-level optimization and output to C source code
        + solver has to be implemented in domain-specific AST
    + symbolic inspector
        + specific w.r.t. 
            + sparse algorthm
            + transformation
        + steps
            + creates inspection graph
                + e.g. directed dependency graph 
            + traverse during inspection strategy
                + e.g. dfs over inspection graph
            + inspection set guides transformations
                + e.g. reach-set guides loop prune
    + inspector-guided transformations
        + variable iteration space pruning (VI-Prune) 
            + prunes iteration space of a loop using info about sparsity info
        + 2D variable-sized blocking
            + _44_: supernodes
            + challenges
                + block sizes variable
                    + keep record of the size of each block in inspection set
                + CSC formats, block elements maynot be consecutive memory locations
                    + runtime: transformed code allocate temp block storage and coppies data prior to operating on the block
                + numerical method used may be changed later
                    + different procedure chosen for each loop/instruction
    + enabled conventional low-level transformations
        + loop peeling
        + removes indirect memor accesses
        + vectorization of loops (since loop boundary known at compile time)
    + Sparse triangular solve
    + cholesky factorizatoin 
        + theory
            + _17,14,52_: elimination tree
    + other matrix method
        + sympiler can support many other methods
        + elimination tree, reachset, supernodes, dependency graph are fundamental symbolic analysis for optimized algorithms such as rank update and rank increase method
    + performance 
        + cholesky
            + https://www.sandia.gov/~srajama/publications/cholmod_toms.pdf
            + libraries cannot afford to have separate implementation for each sparsity pattern and do not implement sparsity-specific optimization (?)
        + when taken into account symbolic analysis, the runtime is 
            + triangular: sympiler is worse than eigen
            + cholesky: faster than eigen and CHOLMOD
        + sympiler  
            + can generate specialized and efficient codes for small dense sub-kernels
    + related work  
        + _63,66,68-70_: develoepd compile time techniques for automatically creating inspector/executors for use at runtime
            + both runtime
            + limited to static index arrays
            + OK for sparse LU, ... where additional nonzeros/fill-ins are not introduced during computation; but limited for sparse matrix method that modifies the array dynamically
        + domain specific compiler
            + _1,48_: graphics
                + code generation for pdes and finite element methods
            + _16,56_: sparse operations
                + https://docs.google.com/viewer?url=https%3A%2F%2Fgithub.com%2FIntelLabs%2FSparso%2Fblob%2Fmaster%2Fdoc%2FSparso%2520PACT16.pptx%3Fraw%3Dtrue
                + (taco) http://tensor-compiler.org/kjolstad-oopsla17-tensor-compiler.pdf
                + sympiler is first compiler that does symbolic analysis at compile time and buids specialized sparse solvers


+ potential problems
    + for triangular solver, both `L` and `b` has to be sparse
    + symbolic analysis at compile time, requiring sparse pattern to be fixed, which is a huge constraint
        + maybe the analysis phase is too long to make this feasible at runtime 
        + how do you check that the matrix has the same sparsity pattern, if have to check everytime you run the algorithm, then what is the benefit of symbolic analysis at compile time
        + examples not so convincing
            + powerline: compile everytime there is a new circuit?
            + graphics: probably want algorithm to work for any mesh of say different topological properties, not just one fixed type of object. limits the usage i guess.
    + loop peeling might bloat generated code for large sparse matrices
    + easy to extend to other sparse algorithms, does assumptions about static sparsity pattern remain applicable?
    + method specific symbolic inspector...
    + have to write numerical algorithms on AST, which may not be easy


The main idea of Sympiler is to move symbolic analysis (as part of a sparse algorithm) to compile time. This gies 


strength:

1. Sympiler is a domain-specific compiler that does symobolic analysis for sparse methods at/before compile time. The motivation for doing this is to better utilize existing compilers (e.g. gcc) for portability and low-level optimization. This is a major strength and novelty of the system
    
    * Compared with library implementation that is hand-tuned, often specifically for an architecture, Sympiler leverage compiler for portability. 
    * Compared with methods that uses inspector-executor model at runtime, inspection/transformation for Sympiler is done at compile time, which allows for additional compiler optimizations (loop unrolling, vectorization). This is not possible for methods whereby the inspection is done at runtime.

    The second point is, by looking at Figure 6. and more so at Figure 7, the single biggest reason why Sympiler's generated code outperforms other state-of-the-art methods.

2. The performance gain is sizable. 

weaknesses:

1. The assumption of static sparsity pattern limits the use case for Sympiler. The assumption, to be more precise, is that location for nonzero values has to be known during compile time. This is a big assumption

    * The paper mentions circuit simulation as a motivating example. Although I am not really familiar with hardware, it seems that to re-generate and re-compile the code for each different circuit configuration is not ideal. 
    * Additionally, it is probably not ideal to recompile a graphics program everytime a different mesh, or the same mesh but with different vertex locations are used. 
    
    For the sake of experimentation, this assumption weight heavily against using Sympiler. For the sake of performance in deployed code, the assumption of static sparsity pattern limits its application. For example, MKL underlies Numpy/Scipy stack in python but Sympiler cannot replace MKL in this case due to this assumption. The usage has to be deliberate and so will limit its use case.

2. Generated code might increase instruction cache miss

3. Performance is highly dependent on the sparsity pattern, decreasing the system's robustness

    * VS-BLOCK, VI-PRUNE might have no tangible benefit for specific scenarios, excluding low-level transformations. For example, in Figure 6. VS-BLOCK does not yield performance gain when matrix do not have supernodes; VI-PRUNE does not yield performance benefit if data dependency mandates the computation of majority of column updates in triangular solve.

4. The claimed benefit (little hand-tuning, easy maintenance and implementation) of Sympiler does not hold in certain cases.

    * For triangular solve, the library implementation (Eigen: https://eigen.tuxfamily.org/dox/TriangularSolver_8h_source.html#l00117) is as simple as a `if` statement that skips `i`-th iteration if the value of `rhs[i]` is 0. There is almost no cost to maintaining this code, no cumbersome hand-tuning, and is not architecture specific. Sympiler makes above-mentioned points more difficult, with the added requirement of known sparsity pattern during compilation, and yields a ~1.5x increase in performance when counting just numeric computation and decreases the performance if taken into account both symbolic and numeric computation. The ~1.5x increase in performance is largely due to low-level optimization of existing compilers. This somewhat questions the significance of Sympiler and if the advantages that the paper claim holds more generally. (Or if triangular solve is simply shown for the sake of explanation)

5. more of a question: 

    * The paper states that CHOLMOD uses BLAS for small sub-kernel and so does not perform well compared to Sympiler, which generates dense code for such sub-kernels. The question is then if CHOLMOD, which has knowlege of column counts, simply implements the heuristic that Sympiler does, it would probably have similar performance. It is true that Sympiler generates kernel-specific optimization. But I somewhat get the idea that CHOLMOD has the choice but did not implement the heuristic. Therefore, the paper can't really say that performance gain from such heuristic is because of the fact that symbolic analysis and numerical computation are decoupled. 

extensions: 

1. As noted in the paper, extend Sympiler to support 

    * more sparse matrix formats
    * more sparse methods
    * shared (ParSy) and distributed memory systems

2. From looking at Figure 7., it seems that there is still room for optimization to VS-BLOCK (excluding low-level optimization), as Sympiler:VS-BLOCK is considerably worse than CHOLMOD's implementation that also utilizes supernodes. 

3. cool idea: Sympiler has row-sparsity information in prune-set and therefore does need to compute `A` transpose in Eigen and CHOLMOD. There might be a lot of potential along this line of work: doing symoblic analysis before numerical computation might simplify sparse algorithms algorithmically. So it may be a good idea to survey other sparse methods that can be simplified










+ [2018_parsy_inspection_and_transformation_of_sparse_matrix_computations_for_parallelism](2018_parsy_inspection_and_transformation_of_sparse_matrix_computations_for_parallelism.pdf)
    + abstract
        + parsy is a framework for generating parallel code for sparse matrix computation
        + inspection+transformation
        + optimized for locality and load balance
            + load balance is something that traditional algorithm do not perform well
        + task coarsening strategy 
            + for well-balanced tasks for parallel execution 
        + perf
            + 3x over MKL pardiso and PaStiX
    + introduction 
        + previous work
            + _54,45,59_: runtime inspector - executor
            + wavefront parallelism:
                + traverse DAG create level sets, i.e. iterations that can be executed in parallel
                + synchronization between level sets ensure data dependencies 
                + however high overheads with number of levels: lead to load imbalance for cholesky, with non-uniform workload
        + parsy
            + new _load-balanced level coarsening (LBC) algorithm_ 
                + on dependence graph to create well-balanced coarsened level sets: hierarchical level set (H-level set)
            + a novel proportional cost model including LBC
                + tradeoff between locality, load balance, parallelism
        + algorithm
            + l-partitioning
            + w-partitioning
                + _25,10_: first fit decreasing bin packing 
    + performance
        + parsy has less CPU utilization, but better cache locality
    + dynamic scheduling introduce overhead at runtime; pastix uses static analysis but is architecture dependent; parsy does static analysis at compile time ... so architecture specific; compile time scheduling of the tasks


+ problems
    + numerical algos locality and synchronization vs. CPU utilization






Hi prof. Dehnavi, 

I am interested in joining your group (reason why I put your name in the application). But I wanted to let you know that at the moment I am somewhat leaning towards graphics as I mentioned in my personal statement. However, I do find making fast scientific computing tools impactful and cool, especially related to algorithms and distributed system. Another thing is I applied to MSc over direct-PhD because I saw on the website that it is in general easier to get into UofT this way. I'd be happy with a direct-PhD for sure but I did not apply for one. I just want to let you know the situation since I cannot change the application now.



Also sorry about the late rely. I originally thought I would have to review both papers and just found out earlier today that I needed to review one of them after reading both. I I have pasted a short view in this email and I hope it is not too general. I spent most time understanding the paper and did not read in detail other papers in the area. Regardless the outcome of the review and interview, would you please let me know what you think about the Sympiler paper, briefly. I am genuinely curious on what your opinions are.

Sincerely,

peiqi



Review for the Sympiler paper:


strength:

1. Sympiler is a domain-specific compiler that does symobolic analysis for sparse methods at/before compile time. The motivation for doing this is to better utilize existing compilers (e.g. gcc) for portability and low-level optimization. This is a major strength and novelty of the system
    
    * Compared with library implementation that is hand-tuned, often specifically for an architecture, Sympiler leverage compiler for portability. 
    * Compared with methods that uses inspector-executor model at runtime, inspection/transformation for Sympiler is done at compile time, which allows for additional compiler optimizations (loop unrolling, vectorization). This is not possible for methods whereby the inspection is done at runtime.

    The second point is, by looking at Figure 6. and more so at Figure 7, the single biggest reason why Sympiler's generated code outperforms other state-of-the-art methods.

2. The performance gain is sizable. 

weaknesses:

1. The assumption of static sparsity pattern limits the use case for Sympiler. The assumption, to be more precise, is that location for nonzero values has to be known during compile time. This is a big assumption

    * The paper mentions circuit simulation as a motivating example. Although I am not really familiar with hardware, it seems that to re-generate and re-compile the code for each different circuit configuration is not ideal. 
    * Additionally, it is probably not ideal to recompile a graphics program everytime a different mesh, or the same mesh but with different vertex locations are used. 
    
    For the sake of experimentation, this assumption weight heavily against using Sympiler. For the sake of performance in deployed code, the assumption of static sparsity pattern limits its application. For example, MKL underlies Numpy/Scipy stack in python but Sympiler cannot replace MKL in this case due to this assumption. The usage has to be deliberate and so will limit its use case.

2. Generated code might increase instruction cache miss

3. Performance is highly dependent on the sparsity pattern, decreasing the system's robustness

    * VS-BLOCK, VI-PRUNE might have no tangible benefit for specific scenarios, when low-level transformations are disabled. For example, in Figure 6. VS-BLOCK does not yield performance gain when matrix do not have supernodes; VI-PRUNE does not yield performance benefit if data dependency mandates the computation of majority of column updates in triangular solve.

4. The claimed benefit (little hand-tuning, easy maintenance and implementation) of Sympiler does not hold in certain cases.

    * For triangular solve, the library implementation (Eigen: https://eigen.tuxfamily.org/dox/TriangularSolver_8h_source.html#l00117) is as simple as a `if` statement that skips `i`-th iteration if the value of `rhs[i]` is 0. There is almost no cost to maintaining this code, no cumbersome hand-tuning, and is not architecture specific. Sympiler makes above-mentioned points more difficult, with the added requirement of known sparsity pattern during compilation, and yields a ~1.5x increase in performance when counting just numeric computation and decreases the performance if taken into account both symbolic and numeric computation. The ~1.5x increase in performance is largely due to low-level optimization of existing compilers. This somewhat questions the significance of Sympiler and if the advantages that the paper claim holds more generally. (Or if triangular solve is simply shown for the sake of explanation)

5. more of a question: 

    * The paper states that CHOLMOD uses BLAS for small sub-kernel and so does not perform well compared to Sympiler, which generates dense code for such sub-kernels. The question is then if CHOLMOD, which has knowlege of column counts, simply implements the heuristic that Sympiler does, it would probably have similar performance. It is true that Sympiler generates kernel-specific optimization. But I somewhat get the idea that CHOLMOD has the choice but did not implement the heuristic. Therefore, the paper can't really say that performance gain from such heuristic is because of the fact that symbolic analysis and numerical computation are decoupled. 

extensions: 

1. As noted in the paper, extend Sympiler to support 

    * more sparse matrix formats
    * more sparse methods
    * shared (ParSy) and distributed memory systems

2. From looking at Figure 7., it seems that there is still room for optimization to VS-BLOCK (excluding low-level optimization), as Sympiler:VS-BLOCK is considerably worse than CHOLMOD's implementation that also utilizes supernodes. 

3. cool idea: Sympiler has row-sparsity information in prune-set and therefore does need to compute `A` transpose in Eigen and CHOLMOD. There might be a lot of potential along this line of work: doing symoblic analysis before numerical computation might simplify sparse algorithms algorithmically. So it may be a good idea to survey other sparse methods that can be simplified

