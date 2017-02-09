[THIS IS A DRAFT OF AN INCOMPLETE MANUSCRIPT, ITS CONTENTS AND ORGANIZATION DO NOT REFELCT THE QUALITY AND ACCURACY OF A PUBLICATION DRAFT, SUCH AS IT IS IT REMAINS AS A HISTORIC RECORD OF THE UNDERTAKING IT DESCRIBES AND MAY BE REVISTED BY THE AUTHOR]

# MortonCodeVsReplication
## System requirements:
### X11 development libraries
### Opengl development libraries glu, glew, gl
### Opengl 4.3 or greater 
### GPU supporting glsl 430 or greater

## Description
Buiding upon LocalDuplicationMwe, this project compares the effect of using Morton Code aka 'z-ordering' against local data replication for enhancing index locality on kernel evalution performance for n-body and multibody simulation codes. Both serial and shared memory parallel performance are compared. The shared memory implemtations are written using OpenCL 1.2 and are performed on a CPU system and a GPU system. No strong attempt at optimising the serial or parallel implementations has been done except enforcing atomic read/write operations when necessary for test subjects lacking local data replication. Obviously this harms performance for those subjects, however, this is a key point in the present discussion. Additionally the relationship between cache miss rate and execution time is analyzed.

## Background
Morton Code is the name for a method invented by an IBM researcher G.M. Morton:
  Morton, G. M. (1966), A computer Oriented Geodetic Data Base; and a New Technique in File Sequencing, Technical Report, Ottawa, Canada: IBM Ltd.
It is essentially to take the binary representations of n-dimensional geometric coordinates, interleave them, and interpret the interleaved result as the order number of the node at those coordinates. This method is used in high performance computer graphics, for example to assist in constructing 'Sparse Voxel Octrees':
  http://www.forceflow.be/2013/10/07/morton-encodingdecoding-through-bit-interleaving-implementations/
J. Baert, in above details several fast methods for computing or looking up Morton Codes, including giving links to purported SIMD implementations.
Recently there has been interest in taking advantage of how Morton Code Ordering aka "z-ordering" enhances the 'index locality' of multidimensional data when it is stored in a 1-dimensional array:
  E. Wes Bethel. zorder-lib:Library API for Z-Order Memory Layout.  Lawrence Berkeley National Laboratory. Berkeley, CA,   USA, 94720 April, 2015
Bethel, in above dicusses Morton Code Ordering (MCO), offers a library for computing such ordering, and points to research using similar techniques applied to computer math operations that indicate some performance advantage can be gained over standard linear "a-ordering" methods. Bethel also suggests that MCO is only easily applicable to cubic grids of dimensions 2 to the power n.
It is suggested in the literature that enhancing index locality decreases the rate of cache misses and is used to good effect for at least 'cache blocking' matrix operation  as in the papers cited by Bethel.

## Literature
(start with Bethel, Morton, Bethel's cites and Baert)


## Motivation
Multibody and n-body numerical simulations codes commonly feature computational kernels which operate over geometrically arranged data points, called nodes, which interact pairwise according to some some geometric rule. In an informal study, LocalDuplicationMWE, it has been suggested that for one pairwise kernel, where the geometric rule is that the radius of interaction between nodes is finite that a slight performance increase during kernel evaluation may be obtained by enhacing the index locality of node data through local data replication. 
Measurement of cache miss rates and execution times during kernel evaluation using traditional ordering and alternatively local data replication show that for small problems, a small advatange lies with local data replication.
It is speculated that as problem size grows, and radius of interaction with it, that lack of index locality becomes a greater computational burden, to be measured as cache miss rate, such that the performance improvement to be gained from local data replication may become appreciable.
certain diciplines that allow a geometric rule such as a finite radius of interaction may no see much improvement from local data replication, but it stands to reason that disciplines where such a rule is not admissible, such that more pairwise interactions over greater distances must be computed, would benifit from a way to alleviate a possible index locality issue.

## Intro
Since MCO is seen in the literature as an accepted method for enhancing index locality for improving computational performance, it provides a logical choice for comparison for determining the relative value of the local data replication method. Furthermore, the applicability of both methods, as interpreted by one of ordinary skill in computational research, to a variety of problems of interest from multiple disciplines may serve as a prediction of the performance others of similar background would obtain from employing the methods in their codes. 
One additional unique advantage offered by local data replication is the opportunity for improving shared memory parallelism since replications allows single array elements per pairwise interaction, obviating the need for atomic read/write in the shared memory parallel implementation. 

## Methods
To reduce the influences of automatic OS processes use of system resources or automatic microcode optimization for long running programs on subject evaluation, a block randomized factorial DOE with replication has been selected to generate statistically useful performance data for the comparison. 
Treatments are a-ordering (AO), z-ordering (MCO) aka Morton Code Ordering and local data replication (LDR). Problem sizes are cubic grids containing N:={2^(3*3), 2^(3*5), 2^(3*7), 2^(3*9)} nodes as well as hypercubic grids N:={2^(5*3), 2^(5*4), 2(5*5)} to demonstrate that the ordering scehemes are work for a multidimensional problem. Powers of 2 dimensions are chosen for compatibility with MCO but are not needed for AO or LDR. 
The Kernel includes an unstencyled, gaussian weighted Silling 'vector state' as in peridynamics, evaluated atypically here in that the radius of interaction is chosen such that all N nodes may interact. A constant stencyled kernel such as in image processing. taking advantage of regular grid spacing, such as would be most suitable for GPU implementation is considered beyond the scope of the present work but may prove of interest. This kernel was chosen for its assumed high computational intensity and similarity to other kernels used in n or m body problems in physics, chemistry and artificial intelligence. 
Paired t-tests are used to infer differences between the treatments in terms of their cache miss rates and execution times. The assumption of near normality of the data is bolstered by randomization, replication and the use of the lest strict t statistic. The null hypothesis in each case is that there is no difference between the treatments, problem size groups, or treament- problemsize groups. Here data is grouped first by problem size alonr, and then also by treatment and then by treatment alone without discriminating by problem size.
Correlation is computed between cache miss rate and normalized execution time. Execution time is normalized to percentile of execution times. The null hypothesis is that there is no relationship between cache miss rate and normalized execution time. Here data is grouped by problem size only.
Correlation is computed between cache miss rate amd problem size in number of nodes, N. The null hypothesis is that there is no relationship. Here data is in one encompassing group.
Correlation is computed as in both above but also in that data is additionally grouped by treatment. 
The purpose of multiple differing group comparisons is to detect whether cache miss rate and execution time are more strongly related to problem size than to specific treatment.
The program Valgrind with the tool cachegrind has been selected for measuring cache miss rate for the serial implementations. Timings are measured separately from cache miss rates as Valgrind may skew the timing results.
The implementations have been compiled into separate executables, are self-timing, and are called by a Bash script generated by a RNG code to schedule the test runs. Each test run for timing or for cache miss rate measurement consists of generating a discretization, an ordering, and evaluating the kernel for several thousand iterations. The timing or cache miss datam is then the average of what from these iterations. Timing is performed using the standard C++ 2011 monotonic clock which is suitable for measuring time intervals.
Generation of the discritization is performed by basic but suitable means commonly known in the field.
The evaluation of the shared memory implementations are similar to the serial implementations, excepting that data transfer time to and from the OpenCL device is added to execution time since it is a relevant factor. Cache miss rates are not measured for the shared memory implementations since doing so is beyond the present scope but may provide statistically verifyable explanation for timing differences.

## Results and analysis

## Future work
(If local replication is good) Having established LDR as a viable data layout strategy for simlation codes for the sake of improving performance on shared memory systems relative to AO and MCO, a parallel communication strategy will be refined for employing LDR on a heterogeneous hybrid shared-distributed memory code.

