<h3>
	<center>The Most Influenced Community Search on Social Networks</center>
</h3>

### Introduction

1. This repository contains the technical report of our paper.
2. This repository contains the codes and datasets used in our paper.

### Datasets

We use 8 publicly available real-world networks, including EmailCore, WikiVote, Epinions, Slashdot, EmailEuall, Pokec,  LiveJournal, and WikiLink.

WikiLink is available on Konect [1].

The Others are available on SNAP [2].

*[1] Jérôme Kunegis. 2013. Konect: the koblenz network collection. In Proceedings of the 22nd international conference on world wide web. 1343–1350.*

*[2] Jure Leskovec and Andrej Krevl. 2014. SNAP Datasets: Stanford large network dataset collection. http://snap.stanford.edu/data.*

### Algorithms

The following files are the codes for our proposed algorithms. All methods are implemented in C++ and executed on a server with Intel(R) Core(TM) 3.70GHz CPU and 128GB of RAM.

#### In InfExpectation_computation folder 

1. infseed.cpp: computing singleton influence of each node, and sort them in descending order of influence.

```shell
./InfExpectation_computation --r 10000 --dataset <dataset root> --model IC
```

2. tim.cpp: selecting k seed nodes with the maximum influence spread using TIM algorithm [3]

*[3] Youze Tang, Xiaokui Xiao, and Yanchen Shi. 2014. Influence maximization: Nearoptimal time complexity meets practical efficiency. In Proceedings of the 2014 ACM SIGMOD International Conference on Management of Data. 75–86.*

```shell
./InfExpectation_computation -dataset <dataset root> -model IC  -seedset <seedset output root>  -epsilon <sampling approximation loss> -k <number of seeds>
```

3. seed_gen.cpp: generating seed nodes using Random/Inf methods

```shell
./InfExpectation_computation --method random --k <number of seeds> --dataset <dataset root> --seedset <seedset output root>
```

```shell
./InfExpectation_computation --method inf --k <number of seeds> --f <fraction> --dataset <dataset root> --seedset <seedset output root>
```

4. InfCompute.cpp: Compute influenced expectation of each node under a seed node set (**S-InfExp algorithm**).

```shell
./InfExpectation_computation --theta 10000 --dataset <dataset root> --model IC --mode <IM/Random/Inf> --seedset <seedset root>  --seednum <number of seeds>
```

#### In MLPred folder

Predicting influenced expectations (**L-InfExp algorithm**).

#### In InfCS folder

Searching for the MIC with **GlobalSearch** and **LocalSearch** algorithm. (SFIS: S-InfExp algorithm, ONPR: L-InfExp algorithm)

```shell
./InfCS -dataset=<dataset root> -subgraph=<subgraph root for scalability test> -model=IC -mode=<IM/Random/Inf> -func=<avg/min/max/sum> -seed=<number of seeds> -k=<indegree constraint> -l=<outdegree constraint> -s=<size constraint of MIC> -alg=approx -med=<global/local> -inf=<SFIS/ONPR> -p=<percentage of nodes for scalability test>
```


### Running Environment

A 64-bit Linux-based OS. 

GCC 4.7.2 and later.
