# Making Temporal Betweenness Computation Faster and Restless

## Network file format

The temporal network file must contains one line with the number of nodes and one line for each temporal edge. A temporal edge is a tuple `(u,v,t,l)`, where `u` and `v` are two nodes, `t` is a time in which the edge from `u` to `v` is available, and `l` is the postive traversal time of the temporal edge (if `l` is 1 for all temporal edges, then it can be omitted in the file). All networks are considered as directed: hence, if the graph is undirected, then both the temporal edge `(u,v,t,l)` and the temporal edge `(v,u,t,l)` have to be included in the file. Nodes are identified by numbers starting from 1 and times have to be integer values. The three or four elements of a temporal edge can be separated by any string, which can be specified while using the reading function (see below). In the following, we assume that this string is just a space.

The file *cannot* include duplicate lines (that is, two identical temporal edges) or self-loop lines (that is, temporal edges from a node to itself).

# How to use the software

Once the repository has been downloaded and extracted, we assume that the `julia` REPL has been started within the directory `TWBC-main` and the following command has been already executed (after having installed all the required packages).

```
include("src/Main.jl");
```

## Reading the temporal graph

The following instruction load the first temporal network of the dataset described in the first experiment of the paper, with no waiting constraint.

```
patg = read_patg("graphs/1_01_hypertext.patg", " ", β=typemax(Int64));
```

The second parameter specifies the string separating the elements of a temporal edge. By default `β` is the maximum integer value, and can be not specified (if the non-restless case is considered).

We can obtain some basic statistics about the temporal network by executing the following command.

```
print_patg_stats(patg);
```

The result in this case should be the following output.

```
====================================================
Point availability temporal graph
====================================================
Number of nodes: 113
Number of edges: 4392
Number of temporal edges: 41636
Number of distinct time steps: 5246
====================================================
```

## Computing the non-restless Sh and SFo betweenness

The values of the non-restless Sh betweenness of the hypertext temporal graph can be computed by executing the following command.

```
b, t = algorithm2("graphs/1_01_hypertext.patg", " ", 10);
```

The third parameter specifies after how many processed nodes a message has to be printed on the console (in order to verify the status of the computation). If this parameter is `0`, then no ouptut is produced. The execution of the above command should require less than one second. The values returned are the array of the SH betweenness values and the execution time.

The values of the temporal shortest betweenness can be saved as follows.

```
save_centrality_values("nrshb.txt", b);
```

Analogously, the values of the non-restless SFo betweenness of the hypertext temporal graph can be computed and saved by executing the following commands.

```
b, t = algorithm1("graphs/1_01_hypertext.patg", " ", 10);
save_centrality_values("nrsfob.txt", b);
```

## Computing the SFo betweenness 

The values of the (restless) SFo betweenness of the hypertext temporal graph with waiting constraint β equal to 600 can be computed and saved by executing the following commands.

```
b, t = algorithm3("graphs/1_01_hypertext.patg", " ", 10, _β=600);
save_centrality_values("sfob_600.txt", b);
```
## Analysing the ranking correlations

The (weighted) Kendall tau correlation of two rankings can be computed by using the [https://law.di.unimi.it](Crawdad library). The top-k intersection index can be easily computed starting form the centrality values saved in text files. 
