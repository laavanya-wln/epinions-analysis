rm(list=ls())
#install.packages("igraph")
library("igraph")
library("ggplot2")

### Original file (soc-Epinions1.txt from "https://snap.stanford.edu/data/soc-Epinions1.html") is plain text format
#data <- read.table("/Volumes/Data/Annelien/KUL-MAI/Analysis of large scale social networks/Project/soc-Epinions1.txt", header=F, colClasses=c("character", "character"))
#epinions <- graph.data.frame(data, directed=TRUE)

### For later use with Gephi and for consistency, the original file was converted into a pajek format using a Python script.
epinions <- read.graph("/Volumes/Data/Annelien/KUL-MAI/Analysis of large scale social networks/Project/soc-Epinions1.net", format="pajek")

### General info
is.directed(epinions) # Network is directed
is.weighted(epinions) # Strenght of the link: higher value=stronger link -> here, unweighted!

#extract Vertices
V(epinions)$id
V(epinions)[1] #0
vcount(epinions) #75879

#extract edges
E(epinions)[1:10]     # 0->4  0->5  0->7  0->8  0->9  0->10 0->11 0->12 0->13 0->14
# E(epinions)[1:10]$weight # unweighted network!
ecount(epinions)      # 508837
ends(epinions,10)     # incident vertices of some graph edges (for edge 10: 0 and 14)

head_of(epinions,10)  # 'head' (head of the arrow), 'tail' (the other end)
tail_of(epinions,10)
ends(epinions,1:10)   # incident vertices of some graph edges (here edges 1 to 10)

#Degree
ep_deg <-degree(epinions) # degree of a vertex: the number of adjacent edges of the vertex
#ep_indegree <- degree(epinions,mode = "in") #calculates in-degree of nodes

mean(ep_deg)  #13.4118
max(ep_deg)   #3079
min(ep_deg)   #1
median(ep_deg)#2

mean(degree(epinions,mode="out")) #6.7059
max(degree(epinions,mode="out")) #1801
min(degree(epinions,mode="out")) #0
median(degree(epinions,mode="out")) #1

mean(degree(epinions,mode="in")) #6.7059
max(degree(epinions,mode="in")) #3035
min(degree(epinions,mode="in")) #0
median(degree(epinions,mode="in")) #1

V(epinions)$name[degree(epinions,mode="in")==max(degree(epinions,mode="in"))] #To get the node with the highest in-degree - "18" (most prominent)
degree(epinions, v=which(V(epinions)$name == "18"), mode="in")  #3035
degree(epinions, v=which(V(epinions)$name == "18"), mode="out") #44

V(epinions)$name[degree(epinions,mode="out")==max(degree(epinions,mode="out"))] #Node with max out-degree - "645" - most central
degree(epinions, v=which(V(epinions)$name == "645"), mode="in")  #408
degree(epinions, v=which(V(epinions)$name == "645"), mode="out") #1801

V(epinions)$name[degree(epinions)==max(degree(epinions))] #To get the node with the highest overall degree - "18"

#Degree Distribution
deg.dist <- degree_distribution(epinions, cumulative=T, mode="all")  # numeric vector of same length as maximum degree plus one. 
plot( x=0:max(ep_deg), y=1-deg.dist, pch=19, cex=1.2, col="orange",
      +       xlab="Degree", ylab="Cumulative Frequency")
# first element: rel.freq.zero degree vertices, second: vertices with degree one, etc.
plot(degree_distribution(epinions),type="l",xlab="Degree", ylab="Degree Distribution") 
#Degree Distribution zoom
plot(degree_distribution(epinions),type='l',xlim=c(0,100),xlab="Degree", ylab="Degree Distribution")


#Closeness
ep_incloseness <- closeness(epinions,mode="in",weights=NA) #inverse of the average length of shortest paths to get TO a given node FROM all other reachable nodes in the network.
range(ep_incloseness)   #1.736851e-10 6.808817e-10
#ep_incloseness_est <- estimate_closeness(epinions,mode="in",weights=NA,cutoff=3, normalized=TRUE)
median(ep_incloseness)  #6.784961e-10
#normalised


ep_outcloseness <- closeness(epinions,mode="out",weights=NA)# #inverse of the average length of shortest paths to get FROM a given node TO all other reachable nodes.
range(ep_outcloseness)  #1.736851e-10 4.676034e-10
median(ep_outcloseness) #4.672405e-10

#Betweenness
ep_betunw<-betweenness(epinions, directed=FALSE, weights = rep(1,ecount(epinions)))
ep_betw<-betweenness(epinions, directed=TRUE)
range(ep_betw)          #0 95831448
median(ep_betw)         #0
quantile(ep_betw,probs=seq(0,1,0.25)) 
#0%      25%      50%      75%     100% 
#0        0        0    12405   95831448 
cor(ep_deg, ep_betw)    #0.7546928
plot(ep_deg, ep_betw, xlab="degree",ylab="betweenness")

#Eigenvector Centrality
epinions_undirected <- as.undirected(epinions, mode='collapse') #Eigenvector centrality gives greater weight to a node the more 
                                                                # it is connected to other highly connected nodes
                                                                # node's network importance.
ev_obj_ep <- evcent(epinions_undirected)
ep_eigen <- ev_obj_ep$vector
range(ep_eigen)   #0 1
median(ep_eigen)  #0.0003800636
quantile(ep_eigen,probs=seq(0,1,0.25)) 
#0%          25%          50%          75%         100% 
#0.000000e+00 3.968753e-05 3.800636e-04 2.094332e-03 1.000000e+00 
cor(ep_deg, ep_eigen) #0.9373242
plot(ep_deg, ep_eigen, xlab="degree",ylab="eigenvector centrality")

#Page Rank
ep_prk <- page.rank(epinions, directed=TRUE)
ep_pagerank <- ep_prk$vector
range(ep_pagerank)   #2.758026e-06 4.535031e-03
median(ep_pagerank)  #3.649017e-06
quantile(ep_pagerank,probs=seq(0,1,0.25)) 
#0%          25%          50%          75%         100% 
#2.758026e-06 2.758026e-06 3.649017e-06 7.095023e-06 4.535031e-03 

cor(ep_deg, ep_pagerank) #0.853307
cor(ep_betw, ep_pagerank) #0.5590658
plot(ep_deg, ep_pagerank, xlab="degree",ylab="PageRank")
plot(ep_pagerank,ep_betw, xlab="PageRank", ylab="betweenness")


#Table with all centralities for all vertices
centrality_epinions <- data.frame(V(epinions)$name, ep_indegree, ep_outdegree, ep_incloseness, ep_outcloseness, ep_betw, ep_eigen, ep_pagerank)

#Highest In-degree nodes
centrality_epinions[order(-centrality_epinions$ep_indegree),] #Top actors - 18(3035), 143(1521), 737(1317)

#Highest Out-degree nodes
centrality_epinions[order(-centrality_epinions$ep_outdegree),] #Top actors - 645(1801), 763(1669), 634(1621)

#Highest In-closeness nodes
centrality_epinions[order(-centrality_epinions$ep_incloseness),] #Top actors - 7437, 62983, 31013

#Highest Out-closeness nodes
centrality_epinions[order(-centrality_epinions$ep_outcloseness),] #Top actors - 69199, 69200, 69201

#Highest Betweenness nodes
centrality_epinions[order(-centrality_epinions$ep_betw),] #Top actors - 44, 763, 634

#Highest eigen centrality
centrality_epinions[order(-centrality_epinions$ep_eigen),] 
#Top actors - 18, 645, 634, 44... check other nodes in top 10

#Highes Pagerank centrality
centrality_epinions[order(-centrality_epinions$ep_pagerank),]
#Top actors - 18, 737, 118

# Generate a table of pairwise correlations of all centrality measures
correlations_centrality <- cor(centrality_epinions[,2:7])


#Density
graph.density(epinions) #8.83774e-05-- very sparse

#Global Clustering Coefficient
transitivity(epinions) #0.06567883

#Average Clustering Coefficient
transitivity(epinions, type = "average") #0.2605128

#Average Degree
mean(ep_deg) #13.4118
mean(ep_indegree) #6.7059
mean(ep_outdegree) #6.7059

#Average Path Length
mean_distance(epinions, directed=T) #4.754723


### Clustering
clu_ep_weak = clusters(epinions, mode="weak")
clu_ep_weak$no #only 2 clusters are formed
max(clu_ep_weak$csize) #75877
min(clu_ep_weak$csize) #2 (the second cluster only has two members, all the other nodes are in the first cluster)
which(clu_ep_weak$csize==max(clu_ep_weak$csize)) #cluster 1 is the largest cluster
index_largest_component=which(clu_ep_weak$csize==max(clu_ep_weak$csize))
which(clu_ep_weak$membership==index_largest_component) #which nodes are in the largest cluster
index_smallest_component=which(clu_ep_weak$csize==min(clu_ep_weak$csize))
which(clu_ep_weak$membership==index_smallest_component) 
# Not a very useful clustering (only 2 clusters, one of which has only 2 members)

clu_ep_strong = clusters(epinions, mode="strong") #based on strongly connected components (i.e. nodes that are fully connected with each other)
clu_ep_strong$no #42176
max(clu_ep_strong$csize) #32223
min(clu_ep_strong$csize) #1
which(clu_ep_strong$csize==max(clu_ep_strong$csize)) #cluster 27256

test <- data.frame(numberclu=clu_ep_strong$no, sizeclu=clu_ep_strong$csize)
ggplot(test, aes(x=numberclu,y=sizeclu)) +
  geom_point()
table(test$sizeclu) # A lot of clusters with only one member, 1 cluster with (32223/75879) 42.5% of the nodes (that means that all these nodes are fully connected with each other! huge!)
#1      2     3     4     5     6     7     8     9    15   32223 
#41112  813   164   47    24    5     1     6     2     1     1 

index_largest_component=which(clu_ep_strong$csize==max(clu_ep_strong$csize))
which(clu_ep_strong$membership==index_largest_component)

which(clu_ep_strong$csize==2)
sizes(clu_ep_strong)[sizes(clu_ep_strong)==2]

#red_ep=delete.vertices(epinions, V(epinions)[clu_ep_strong$membership!=index_largest_component])
#vcount(red_ep) #32223
#ecount(red_ep) #443506

red_ep=induced.subgraph(epinions, which(clu_ep_strong$membership==index_largest_component))
vcount(red_ep) #32223
ecount(red_ep) #443506

write_graph(red_ep, "Reduced_graph_largestcluster.net", format = "pajek")


### Another kind of clustering for the entire network
LE_ep=cluster_leading_eigen(epinions_undirected)
membership(LE_ep)
length(LE_ep) #5
sizes(LE_ep)
# 1     2     3     4     5 
#15948     2  4182  9512 46235 

sizes(LE_ep)[sizes(LE_ep)>10] #clusters with more than 10 members
modularity(epinions_undirected, membership(LE_ep)) #0.3126266

V(epinions)$color=membership(LE_ep)
plot(epinions, layout=layout_with_drl, vertex.label=NA, vertex.size=4)

### Here we compare two clustering methods on the reduced network (the largest strongly connected cluster from above)
LE_redep=cluster_leading_eigen(red_ep)
membership(LE_redep)
length(LE_redep) #5
sizes(LE_redep)
#Community sizes
#1     2     3     4     5 
#6654  3895  3623 16560  1491 

sizes(LE_redep)[sizes(LE_redep)>10] #all clusters
modularity(red_ep, membership(LE_redep)) #0.2998961

LV_redep=cluster_louvain(redep_undirected)
length(LV_redep) #423
sizes(LV_redep)
modularity(red_ep, membership(LV_redep)) # 0.4251552

compare(LE_redep, LV_redep,'adjusted.rand') # 0.1211313

V(red_ep)$color=membership(LV_redep)
plot(red_ep, layout=layout_with_drl, vertex.label=NA, vertex.size=4)



