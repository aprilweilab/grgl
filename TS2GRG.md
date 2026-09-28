# tskit TreeSequence to GRG Conversion

## The base algorithm

TODO: summarize the algorithm from the paper (with ref to paper for more detail).

## Multiple mutations at the same site

There are three scenarios for multiple mutations associated with the same site.

### Separate trees, separate alleles

The simplest scenario is that at site `s` there are two different alleles that map to two distinct sample sets. This will manifest as separate mutations with separate trees beneath them, for both tskit and GRG. There is no special handling needed. 

### Separate trees, same allele

```
        /           \
       m1           m2
       |            |
   <samples1>   <samples2>
```

A similar scenario, "recurrent mutations", has two mutations having the same allele, each with a distinct sample set. In this case, GRG wants there to be a single Mutation per `(position, allele)` combination, but the TS2GRG algorithm will create two nodes.

For each coverted local tree, we track the nodes for mutations. When we see two (or more) nodes for the same `(position, allele)`, we can simply create a parent node (in the multi-tree, not in the tree) for these nodes and map a single copy of the mutation to this parent node.

### Nested trees

When mutations for the same site, regardless of allele, are reachable as an ancestor or descendant of each other, we need slightly more complex handling. This is sometimes called a "back mutation",
when traversing down the tree one mutation "overrides" the previous one.

```
        m1
       /  \
      A    m2
            \
             B
```

Here `m1` and `m2` are mutations at the same site, and `A` and `B` both represent subtrees (and their sample sets). Typically in an ARG we consider the set of samples for a mutation `m1` to be
the sample nodes that are reachable (downwards) from the node associated with `m1`. However, with these nested/back mutations that property no longer holds. In the GRG, we require that this
property _does_ hold. We then recursively apply this procedure: 

# Converting from TS to GRG can introduce some wackiness with mutation parents. There are two scenarios
# we need to consider where there are multiple mutations at the same site.
# 1. Back mutations
#           m1
#            \
#            m2
#    Here the set of samples for m1 is samples_beneath(m1) - samples_beneath(m2). If there are further
#    child mutations below m2 then we also need to consider them. It should be true that any mutation
#    that is a child of another mutation has a subset of its samples.
#
# 2. Recurrent mutations
#
#    Here, we have two independent sets of samples that have "the same" mutation.
#
# Scenario #2 is easy: we can just create a new node that has edges to only m1, m2, and create our new
# single mutation on that node. We just have to detect this occurrence.
# Scenario #1 is trickier, but the tskit.Mutation.parent property makes it relatively simple. Whenever
# we have a 
