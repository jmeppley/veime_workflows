from itertools import chain
import json
import pandas
import igraph
import leidenalg
output_data = {}

input = snakemake.input
output = snakemake.output

# load gene data, but only keep annotated genes
annot_col = "annot_e"
gene_table = pandas.read_csv(str(input), sep='\t', index_col=0)
gene_table = gene_table[gene_table[annot_col].notna()]

# get the annotations
n_uniq_annots = None

last_annot_count = 0
last_genome_count = 0
annots = set(gene_table[annot_col].values)
genomes = set(gene_table.genome.values)

# iteratively filter genomes and genes until we stabilize
while (
    len(annots) != last_annot_count 
    and len(genomes) != last_genome_count
):
    last_annot_count = len(annots)
    last_genome_count = len(genomes)

    # require annots to be in at least 2 genomes
    # require genomes to have at least 2 annots found in other genomes
    annot_counts = \
        gene_table[[g in genomes for g in gene_table.genome]] \
            .groupby(['genome', annot_col]) \
            .agg({'gene_no': len}) \
            .reset_index() \
            .groupby(annot_col) \
            .agg({'genome':len}) \
            .sort_values('genome', ascending=False)

    if n_uniq_annots is None:
        n_uniq_annots = annot_counts.shape[0]

    # require annots to be in at least 2 genomes
    annots = set(annot_counts.query('genome > 1').index)

    # require genomes to have at least 2 annots found in other genomes
    genomes = set(g for g, df 
                  in gene_table[[(a in annots) 
                                 for a 
                                 in gene_table[annot_col].values]] \
                          .groupby('genome')
                  if df.shape[0] > 1)

# only keep genomes with 2 annots found elsewhere
named_edges = set(tuple(v) 
                  for v in gene_table[['genome',annot_col]].values
                  if v[1] in annots and v[0] in genomes)

output_data['named_edges'] = tuple(named_edges)

# count some things
graphed_seqs = set(n[0] for n in named_edges)
graphed_annots = set(n[1] for n in named_edges)
n_seqs, n_annots = len(graphed_seqs), len(graphed_annots)

n_genomes = len(set(gene_table.genome))
n_annotated = len(set(gene_table.genome))

output_data['counts'] = dict(
    annot_col=annot_col,
    starting_genomes=n_genomes,
    annotated_genomes=n_annotated,
    used_genomes=n_seqs,
    starting_annotations=n_uniq_annots,
    used_annotations=n_annots,
)

print(f'''
column: {annot_col}
# genomes: {n_genomes}
 # annotated: {n_annotated}
 # in edges: {n_seqs}
# annots {n_uniq_annots}
 # in edges {n_annots}

Did we lose any genomes with non-null annotations?
 == {n_genomes - n_annotated} lost PRs because they don't have any annotated genes
 == {n_annotated - n_seqs} lost PRs because none of their annoated genes appear in any other genomes
''')

## Build graph and calcuate modules
# save lookup dict for sving module info 
node_ids = {i:n for i,n in enumerate(chain(graphed_seqs, graphed_annots))}
node_id_dict = {n:i for i,n in node_ids.items()}

#### build graph
numbered_edges = [tuple(node_id_dict[n] for n in e) for e in named_edges]
graph = igraph.Graph(edges=numbered_edges)

#### Label layers
for btype, nameset in ((0, graphed_annots), (1, graphed_seqs)):
    for name in nameset:
        node_id = node_id_dict[name]
        graph.vs[node_id]['type'] = btype

if not graph.is_bipartite():
    raise Exception("Graph is not bipartite somehow!!")


#### Get largest connected subgraphs
ccss = graph.components(mode='strong')
output_data['counts']['n_connected_components'] = len(ccss)
if not graph.is_connected():
    print('subgraph sizes: ', [len(sg.vs) for sg in ccss.subgraphs()])

    core_graph = ccss.subgraphs()[0]
else:
    print('already connected!')
    core_graph = graph

#### Find modules
optimiser = leidenalg.Optimiser()
memberships_by_resolution = {}
for r in [2, 1.5, 1, .75, .5, .25]:
    p_01, p_0, p_1 = leidenalg.CPMVertexPartition.Bipartite(core_graph, 
                                                            resolution_parameter_01=r / len(numbered_edges), 
                                                            degree_as_node_size=True)

    diff = optimiser.optimise_partition_multiplex([p_01, p_0, p_1], layer_weights=[1, -1, -1], n_iterations=-1)
    memberships_by_resolution[r] = list(p_01.membership)
    print(r, len(set(memberships_by_resolution[r])))


#Connect names to module indexing

# the subgraph is reindexed from the larger graph
# build map from subgraph index (j) to full graph index (i)
core_node_to_node = {j:i for j, i in enumerate(i for i,g in enumerate(ccss.membership) if g == 0)}

# revese the dict we made to turn nodes into numbers
node_name_dict = {i:n for n,i in node_id_dict.items()}

# combine the dict to map subgaph numbers back to graph numbers and then to names
core_node_to_name = {j:node_name_dict[i] for j, i in core_node_to_node.items()}

# Make a map from node to subgraph for the n small subgraphs. 
# Let's call these clusters 1,2,... and then add n to all the module indexes. These n were set aside before the module fining code was run, 
# so they are the same for all param values (but may differ by annot col).
n = len(ccss.subgraphs()) - 1

# in the ccss membership, 0 is the big subgraph we found modules in and 1,2,.. are the other much smaller subgraphs
name_to_subgraph = {node_name_dict[i]:sgi for i, sgi in enumerate(ccss.membership)}

# Make a map from node name to module ID for each parameter value.
output_data['membership_by_resolution'] = {
    res:{core_node_to_name[i]:(mi + n + 1) for i, mi in enumerate(membs)}
    for res, membs in memberships_by_resolution.items()
}

# update dicts with unconnected subgraphs as modules 1,2...
for node_name, subgraph in name_to_subgraph.items():
    if subgraph > 0:
        for res in memberships_by_resolution.keys():
            output_data['membership_by_resolution'][res][node_name] = \
                subgraph

with open(str(output), 'wt') as json_data_out:
    json.dump(output_data, json_data_out)

