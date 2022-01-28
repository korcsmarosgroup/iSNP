# Enrich network with additional links / neighbours

**Description:** 

This tool will take a single network file as input, and enrich it with additional links and
neighbours using some reference protein interaction database. Both the input network and 
the reference database is expected to be given in standard NavigOmiX MITAB format.

The tool will not map any protein or gene IDs from the input network or from the reference
interaction database. The user must make sure these identifiers are aligned before starting 
the tool.

The output of the tool is a standard MITAB file (using standard NavigOmiX identifiers). The 
output will not contain self-loops or duplicate links. All the metadata from the input or
from the reference database will be copied to the output file as well. If an interaction is 
present both in the input network and in the reference database, it will be copied from the
input network, as it is (with all its metadata, without any change).

The tool will take both the input network and both the reference interaction database as
undirected. 

The tool will take any additional information (including alternative tags) within the enriched
MITAB formats. Therefore, reference networks should be complete and loaded through the `MiTabHandler`
if mapping beyond interactor-a and interactor-b.

If the input mitab file is empty, the tool will always return with an empty file. If the
reference interaction database is empty, or if there is no overlap between the input network
and the reference interaction database, then the output will be same as the input file.
 
    
**Parameters:** 

-i, --inuput <mitab path>  : input MITAB file path [mandatory]

-d, --distance <0..255>    : number of hops for the enrichment [optional, default: 1]
                             0: no new node, just adding links between existing nodes
                             1: adding first neighbour nodes, then adding links
                             2: adding first and second neighbour nodes, then adding links 
                             ...

-r, --reference-net <path> : path to the reference interaction database MITAB file [mandatory]

-o, --output <path>        : output MITAB file [mandatory]
