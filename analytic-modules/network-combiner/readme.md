# Network Combiner 

**Description:** 

The network combiner takes a set of networks as input, and delivers a single network by 
applying an operation on the input networks. It can either calculate the union or the 
intersection of the networks, or the set difference of the arguments.

At least one input network file must be provided. If all the input files are representing 
empty networks, then an output file will be also empty. The output network will contain no 
duplicate links. We are treating all the links in the input files as undirected.

When the 'union' method is chosen, then all those links will appear in the output networks,
which were present in any of the input networks. For the 'intersection' method, only the links
present in all input networks should be added to the output. If the 'difference' method was selected, 
then only those links will appear in the output, which were not present in all the input files. 
(So in this definition: difference = union - intersection)

If there is only a single network as an input, then the union and the intersection will be equal 
to the input file, while the difference is an empty network.

The metadata stored for the input files for the nodes or links will not be preserved in the 
output file.

It does not matter what method will be used, because every SNP metadata information will be copied to the output file
in a deduplicated way, separated by | delimiter.


**Parameters:**

--input-files <comma separated list of network files from each iteration> [mandatory]

--output-file <path to an output network set file> [mandatory]

--method <method name: union|intersection|difference> [mandatory]


**Exit codes:**

Exit code 1: One of the specified input file doesn't exists!
Exit code 2: The method doesn't exists!

