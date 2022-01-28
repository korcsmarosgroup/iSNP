# Network collector 

**Description:** 

The collector modules are the output points of iterations in the workflow. These
modules can be used to propagate the results from the iterations to the rest of
the workflow.

This collector takes a single network file as an input from each iteration, 
and producing a single network set file as an output.

The input format is a standard NavigOmiX MITAB file. The output is a NavigOmiX 
network set file, a tar.gz archive with all the network files collected in 
each iteration.

If the input file list is empty, then an output tar.gz file will be created with no content.

**Parameters:**

--input-files <comma separated list of network files from each iteration> [mandatory]

--output-file <path to the tar.gz network set file> [optional]   
