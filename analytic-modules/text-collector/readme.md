# Text collector 

**Description:** 

The collector modules are the output points of iterations in the workflow. These
modules can be used to propagate the results from the iterations to the rest of
the workflow.

This collector takes NavigOmiX primitive files (text files) as input from 
each iteration, and combining them into a single primitive output variable.

The tool will make sure that the input primitives from each iteration will be 
started in a new line (adding extra new line characters if needed).

If the input file list is empty, then an empty output file will be created.

**Parameters:**

--input-files <comma separated list of primitive files from each iteration> [mandatory]

--output-file <path to the output primitive file> [mandatory]   
