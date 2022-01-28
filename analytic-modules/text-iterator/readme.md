# Text iterator

**Description:** 

The iterator modules are the entry points of iterations in the workflow. This 
iterator takes a NavigOmiX primitive (text file) as an input, and then 
iterate over each line of the input file, providing single NavigOmiX primitive 
variable for each iteration, represented as a separate text file.

The output is an ordered set of NavigOmiX primitive files in a given output folder. 
The file names will be `iteration-1`, `iteration-2`, ... The ordering of the 
output files is deterministic, following the order of lines in the input file.

**Parameters:**

--input-file <path to the input file> [mandatory]   

--output-folder [mandatory]