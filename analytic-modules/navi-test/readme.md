# NaviTest

**Description:** 

The tool will be used in integration tests, simulating the different behaviours of the analytical tools. It is basically repeating in the output file what it reads from the input file. Before creating the output it is optionally waiting for the specified time. After producing the output file, it will exit with the exit code given as a parameter. It will create the output file even if the exit code is non-zero. 

```python3 navi_test.py -i inputfile.csv -o outfile.csv -w 3 -e 0```

**Parameters:**

--input -i <path to an existing file> [mandatory]

--output -o <path to the new output file> [mandatory]

--waitsec -w <0..255> [optional, default is 0]

--exitcode -e <0..255> [optional, default is 0]


** execute tests: **

```
cd analytic-modules/navi-test/test
PYTHONPATH=$PYTHONPATH:`pwd`/.. pytest
```