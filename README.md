# RNAfbinv 2.0

RNAfbinv is a fragment based RNA design tool. It uses a simulated annealing process to optimize a 2D RNA structure.
The similarity is based on fragment based design. A tree alignment is done based on nodes (structural motifs).
Nodes are comparable if they are both bounded motifs (stems) or unbounded motifs (multi loop, interior loops, bulges ect...).
Each iteration the target motif tree will be aligned to the current candidate tree.
The best alignment with the addition of other valuable features will generate a design score.
Design score of 0 is exact fit but even higher scores can generate a good candidate.

## Attaching Vienna RNA

Vienna RNA package is required for RNAfbinv to work. This must be installed separately.
Current version was tested with Vienna 2.4 and above.
RNAfbinv will identify Vienna package if it's bin directory is in PATH or if the location is given via parameter.

## Usage

The design process can be ran using the following code:
```python
import rnafbinv
rnafbinv.main(command_line_arguments)
```

### Command line arguments:

-h : Shows usage text
-i <number of iterations> : sets the number of simulated annealing iterations (default is 100)
--seed <random number generator seed> : a long number that is used by the random number generator
-t <number of look ahead attempts> : number of look head mutation attempts for each iteration (default is 4)
-e : designs a circular RNA (default is False)
-m <motif[,...]> : comma separated list of motifs to preserve,
-s <starting sequence> : the initial sequence for the simulated annealing process
-r : force starting simulated annealing with a random sequence
-p <MFE|centroid> : uses RNAfold centroid or MFE folding. (default is MFE)
--verbose : Additional info message on simulation process
--debug : Debug information
-l <log file path> : Logging information will be written to a given file path (rewrites file if exists)
--length <length diff> : The resulting sequence size is target structure length +- length diff (default it 0)

-f <input file path> : Path of ini file that includes mandatory information. Some options can also be set via file.
                       command line options take precedence.

### File input:

```
# mandatory
TARGET_STRUCTURE=<target structure>
TARGET_SEQUENCE=<target sequence>
# optional
TARGET_ENERGY=<target energy>
TARGET_MR=<target mutational robustness>
SEED=<random seed>
STARTING_SEQUENCE=<starting sequence>
ITERATION=<number of simulated annealing iterations>
```

## GUI / Command line

You can download the RNAfbinv wrapper from <MISSING LINK>
It includes python code to run the GUI / command line and the [VARNA](http://varna.lri.fr/ "VARNA rna homepage") java package to generate 2D images.
If you remove the VARNA jar or do not have java installed, images will not be generated but the design process will proceed normally.

To run command line just use the '-cl' flag for RNAfbinv.py. by default a GUI version will run.

## Webserver
