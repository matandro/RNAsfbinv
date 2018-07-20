# RNAfbinv 2.0

RNAfbinv is a fragment based RNA design tool. It uses a simulated annealing process to optimize a 2D RNA structure.<br/>
The similarity is based on fragment based design. A tree alignment is done based on nodes (structural motifs).<br/>
Nodes are comparable if they are both bounded motifs (stems) or unbounded motifs (multi loop, interior loops, bulges ect...).<br/>
Each iteration the target motif tree will be aligned to the current candidate tree.<br/>
The best alignment with the addition of other valuable features will generate a design score.<br/>
Design score of 0 is exact fit but even higher scores can generate a good candidate.

## Attaching Vienna RNA

[Vienna RNA package](https://www.tbi.univie.ac.at/RNA/ "Vienna RNA home") is required for RNAfbinv to work. This must be installed separately.<br/>
Current version was tested with Vienna 2.4 and above. RNAfbinv will identify Vienna package if it's bin directory is in PATH.<br/>
If you wish to link a specific installation of Vienna set the VIENNA_PATH environment variable to the correct bin directory.

You can set Vienna location in python
```python
import os
os.environ['VIENNA_PATH'] = "VIENNA_BIN_DIR_PATH"
```

or directly via the vienna script
```python
from rnafbinv import vienna
vienna.set_vienna_path("VIENNA_BIN_DIR_PATH")
```

## Usage

The design process can be ran using the following code:
```python
from rnafbinv import RNAfbinvCL
RNAfbinvCL.main(command_line_arguments)
```

To generate a tree for a specific sequence / structure:<br/>
Structure is a dot bracket notation structure and sequence is an IUPAC string with the same length
```python
from rnafbinv import shapiro_tree_aligner
shapiro_tree_aligner.get_tree(sructure, sequence)
```

To compare two trees and score them:
alignment_rules has a default value and is optional
```python
from rnafbinv import shapiro_tree_aligner
shapiro_tree_aligner.align_trees(source_tree, tree_target, alignment_rules)
```

### Command line arguments:

```
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
--cl : Uses command line version (Only relevant for wrapper)
```

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

You can download the RNAfbinv wrapper from [RNAfbinv2.0 user application]()<br/>
It includes python code to run the GUI / command line and the [VARNA](http://varna.lri.fr/ "VARNA rna homepage") java package to generate 2D images.<br/>
If you remove the VARNA jar or do not have java installed, images will not be generated but the design process will proceed normally.<br/><br/>

To run command line just use the '-cl' flag for RNAfbinv.py. by default a GUI version will run.

To specify Vienna package location please update the 'VIENNA' parameter in config.ini (or set VIENNA_PATH environment variable)<br/>
To specify Java location please update the 'JAVA' parameter in config.ini (or set JAVA_PATH environment variable)<br/>
In both cases you need to specify the folder where the binary file sits. Only needed if not specified in PATH.

## Webserver

RNAfbinv2.0 can be found in a web server combined with incaRNAtion. The webserver generates starting seeds using incaRNAtion global sampling algorithm.<br/>
Te seed sequences are then sent to RNAfbinv2.0 for design. [incaRNAfbinv web server](https://www.cs.bgu.ac.il/incaRNAfbinv/ "incaRNAtion & RNAfbinv")