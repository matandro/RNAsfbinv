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

## GUI / Command line

You can download the RNAfbinv wrapper from [RNAfbinv2.0 git repository](https://github.com/matandro/RNAsfbinv/)<br/>
The main folder includes python code to run the GUI / command line and a configuration file:
* RNAfbinv.py - A GUI wrapper for RNAfbinv2.0
* RNAfbinvCL.py - A command line wrapper for RNAfbinv2.0
* **Required** varna_generator.py - Used to generate images based on [VARNA](http://varna.lri.fr/ "VARNA rna homepage")
* **Required** config.ini - Configuration file with paths to required software (information below).

If you remove the VARNA jar or do not have java installed, images will not be generated but the design process will proceed normally.<br/><br/>

To specify [vienna package](https://www.tbi.univie.ac.at/RNA/ "The ViennaRNA Package homepage") binary folder please update the 'VIENNA' parameter in config.ini (or set VIENNA_PATH environment variable)<br/>
To specify Java binary folder please update the 'JAVA' parameter in config.ini (or set JAVA_PATH environment variable)<br/>
To specify [VARNA](http://varna.lri.fr/ "VARNA rna homepage")'s jar file please update the 'VARNA' parameter in config.ini (or set VARNA_PATH environment variable)<br/>
Note that if the java or vienna package binaries are in your environment variables you may leave it empty.

Example to a valid config.ini file which has java installed and within the system's path:
```
[PATH]
VIENNA=~/ViennaRNA/bin/
#JAVA=
VARNA=~/VARNA/VARNAv3-93.jar
```

### Command line arguments:

```
usage: RNAfbinvCL.py [-h] [-l LOG_OUTPUT] [--verbose | --debug]
                     [-p {MFE,centroid}] [-i ITERATIONS] [--seed SEED]
                     [-t LOOK_AHEAD] [--reduced_bi REDUCED_BI] [-e]
                     [--seq_motif] [-m MOTIF_LIST] [-s STARTING_SEQUENCE | -r]
                     [--length LENGTH] [-f INPUT_FILE]

optional arguments:
  -h, --help            show this help message and exit
  -l LOG_OUTPUT, --log_output LOG_OUTPUT
                        Path to output log file. (default: None)
  --verbose             Increase output verbosity. (default: False)
  --debug               Debug level logging. (default: False)
  -p {MFE,centroid}, --structure_type {MFE,centroid}
                        uses RNAfold centroid or MFE folding. (default: MFE)
  -i ITERATIONS, --iterations ITERATIONS
                        Sets the number of simulated annealing iterations.
                        (default: 100)
  --seed SEED           Random seed used in the random number generator.
                        (default: None)
  -t LOOK_AHEAD, --look_ahead LOOK_AHEAD
                        Number of look head mutation attempts for each
                        iteration. (default: 4)
  --reduced_bi REDUCED_BI
                        Remove extra penalty for removal or addition of bulges
                        and interior loops under the given size. Alignment
                        penalties still occur. (default: 0)
  -e, --circular        Designs a circular RNA. (default: False)
  --seq_motif           Enables increased penalty for insertion or deletions
                        within marked regions (lower case characters in
                        sequence constraint). The feature was added to control
                        multi base sequence constraints (sequence motifs).
                        Only valid within a specific structural motif.
                        (default: False)
  -m MOTIF_LIST, --motif_list MOTIF_LIST
                        A comma separated list of motifs that are targeted for
                        preservation with size.Single motif format: <motif
                        No>[M|H|E|I|S|B]<motif No of bases>. Use
                        rnafbinv.ListMotifs.list_motifs(structure) to retrieve
                        a list of legal motifs for a given structure.
                        (default: [])
  -s STARTING_SEQUENCE, --starting_sequence STARTING_SEQUENCE
                        The initial sequence for the simulated annealing
                        process in IUPAC nucleotide codes. (default: None)
  -r, --random_start    Start simulated annealing with a random sequence.
                        (default: False)
  --length LENGTH       Maximum variation in result length compared to target
                        structure. (default: 0)
  -f INPUT_FILE         Path of ini file that includes mandatory information.
                        Some options can also be set via file. command line
                        options take precedence. (default: None)
```

### Input file format (the '-f' parameter):

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

## Webserver

RNAfbinv2.0 can be found in a web server combined with incaRNAtion. The webserver generates starting seeds using incaRNAtion global sampling algorithm.<br/>
Te seed sequences are then sent to RNAfbinv2.0 for design. [incaRNAfbinv web server](https://www.cs.bgu.ac.il/incaRNAfbinv/ "incaRNAtion & RNAfbinv")