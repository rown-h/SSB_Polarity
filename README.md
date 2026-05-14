# From PRE–NMR to coloured PyMOL structures!
This repository contains scripts that achieve the following:

* Produce CSVs and column graphs of PRE signal intensity ratios
* Exports the data to a .txt format readable by PyMOL

## Acknowledgements
Thank you to those who developed helper PyMOL scripts, included here in the `Scripts/PyMOL_Scripts` folder:

* Campbell, R.L., Holder, T., Asai, S., 2013. data2bfactor.
* Holder, T., 2010. spectrumany.


## Software requirements
| Software             | Version tested                | Citation                                                                 |
|----------------------|-------------------------------|--------------------------------------------------------------------------|
| RStudio              | 2026.01.1 Build 403 for macOS | [Posit Software, PBC (2026)](https://posit.co/download/rstudio-desktop/) |
| R                    | 4.3.3 for macOS               | [R Core Team (2024)](https://www.R-project.org/)                         |
| R package: here      | 1.0.1                         | [Müller K (2020)](https://CRAN.R-project.org/package=here)               |
| R package: tidyverse | 2.0.0                         | [Wickham et al. (2019)](https://doi.org/10.21105/joss.01686)             |
| Python               | 3.13.9                        | [Python Software Foundation (2026)](https://www.python.org/)             |
| PyMOL                | 3.1.5.1                       | [Schrödinger, LLC (2026)](https://www.pymol.org/)                     |


## To use:

1. Export .list files from [POKY](https://sites.google.com/view/pokynmr/home) (or its predecessor, [SPARKY](https://www.cgl.ucsf.edu/home/sparky/))
  * Ensure peaks displayed in the form "A87N-H"
  * Only peak, coordinates (N and H) and data height columns selected
  * Save to `Data/List_Files/YOUR_PROTEIN` with the following format:
  * DNA_ox.list (oxidised, active PRE label) or DNA_red.list (reduced, inactive PRE label)

 \
2. Save your protein name and amino acid sequence to `Data/protein_sequences_position_shifts`
  * Optionally include a position shift value here, i.e., the difference between the numbering of the amino acid sequence provided and the residues from POKY.
    This option may be useful if your protein has a non-native N-terminus that is included in the NMR spectrum, but you do not want in your output.

 \
3. Create a folder for each protein in `PyMOL_Structures`

 \
4. Populate the protein-DNA complex 3D structures.
  * Return to the `PyMOL_Structures/YOUR_PROTEIN` folder, and make the PyMOL structure .pse file.
      * You can use `fetch <PDB ID>` in [PyMOL](https://www.pymol.org/)
  * Isolate the protein to a selection called 'prot'
  * Isolate the DNA to a selection called 'DNA'
      * You can do this by toggling on the sequence, highlighting all the relevant residues there, and typing `select DNA, sele`

 \
5. **[OPTIONAL]** Populate `Data/Distance_Data` with the distances between the DNA ends and amino acid backbones.
   * This can be used to graph a distance curve, which should broadly correspond with the signal intensity ratio if a protein binds unidirectionally.
   * Open the PyMOL structure and load the provided script in PyMOL: `run /path/to/SSB_Polarity/Scripts/PyMOL_Scripts/dna_distances.py`
   * Run the custom function: `dna_distances(out_dir = '/path/to/SSB_Polarity/Data/Distance_Data', protein = 'my_protein')`
      * Adjust other parameters if necessary for your PyMOL file.

 \
6. Run the R script, `Scripts/PRE_ratio.R`
   * This has toggles for all the outputs you may want to make.

 \
7. Colour the PyMOL structure
   * Open the PyMOL structure .pse file.
   * Copy all the text from its corresponding colouring.txt file, and paste in the PyMOL command line.
       * The structure will be coloured!

## Example output
The default colour scale is:
| Colour       | Meaning           | Detail                                                                                                     |
|--------------|-------------------|------------------------------------------------------------------------------------------------------------|
| White        | Unassigned        | Peak is not annotated in either the _ox or _red .list file                                                 |
| Grey         | Unaffected        | Signal intensity ratio ≥ mean                                                                              |
| Blue         | Affected          | Signal intensity < mean                                                                                    |
| Green        | Strongly affected | Signal intensity < mean - standard deviation                                                               |
| Green circle | Peak disappeared  | Signal only visible in the _red.list file, completely disappeared after complexing with oxidised substrate |

![Example of the column graphs, PyMOL structure and colour scale produced by this workflow](Images/example_figures.svg)
