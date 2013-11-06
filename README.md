compute-accuracy
================

Compute the accuracy of different RNA alignment methods

All scripts and sample data files are located in tools/

The .gtf file for the drosophila genome used is in drosophila.gtf
Flux simulated strand data is located in simulation.bed
The tophat output is contained in accepted_hits.sam
Cufflinks output is in cufflinks.gtf

To compare genome coverage, execute:
  ./comp_coverage actual predicted
e.g. to compare the tophat alignments to the actual flux simulation data, run
  ./comp_coverage simulation.bed accepted_hits.sam

To compare junction locations, execute
  ./comp_junctions actual predicted
e.g. to compare the cufflinks results to the reference genome, run
  ./comp_junctions drosophila.gtf cufflinks.gtf

find_junctions.py uses a naive method of computing junctions by looking for local minima and maxima in the derivate of the coverage function. To compare the results of this method to the actual junctions, run
  ./find_junctions drosophila.gtf accepted_hits.sam

