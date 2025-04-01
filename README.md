Custom code associated with "Bacterial peptide deformylase inhibitors induce prophages in competitors", written by Julie Chen.<br>
To cite the code, please use the most recent copy of the manuscript: TBA

Related kChip image analysis code: https://github.com/megantse/kchip  
<hr>

**Included:**
- **kChip dataset analysis**: combination of python scripts & jupyter notebooks <br>
- **microtitre plate hit calling analysis**: jupyter notebook<br>
  *individual plate analysis not provided as it is unique to BioTek Cytation 5, available upon request*
- **calling prophage induction**: python script 
<hr>

**kCHIP DATA ANALYSIS**<br>
Analysis for lysis time courses from pairwise co-cultures.<br><br>
*pre-requisite: image analysis completed (see above for code). The code provided is only for Steps 2-5.*<br>
*note: scripts and notebooks rely on source scripts found in "sytox_scripts"*<br><br>
**Step 1:** process raw images to deconvolute microwells to pairwise co-cultures or monocultures (e.g. a pair of media + monoculture)<br>
**Step 2:** run QC on the dataset & blank media values from all co-cultures<br>
**Step 3:** aggregate replicates into summary data (e.g. mean)<br>
*note: this can be a very slow step from the bootstrapping*<br>
**Step 4:** compile co-culture values with respective monocultures & concatenate all kChips into one batch dataset<br>
**Step 5:** call hits on the entire batch<br>
<hr>

**PLATE DATA ANALYSIS**<br>
Analysis for lysis time courses from pairwise co-cultures.<br><br>
Jupyter notebook analogous to Step 5 from kChip data analysis to call hits (e.g. stats on the null population, thresholds)
<hr>

**AUTOMATED PROPHAGE INDUCTION DETECTION** <br>
Calls putatively induced prophages from mapped Illumina reads.<br>
*pre-requisite: reads have been mapped and basepair-resolved depth as been calculated<br><br>*
For simplicity & efficiency, the associated work was mapped with Bowtie2, but feel free to use bbmap on shorter references to not split reads between matching sites (which ultimately affects coverage, the code is defaulted to searching for at least 10X coverage for a putative prophage).<br><br>
After mapping, I recommend using Samtools to convert from .sam to .bam (view -bS -o), sorting (sort), indexing (index), and calculating depth (depth -a -H and produce a .txt).
optional: de-duplicating (rmdup) to remove PCR artefacts<br>

=========================================<br>
Once you have the depth calculated, run the script as:<br>
<code>/path_to_script_from_current_directory/phage_call.py --date "YYYYMMDD" --source_dir "/path_to_txt/" --save_dir "/path/" --clr "optional_HEX"</code>

This script will produce .csv files with information on putatively induced regions (in case you want to review manually) and called phages, along with background coverage and bp-resolution which is used to produce a coverage plot for each called phage. Ensure you have the source script too when running phage_call.py. <br>


