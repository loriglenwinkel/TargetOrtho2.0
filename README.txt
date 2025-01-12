Here's the spell-checked and corrected version:

TargetOrtho: A Phylogenetic Footprinting Tool to Identify Transcription Factor Targets.

System Requirements
macOS X High Sierra version 10.13.3 
python2.7 (requires sklearn and pandas modules)
Xcode 9.3
Xcode Command Line Tools 9.3
MEME version 4.12.0 (see MEME installation notes below)
BEDOPS 2.4.30 
bedtools version v2.27.1

macOS X El Capitan version 10.11.6
python2.7 (requires sklearn and pandas modules)
Xcode 7.3
Xcode Command Line Tools 7.3
MEME version 4.12.0 (see MEME installation notes below)
BEDOPS version 2.4.30
bedtools version v2.27.1

Ubuntu 14.04 (command line version only)
python2.7 (requires scikit-learn version > 0.18.1 and pandas modules - can all be installed with anaconda2)
MEME version 4.12.0 (make sure fimo is executable by copying fimo to /usr/local/bin/) 
BEDOPS version 2.4.35
bedtools version v2.17.0

*MEME installation notes:
Use MacPorts to install dependencies (not Homebrew). Follow MEME install instructions and use the ./configure --prefix=$HOME/meme option.
Make sure fimo is executable by copying to /usr/local/bin/. Open the Terminal App and type: sudo cp $HOME/meme/bin/fimo /usr/local/bin/

INSTALL OS X:
Download and unzip TargetOrtho2.0 directory by double clicking. Open the TargetOrtho2.0 directory and double click the file named "setup.command". Wait until you see the "Process completed" notification before running TargetOrtho.

Running TargetOrtho 2.0 Using the GUI application (user interface that runs the command line tool automatically):
Double click TargetOrtho_App.command and follow the prompt.
When the job is done running, the output directory window will pop up.

Running TargetOrtho from the command line:
Navigate to the TargetOrtho2.0 directory from the command line
command: python targetortho.py -h to see options.

Running TargetOrtho from the command line with an example motif file:
command: python targetortho.py -f data/input_motif_examples/COE_motif_PSSM_meme4.txt -d 500 -p 0.0001
This command runs TargetOrtho with the COE motif using a p-value threshold for fimo of 0.0001 and restricts the motif match search to 500 nucleotides within a gene start position.

INSTALL Linux (Ubuntu):
git clone https://github.com/loriglenwinkel/TargetOrtho2.0.git
tar -xvzf TargetOrtho2.0.tar.gz
cd TargetOrtho2.0
python setup.command
command line: python targetortho.py -options

Changes in version 2.0:
* The entire job runs in less than 10 to 45 minutes (instead of days!).
* Command line version or user-friendly double-clickable application available for OS X.
* All genomes are updated to the latest version available from WormBase ParaSite as of May 2018.
* Additional nematode genomes are included (P. pacificus, P. exspectatus, and A. lumbricoides)
* The main result file of interest is the "TargetOrtho2_ranked_Genes_summary.csv" file.
* Results are ranked per reference species gene (instead of per binding site in the reference species genome).
* Results are expressed as prediction probabilities where values between 0 and 1 are increasingly likely to be real transcription factor targets (0.99 = 99% probability of being a real TF target gene) and results between -1 and 0 are unlikely to be real targets (-0.99 = 99% probability that the gene is not a real TF target gene).
* Predictions are generated from modeling of validated transcription factor target genes from UNC-3, ASE, and ttx-3/ceh-10 target genes using a support vector machine to train the model.
* Detailed information about individual motif match results per species are found in separate spreadsheets
* Results from the fimo motif scanner are stored for browsing in the results
