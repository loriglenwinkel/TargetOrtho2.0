TargetOrtho: A Phylogenetic Footprinting Tool to Identify Transcription Factor Targets.

system Requirements

macOS High Sierra version 10.13.3 (not tested on other versions)

Command Line Tools (\macOS 10.13) for Xcode 9.3
Install here: https://developer.apple.com/download/more/

meme version 4.12.0 (fimo tool scans genomes for motif matches) 
Installation guide: http://meme-suite.org/doc/install.html?man_type=web
Download here: http://meme-suite.org/meme-software/4.12.0/meme_4.12.0.tar.gz

Make sure fimo is executable by copying to /usr/local/bin/. Open the terminal and enter this command with "user_name" replaced with your personal user name.
command: sudo cp /Users/"user_name"/meme/bin/fimo /usr/local/bin/

BEDOPS 2.4.30 (This tool associates motif match coordinates with adjacent gene features.
Installation guide: http://bedops.readthedocs.io/en/latest/content/installation.html#mac-os-x

INSTALL:
Download and unzip TargetOrtho2.0_download directory by double clicking.
Open the directory and double click the file named "setup.command". Wait until you see the "Process completed" notification before running TargetOrtho.


Running TargetOrtho 2.0 Using the GUI application (user interface that runs the command line tool automatically):
1. Double click TargetOrtho_App.command and follow the promt.
2. when the job is done running the output directory window will pop up.

Running TargetOrtho from the command line:
1. navigate to the TargetOrtho2.0 directory from the command line
2. command: python taragetortho.py -h to see options.

Running TargetOrtho from the commane line with an example motif file:

command: python targetortho.py -f data/input_motif_examples/COE_motif_PSSM_meme4.txt -d 500 -p 0.0001

This command runs TargetOrtho with the COE motif using a p value threshold for fimo of 0.0001 and resticts the motif match search to 500 nucleotides within a gene start position.

#changes to version 2.0

*the entire job runs in less than 10 to 45 minutes (instead of days!).

*command line version or user friendly double clickable application available of OS X.

*all genomes are updated to the latest version available from wormbase parasite as of May 2018.

*Additional nematode genomes are included (P. pacificus, P. exspectatus, and A. lumbricoides)

*the main result file of interest is the "TargetOrtho2_ranked_Genes_summary.csv" file.

*Results are ranked per reference species gene (instead of per binding site in the reference species genome).

*Results are expressed as prediction probablies where values between 0 and 1 are increasingly likely to be real transcription factor targets (0.99 =99% probabiliy of being a real TF target gene) and results between -1 and 0 are unlikely to be real targets (-0.99=99% probabiliy that the gene is not a real TF target gene).

*Predictions are generated from modeling of validated transcription factor target genes from UNC-3, ASE, and ttx-3/ceh-10 target genes using a support vector machine to train the model.

*detailed information about individual motif match results per species are found in separate spreadsheets

*results from the fimo motif scanner are stored for browsing in the results 

# TargetOrtho2.0
