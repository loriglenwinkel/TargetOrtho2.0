# TargetOrtho2.0

#this version is not ready for deployment as of May 4, 2018. Please check back soon
system Requirements

system Requirements

macOS High Sierra version 10.13.3 (not tested on other versions)

Command Line TOols (\macOS 10.13) for Xcode 9.3
Install here: https://developer.apple.com/download/more/

meme version 4.12.0 (fimo tool scans genomes for motif matches)
Installation guide: http://meme-suite.org/doc/install.html?man_type=web
Download here: http://meme-suite.org/meme-software/4.12.0/meme_4.12.0.tar.gz

Make sure fimo is executable by copying to /usr/local/bin. Open the terminal and enter this command with "user_name" replaced with your personal user name.
command: sudo cp /Users/"user_name"/meme/bin/fimo /usr/local/bin/

BEDOPS 2.4.30 (This tool associates motif match coordinates with adjacent gene features.
Installation guide: http://bedops.readthedocs.io/en/latest/content/installation.html#mac-os-x

INSTALL:
Download and unzip TargetOrtho2.0_download directory	by double clicking.
Open the directory and double click the file named "setup.command". Wait until you see the "Process completed" notification before running TargetOrtho.


Running TargetOrtho 2.0 Using the GUI application (user interface that runs the command line tool automatically):
1. find TargetOrtho_App.command file within the TargetOrtho2.0 directory.
2. Double click TargetOrtho_App.command and follow the promt.

Running TargetOrtho from the command line:
1. navigate to the TargetOrtho2.0 directory from the command line
2. command: python taragetortho.py -h to see options.

Running TargetOrtho from the commane line with an example motif file:
This command runs the TargetOrtho with the COE motif using a p value threshold for fimo of 0.0001 and resticts the motif match search to 500 nucleotides within a gene start position.
command: python targetortho.py -f data/input_motif_examples/COE_motif_PSSM_meme4.txt -d 500 -p 0.0001
