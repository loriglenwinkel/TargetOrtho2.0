
# coding: utf-8


import subprocess
from threading import Thread
import pandas as pd
from functools import reduce
import os
import numpy as np
import scipy.stats
from sklearn import linear_model
from sklearn import svm
import datetime as dt
from optparse import OptionParser
import sys
#supress pandas errant copy warnings
pd.options.mode.chained_assignment = None

#requires meme_4.12.0 command line version (uses fimo). Copied fimo after install from cp /Users/lori/meme/bin/fimo ../scripts/fimo
#os 10.13.3, Xcode command line tools 9.3
#bedops 2.4.30
#citation: http://bioinformatics.oxfordjournals.org/content/28/14/1919.abstract
#          https://doi.org/10.1093/bioinformatics/bts277
#version:  2.4.30 (typical)
#authors:  Shane Neph & Scott Kuehn, bedtools v2.27.1

#associate genes script
#takes fimo file as input
#references genes.bed files and exons.bed file

#need fimo to be executable
#cp /Users/lori/meme/bin/fimo /usr/local/bin/


parser = OptionParser()
parser.add_option( '-f', '--file', dest='infile', help='-f filepath/motif_input.txt' )
parser.add_option('-o', '--output_dir',dest='output_dir',help='-o my_results_dir, default output is to www.alice.cumc.columabie.du/TargetOrtho/jobID. A link to this webpage is printed at the end of the job')
parser.add_option('-j','--jobID',dest='jobID',help='-g 123')
parser.add_option( '-d', '--max_dist', dest='max_dist', help='-x None #maximum distance to search for motif match from start codon')
parser.add_option( '-y', '--ordered', dest='motif_order', help='True or False: -o True means motif matches must be same order as input file on DNA')
parser.add_option('-p','--p', dest='p_cutoff',help='--p_cutoff 1e-6')
parser.add_option('-s','--species_file',dest='species_file',help='-s  species.txt -text file with list of species to include in motif searches')
parser.add_option('-r','--ref_species',dest='ref_species',help='-r C.elegans or -r P.pacificus')
( options, args ) = parser.parse_args()

print('sys.argv',sys.argv)


#TargetOrtho_path=os.path.abspath(".")
TargetOrtho_path=os.path.dirname(os.path.abspath(__file__))

now = dt.datetime.now()
if options.jobID==None:
    jobID="j%s%s%s%s%s%s" %(now.year,now.month,now.day,now.hour,now.minute,now.second)        
else:
    jobID=options.jobID

max_search_distance=options.max_dist

if options.p_cutoff==None:
    p_value="0.0001"
else: 
    p_value=options.p_cutoff

    
motif_file_path=options.infile


if options.output_dir==None:
    output_dir="%s/%s_TargetOrtho2.0_Results" %(TargetOrtho_path,jobID)
else:
    output_dir=options.output_dir

    
if options.species_file==None:
    speciesList=["c_eleg","c_brig","c_bren","c_rema","c_japo"]
else:
    speciesList=list(pd.read_csv(options.species_file,sep='\t',header=None)[0].values)
    
if options.ref_species=="C.elegans":
    speciesList=["c_eleg","c_brig","c_bren","c_rema","c_japo","p_paci","p_exsp","a_lumb"]
if options.ref_species=="P.pacificus":
    speciesList=["p_paci","c_eleg","c_brig","c_bren","c_rema","c_japo","p_exsp","a_lumb"]
 
    
params=[TargetOrtho_path,jobID,max_search_distance,p_value,motif_file_path,output_dir,speciesList]
for p in params:
    print(p)


genomesDic={"c_eleg":"caenorhabditis_elegans.PRJNA13758.WBPS10.genomic.fa",
            "c_brig":"caenorhabditis_briggsae.PRJNA10731.WBPS10.genomic.fa",
           "c_bren": "caenorhabditis_brenneri.PRJNA20035.WBPS10.genomic.fa",
           "c_japo": "caenorhabditis_japonica.PRJNA12591.WBPS10.genomic.fa",
           "c_rema": "caenorhabditis_remanei.PRJNA53967.WBPS10.genomic.fa",
           "p_paci": "pristionchus_pacificus.PRJNA12644.WBPS10.genomic.fa",
            "p_exsp": "pristionchus_exspectatus.PRJEB6009.WBPS10.genomic.fa",
           "a_lumb": "ascaris_lumbricoides.PRJEB4950.WBPS10.genomic.fa"}
versionDic={"c_eleg":"WS258",
            "c_brig":"WS258",
           "c_bren": "WS258",
           "c_japo": "WS258",
           "c_rema": "WS258",
           "p_paci": "WS263",
            "p_exsp": "2014-12-WormBase",
           "a_lumb": "2014-06-50HGPpatch"}

def system_call(command):
    p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
    return p.stdout.read()

def mkDirs(jobID,output_dir):
    command="mkdir /'%s/'" %(output_dir)
    system_call(command)
    command="mkdir \'%s/temp_files/%s\'" %(TargetOrtho_path,jobID)
    system_call(command)
    command="mkdir \'%s/fimo_out\'" %(output_dir)
    system_call(command)
    command="mkdir \'%s/motif_match_data_per_species\'" %(output_dir)
    system_call(command)

def clear_temp():
    command="rm -r \'%s/temp_files/%s\'" %(TargetOrtho_path,jobID)
    system_call(command)
    for species in speciesList:
        command="rm  \'%s/fimo_out/%s_%s.bed\'" %(output_dir,jobID,species)
        system_call(command)
    for species in speciesList:
        command="rm  \'%s/fimo_out/%s_%s.bed.sorted\'" %(output_dir,jobID,species)
        system_call(command)
       
    
def clear_error():
    command="rm -r \'%s/temp_files/%s\'" %(TargetOrtho_path,jobID)
    print(command)
    system_call(command)
    command="rm -r \'%s\'" %output_dir
    print(command)
    system_call(command)
    
def writeParams():
    outfile=open("%s/%s_input_summary.txt" %(output_dir,jobID),'w')
    outfile.write("job ID: %s\n" %(jobID))
    outfile.write("motif file path: %s\n" %motif_file_path )
    
    #all params:
    outfile.write("species: ")
    for species in speciesList:
        outfile.write("%s, " %species)        
    outfile.write("\n")
    
    outfile.write("fimo p value threshold: %s\n" %p_value)
    outfile.write("max motif match search distance from gene start: %s\n" %max_search_distance)
    
    #full TargetOrtho command used
    outfile.write("command: ")
    for a in sys.argv:
        outfile.write("%s " %a)
    outfile.write("\n")
    

def runFIMO(species,jobID):
    command="fimo --oc \'%s/fimo_out/%s_%s\' --max-strand --thresh %s \'%s\' \'%s/genomes/%s\'" %(output_dir,jobID,species,p_value,motif_file_path,TargetOrtho_path,genomesDic[species])
    print(command)
    print("\n")
    system_call(command)

def runFIMO_threaded(speciesList,jobID):
    #run fimo in parallel for each species (results written to fimo_out/jobIDspecies/ dir)
    queue=[]
    for species in speciesList:
        t = Thread(target=runFIMO, args=(species,jobID,))
        queue.append(t)
        t.start()
    for q in queue:
        q.join()
        
def convertFimo(species,bed_path):
    fimo=pd.read_csv("%s/fimo_out/%s_%s/fimo.txt" %(output_dir,jobID,species),sep='\t')

    #convert fimo output to bed format
    fimo = fimo[['sequence_name', 'start', 'stop','# motif_id','score','strand']]
  
    fimo.loc[:,"# motif_id"]=list(fimo.index.values)
    fimo.to_csv(bed_path,header=False, index=False, index_label=None,sep='\t')
    
    #sort the bed file for bedops
    command="bedtools sort -i \'%s\' > \'%s.sorted\'" %(bed_path,bed_path)
    system_call(command)

def removeNA(species):
    #results with no left or right feature assocation have fewer columns in output file so can't be imported by pandas.read_csv
    c1="grep -v \'\tNA\'  \'%s/temp_files/%s/%s_associated_genes.bed\' >\'%s/temp_files/%s/%s_associated_genes_no_NA.bed\'" %(TargetOrtho_path,jobID,species,TargetOrtho_path,jobID,species)
    
    c2="grep \'\tNA\' \'%s/temp_files/%s/%s_associated_genes.bed\' >\'%s/temp_files/%s/%s_associated_genes_NA_all.bed\'"%(TargetOrtho_path,jobID,species,TargetOrtho_path,jobID,species)
    #get rid of line with no gene association data (4 NAs)
    c2b="grep -v \'\tNA\tNA\tNA\tNA\' \'%s/temp_files/%s/%s_associated_genes_NA_all.bed\' > \'%s/temp_files/%s/%s_associated_genes_NA.bed\'"  %(TargetOrtho_path,jobID,species,TargetOrtho_path,jobID,species)
    
    c1n="grep -v \'\tNA\' \'%s/temp_files/%s/%s_associated_genes.bed--no-overlaps\' >\'%s/temp_files/%s/%s_associated_genes_no_NA_no_overlap.bed\'"%(TargetOrtho_path,jobID,species,TargetOrtho_path,jobID,species)
    
    c2n="grep \'\tNA\' \'%s/temp_files/%s/%s_associated_genes.bed--no-overlaps\' >\'%s/temp_files/%s/%s_associated_genes_NA_all_no_overlap.bed\'"%(TargetOrtho_path,jobID,species,TargetOrtho_path,jobID,species)
    c2nb="grep -v \'\tNA\tNA\tNA\tNA\' \'%s/temp_files/%s/%s_associated_genes_NA_all_no_overlap.bed\' > \'%s/temp_files/%s/%s_associated_genes_NA_no_overlap.bed\'"  %(TargetOrtho_path,jobID,species,TargetOrtho_path,jobID,species)
    
    
    c3="grep -v \'\tNA\' \'%s/temp_files/%s/%s_associated_exons.bed\' >\'%s/temp_files/%s/%s_associated_exons_no_NA.bed\'"%(TargetOrtho_path,jobID,species,TargetOrtho_path,jobID,species)
    
    c4="grep \'\tNA\' \'%s/temp_files/%s/%s_associated_exons.bed\' >\'%s/temp_files/%s/%s_associated_exons_NA_all.bed\'"%(TargetOrtho_path,jobID,species,TargetOrtho_path,jobID,species)
    c4b="grep -v \'\tNA\tNA\tNA\tNA\' \'%s/temp_files/%s/%s_associated_exons_NA_all.bed\' > \'%s/temp_files/%s/%s_associated_exons_NA.bed\'" %(TargetOrtho_path,jobID,species,TargetOrtho_path,jobID,species)
    
    c3n="grep -v \'\tNA\' \'%s/temp_files/%s/%s_associated_exons.bed--no-overlaps\' >\'%s/temp_files/%s/%s_associated_exons_no_NA_no_overlap.bed\'"%(TargetOrtho_path,jobID,species,TargetOrtho_path,jobID,species)
    
    c4n="grep \'\tNA\' \'%s/temp_files/%s/%s_associated_exons.bed--no-overlaps\' >\'%s/temp_files/%s/%s_associated_exons_NA_all_no_overlap.bed\'"%(TargetOrtho_path,jobID,species,TargetOrtho_path,jobID,species)
    c4nb="grep -v \'\tNA\tNA\tNA\tNA\' \'%s/temp_files/%s/%s_associated_exons_NA_all_no_overlap.bed\' > \'%s/temp_files/%s/%s_associated_exons_NA_no_overlap.bed\'" %(TargetOrtho_path,jobID,species,TargetOrtho_path,jobID,species)
    
    system_call(c1)
    system_call(c2)
    system_call(c2b)
    system_call(c1n)
    system_call(c2n)
    system_call(c2nb)
    system_call(c3)
    system_call(c4)
    system_call(c4b)
    system_call(c3n)
    system_call(c4n)
    system_call(c4nb)

def runBEDOPS(species):
    fimo_bed_path="%s/fimo_out/%s_%s.bed"%(output_dir,jobID,species)
    fimo=convertFimo(species,fimo_bed_path)
    
    #associate genes
    gene_bed_path="%s/gene_coords/mart_export_%s_protein_coding_genes%s.bed" %(TargetOrtho_path,species,versionDic[species])
    associated_gene_out_path="%s/temp_files/%s/%s_associated_genes.bed" %(TargetOrtho_path,jobID,species)
  
    overlap="--no-overlaps"
    command="closest-features --dist --delim '\\t' %s  \'%s.sorted\' \'%s\' > \'%s%s\'" %(overlap,fimo_bed_path,gene_bed_path,associated_gene_out_path,overlap)
    print(command)
    system_call(command)
    
    
    overlap=""
    command="closest-features --dist --delim '\\t' %s  \'%s.sorted\' \'%s\' > \'%s%s\'" %(overlap,fimo_bed_path,gene_bed_path,associated_gene_out_path,overlap)
    print(command)
    system_call(command)
    
    
    #associate exons
    exon_bed_path="%s/gene_coords/mart_export_%s_exons%s.bed" %(TargetOrtho_path,species,versionDic[species])
    associated_gene_out_path="%s/temp_files/%s/%s_associated_exons.bed" %(TargetOrtho_path,jobID,species)
  
    overlap="--no-overlaps"
    command="closest-features --dist --delim '\\t' %s  \'%s.sorted\' \'%s\' > \'%s%s\'" %(overlap,fimo_bed_path,exon_bed_path,associated_gene_out_path,overlap)
    print(command)
    system_call(command)
    
    overlap=""
    command="closest-features --dist --delim '\\t' %s  \'%s.sorted\' \'%s\' > \'%s%s\'" %(overlap,fimo_bed_path,exon_bed_path,associated_gene_out_path,overlap)
    print(command)
    system_call(command)
     
    #separate BEDOPS results that have no left or right feature assocation (fewer columns in results so need to be imported with pandas separately)
    removeNA(species)
    
def runBEDOPS_all(speciesList):
    for species in speciesList:
        runBEDOPS(species)
            
def set_up_gene_associations(speciesList):
    #print("set_up_gene_associations")
    #set up associated gene dfs
    s_dfs=[]
    for species in speciesList:
        s_dfs.append(mk_associated_genes_df(species))
    return s_dfs

def compute_intron_distances(df,species):

    df_gene_coords=pd.read_csv("%s/gene_coords/mart_export_%s_protein_coding_genes%s.bed" %(TargetOrtho_path,species,versionDic[species]),header=None,sep='\t')
    df_gene_coords.columns=["gene_dna","gene_start","gene_end","WBgene","gene_name","gene_strand"]
    df_gene_coords=df_gene_coords[["gene_dna","gene_start","gene_end","WBgene","gene_name"]]
    
    
    df_g=df.merge(df_gene_coords,how="outer", on="gene_name",suffixes=('_hit', '_gene'))
    #merge on index after calcuating distance for gene pos strand, and gene neg strand separately.

    pos_strand=df_g.loc[(df_g["gene_strand"]==1) & (df_g["distance"]==0)]
    neg_strand=df_g.loc[(df_g["gene_strand"]==-1) & (df_g["distance"]==0)]
    intergenic=df_g.loc[(df_g["distance"]>0) ^ (df_g["distance"]<0)]

    pos_strand.loc[:,"distance"]=pos_strand['start'].sub(pos_strand['gene_start'], axis=0)
    neg_strand.loc[:,"distance"]=neg_strand[ 'gene_end'].sub(neg_strand['end'], axis=0)
    #append the updated intragenic distance dfs with the intergenic df
    m=pos_strand.append(neg_strand)
    m=m.append(intergenic)
    #print(m.columns.values,'m.columns.values in compute_intron_distances')

    return m

def mk_associated_genes_df(species):
    genes_no_overlap,genes,exons_no_overlap,exons=import_associated_genes(species)
    genes_no_overlapNA,genesNA,exons_no_overlapNA,exonsNA=import_associated_genesNA(species)

    def reformat_df(df):
        dfL=df[[0,1,2,3,10,11,12]]
        dfR=df[[0,1,2,3,17,18,19]]
        dfL.columns=[["dna","start","end","hit_id","gene_name","gene_strand","distance"]]
        dfR.columns=[["dna","start","end","hit_id","gene_name","gene_strand","distance"]]
        
        df=dfL.append(dfR)

        #set distance dataframe and hit_id to integer
        df.loc[:,"start"]=[str(i) for i in list(df["start"].values)]
        df.loc[:,"end"]=[str(i) for i in list(df["end"].values)]
        df.loc[:,"coord"]=df["dna"].map(str)+ ":" + df["start"] + "-" + df["end"]
        
        #convert start and end coords back to integers.
        df.loc[:,"start"]=[int(i) for i in list(df["start"].values)]
        df.loc[:,"end"]=[int(i) for i in list(df["end"].values)]
        #set distance dataframe and hit_id to integer (take absolute value of distance after upstream and downstream regions are annotated)
        df.loc[:,"distance"]=[int(i) for i in list(df["distance"].values)]
        df.loc[:,"hit_id"]=[int(i) for i in list(df["hit_id"].values)]
        return df
    
    def reformat_dfNA(df):
        #print("reformat_dfNA")
        dfL=df.loc[~df[6].isnull()][[0,1,2,3,10,11,12]]
        dfR=df.loc[~df[14].isnull()][[0,1,2,3,12,13,14]]
        dfL.columns=[["dna","start","end","hit_id","gene_name","gene_strand","distance"]]
        dfR.columns=[["dna","start","end","hit_id","gene_name","gene_strand","distance"]]        
        df=dfL.append(dfR)
        
        #set distance dataframe and hit_id to integer
        df.loc[:,"start"]=[str(i) for i in list(df["start"].values)]
        df.loc[:,"end"]=[str(i) for i in list(df["end"].values)]
        df.loc[:,"coord"]=df["dna"].map(str)+ ":" + df["start"] + "-" + df["end"]
        
        #convert start and end coords back to integers.
        df.loc[:,"start"]=[int(i) for i in list(df["start"].values)]
        df.loc[:,"end"]=[int(i) for i in list(df["end"].values)]
        #set distance dataframe and hit_id to integer (take absolute value of distance after upstream and downstream regions are annotated)
        df.loc[:,"distance"]=[int(i) for i in list(df["distance"].values)]
        df.loc[:,"hit_id"]=[int(i) for i in list(df["hit_id"].values)]
        return df
    
    genes_no_overlap=reformat_df(genes_no_overlap)
    genes=reformat_df(genes)
    exons_no_overlap=reformat_df(exons_no_overlap)
    exons=reformat_df(exons)
    
    genes_no_overlapNA=reformat_dfNA(genes_no_overlapNA)
    genesNA=reformat_dfNA(genesNA)
    exons_no_overlapNA=reformat_dfNA(exons_no_overlapNA)
    exonsNA=reformat_dfNA(exonsNA)
    
    genes_no_overlap=genes_no_overlap.append(genes_no_overlapNA)
    genes=genes.append(genesNA)
    exons_no_overlap=exons_no_overlap.append(exons_no_overlapNA)
    exons=exons.append(exonsNA)
    
    #intragenic_hits from gene df
    intragenic=genes.loc[genes["distance"]==0]
    intragenic_hits_list=set(list(intragenic["hit_id"].values))
       
    #true intergenic only hits list ( excludes hits that are in intron or exon of adjacent gene -intragenic_hits)
    intergenic=genes_no_overlap.loc[~genes_no_overlap["hit_id"].isin(intragenic_hits_list)]

    #upstream and downstream hits
    #dist < 0 and strand =-1, hit is upstream
    #if dist >0, and strand=1, hit is upstream
    upstream1=intergenic.loc[(intergenic["distance"]<0) & (intergenic["gene_strand"]==-1)]
    upstream2=intergenic.loc[(intergenic["distance"]>0) & (intergenic["gene_strand"]==1)]
    upstream=upstream1.append(upstream2)

    #get rid of max_search distance results
    #first make sure all distances are positive
    upstream.loc[:,"distance"]=[abs(int(i)) for i in list(upstream["distance"].values)]
    if max_search_distance!=None:
        upstream=upstream.loc[upstream["distance"]<int(max_search_distance)]
    upstream.loc[:,"region"]="upstream"

    #dist <0 and strand =1, hit is downstream
    #if dist >0, and strand=-1, hit is downstream
    downstream1=intergenic.loc[(intergenic["distance"]<0) & (intergenic["gene_strand"]==1)]
    downstream2=intergenic.loc[(intergenic["distance"]>0) & (intergenic["gene_strand"]==-1)]
    downstream=downstream1.append(downstream2)
    downstream.loc[:,"distance"]=[abs(int(i)) for i in list(downstream["distance"].values)]
    if max_search_distance!=None:
        downstream=downstream.loc[downstream["distance"]<int(max_search_distance)]
    downstream.loc[:,"region"]="downstream"
    #determine which intragenic hits are exonic or intronic (get exon file with and without overlap
    #get list of hit ids that have dist of 0. these are exonic hits. all other intragenic hit ids are intronic

    #exonic hits
    exons_only=exons.loc[exons["distance"]==0]
    exon_hits_list=list(exons_only["hit_id"].values)
    exons_only.loc[:,"region"]="exon"

    #intronic hits (not in exon hits list but intragenic)
    introns=intragenic.loc[~intragenic["hit_id"].isin(exon_hits_list)]
    introns=introns.loc[introns["distance"]==0]
    #introns=intragenic hits (dist=0 in genes df) that are not exonic hits (dist=0 in exons df).
    #get list of intragenic hit ids that are not in exonic hit id list.
    
    introns.loc[:,"region"]="intron"
    #merge all dataframes into one associated_genes dataframe
    dfs = [upstream, downstream, introns, exons_only]
    associated_genes = pd.concat(dfs)
    
    #assigne intron distances
    associated_genes=compute_intron_distances(associated_genes,species)
    
    #set all distance to absolute value
    associated_genes.loc[:,"distance"]=[abs(i) for i in list(associated_genes["distance"].values)]
    #restrict results by max search distance  --need to compute average distances after introns distances are present.
    if max_search_distance!=None:
        associated_genes=associated_genes.loc[associated_genes["distance"]<int(max_search_distance)]
    
    return associated_genes

def is_non_zero_file(fpath):  
    return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 0 else False

def import_associated_genes(species):
    
    f="%s/temp_files/%s/%s_associated_genes_no_NA_no_overlap.bed" %(TargetOrtho_path,jobID,species)
    if is_non_zero_file(f):genes_no_overlap=pd.read_csv(f,sep='\t',header=None)
    else: genes_no_overlap=pd.DataFrame(index=None, columns=[list(range(20))])
        
    f="%s/temp_files/%s/%s_associated_genes_no_NA.bed"%(TargetOrtho_path,jobID,species)
    if is_non_zero_file(f):genes=pd.read_csv(f,sep='\t',header=None)
    else:genes=pd.DataFrame(index=None, columns=[list(range(20))])
        
    f="%s/temp_files/%s/%s_associated_exons_no_NA.bed"%(TargetOrtho_path,jobID,species)
    if is_non_zero_file(f):exons=pd.read_csv(f,sep='\t',header=None)
    else:exons=pd.DataFrame(index=None, columns=[list(range(20))])
    
    f="%s/temp_files/%s/%s_associated_exons_no_NA_no_overlap.bed"%(TargetOrtho_path,jobID,species)
    if is_non_zero_file(f):exons_no_overlap=pd.read_csv(f,sep='\t',header=None)
    else: exons_no_overlap=pd.DataFrame(index=None, columns=[list(range(20))])
    
    
    
    return genes_no_overlap,genes,exons_no_overlap,exons


def import_associated_genesNA(species):
    #print("import_associated_genesNA")

    f="%s/temp_files/%s/%s_associated_genes_NA_no_overlap.bed" %(TargetOrtho_path,jobID,species)
    if is_non_zero_file(f):genes_no_overlapNA=pd.read_csv(f,sep='\t',header=None)
    else: genes_no_overlapNA=pd.DataFrame(index=None, columns=[list(range(20))])
    
    f="%s/temp_files/%s/%s_associated_genes_NA.bed"%(TargetOrtho_path,jobID,species)
    if is_non_zero_file(f):genesNA=pd.read_csv(f,sep='\t',header=None)
    else:genesNA=pd.DataFrame(index=None, columns=[list(range(20))])
    
    f="%s/temp_files/%s/%s_associated_exons_NA.bed"%(TargetOrtho_path,jobID,species)
    if is_non_zero_file(f):exonsNA=pd.read_csv(f,sep='\t',header=None)
    else: exonsNA=pd.DataFrame(index=None, columns=[list(range(20))])
    
    f="%s/temp_files/%s/%s_associated_exons_NA_no_overlap.bed"%(TargetOrtho_path,jobID,species)
    if is_non_zero_file(f):exons_no_overlapNA=pd.read_csv(f,sep='\t',header=None)
    else: exons_no_overlapNA=pd.DataFrame(index=None, columns=[list(range(20))])
   
    return genes_no_overlapNA,genesNA,exons_no_overlapNA,exonsNA


def set_up_merged_gene_fimo_data(speciesList,s_dfs,ortho_dfs):
    fimo_dfs=[]
    s_m_dfs=[]
    
    for species in speciesList:
        fimo=pd.read_csv("%s/fimo_out/%s_%s/fimo.txt" %(output_dir,jobID,species),sep='\t')
        fimo_dfs.append(fimo)

    #merge reference species data
    s1=s_dfs[0]
    s1_m=s1.merge(fimo_dfs[0],how="inner",right_index=True, left_on="hit_id")
    s1_m["ref_gene_name"]=list(s1_m["gene_name"].values)
    s1_m=s1_m[["ref_gene_name","hit_id","region","distance","coord","strand","score","p-value","q-value","matched_sequence"]]
    s_m_dfs.append(s1_m)
    
    #merge using inner join so that only only line per hit_id is present ()
    how_string="inner"
    for i in range(len(s_dfs)-1):
        s_m=mergeOrtho_fimo(s_dfs[i+1],ortho_dfs[i],fimo_dfs[i+1],how_string)
        s_m_dfs.append(s_m)

    for i in range(len(s_m_dfs)):
        print("%s: %s motif matches" %(speciesList[i],len(s_m_dfs[i])))
    return s_m_dfs,fimo_dfs


def mergeOrtho_fimo(df,ortho,fimo,how_string):

    m=df.merge(ortho, how=how_string,on="gene_name")
    m=m.merge(fimo,how=how_string,right_index=True, left_on="hit_id")

    m=m[["hit_id","gene_name","ref_gene_name","ht","region","distance","coord","strand","score","p-value","q-value","matched_sequence"]]
    m.columns=["hit_id","gene_name","ref_gene_name","ht","region","distance","coord","strand","score","p-value","q-value","matched_sequence"]
    return m   

def computeStats(s_m_dfs):
    #print("computeStats")
    s1_m=s_m_dfs[0]
    region="upstream"
    ref_genes_up=list(set(list(set(list(s1_m.loc[s1_m["region"]==region]["ref_gene_name"].values)))))
    l_upstream=float(len(ref_genes_up))
    print('total %s genes with %s motif matches: %s' %(speciesList[0],region,int(l_upstream)))
    
    region="intron"
    ref_genes_intron=list(set(list(set(list(s1_m.loc[s1_m["region"]==region]["ref_gene_name"].values)))))
    l_intron=float(len(ref_genes_intron))
    print('total %s genes with %s motif matches: %s' %(speciesList[0],region,int(l_intron)))
    
    l=l_upstream+l_intron
    
    def computeRegion(region,count,ref_genes,l):
        #make and empty df and fill it line by line
        results_DF=pd.DataFrame()
        for g in ref_genes[:]:
            count+=1
            if count%100==0:
                print("%.2f percent of genes processed" %(count/l*100))
            dfs=[]
            for i in range(len(s_m_dfs)):
                d=getGeneInfo(s_m_dfs[i],g,region)
                dfs.append(d)

            #conservation
            binary=[1 if len(d)>0 else 0 for d in dfs]
            #get species list from binary
            cons_species=[]
            for i in range(len(binary)):
                if binary[i]>0:
                    cons_species.append(speciesList[i])
            conservation=sum((binary))

            #build list of species stats 
            dfs_non_empty=[]
            for df in dfs:
                if len(df)>0:
                    dfs_non_empty.append(getSpeciesStats(df))
                    #get hit info

            g_stats_matrix=np.matrix(dfs_non_empty)

            avg_species_stats=np.sum(g_stats_matrix, axis=0)/len(dfs_non_empty)

            line_list=[g,region,conservation]


            for i in range(5):
                line_list.append(avg_species_stats.item(0,i))
            for sp in range(0,conservation):
                line_list.append(cons_species[sp])
                for i in range(5):
                    line_list.append(g_stats_matrix.item(sp,i))
            #fill in emtpy columns to keep dataframe same size
            for i in range(conservation,len(speciesList)):
                for i in range(6):
                    line_list.append("NA")
            line_dic={}
            for i in range(len(line_list)):
                line_dic[i]=line_list[i]
            
            results_DF=results_DF.append(line_dic, ignore_index = True)
        return results_DF,count
    count=1
    results_DF_upstream,count=computeRegion("upstream",count,ref_genes_up,l)
    results_DF_intron,count=computeRegion("intron",count,ref_genes_intron,l)
        
    return results_DF_upstream,results_DF_intron

def getGeneInfo(df,g,region):
    df=df.loc[(df["ref_gene_name"]==g) & (df["region"]==region)]
    return df

def getMaxAvg(df,header):
    #print("getMaxAvg")
    scores=list(df[header].values)
    if len(scores)>0:
        if header=="distance":
            max_score=min(scores)
        else:
            max_score=max(scores)
        avg_score=sum(scores)/float(len(scores))
    else:
        max_score=np.nan
        avg_score=np.nan
    return max_score,avg_score

def getSpeciesStats(d):
    if len(d)>0:
        max_score,avg_score=getMaxAvg(d,"score")
        max_dist,avg_dist=getMaxAvg(d,"distance")
        site_count=len(d)
    else:
        max_score=np.nan
        avg_score=np.nan
        max_dist=np.nan
        avg_dist=np.nan
        site_count=np.nan
    return np.array([max_score,avg_score,max_dist,avg_dist,site_count])


def merge_upstream_intron(speciesList,results_up,results_intron):
    #merge upstream and intron gene info 
    #rename column headers
    region="upstream"
    region="intron"

    h=["ref_gene_name","region","conservation","avg_species_max_PSSM_score","avg_species_avg_PSSM_score",
       "avg_species_min_dist","avg_species_avg_dist","avg_species_site_count"]
    
    count=0
    for species in speciesList[:]:
        count+=1
        h=h+["species%s" %count,"species%s_max_PSSM_score" %count,"species%s_avg_PSSM_score" %count,
             "species%s_min_dist" %count,"species%s_avg_dist" %count,"species%s_site_count" %count]
  
    results_up.columns=h
    results_intron.columns=h
    
    #merge on gene name only, have region1 and region2 available for each, then sum site counts and avg score to get gene info
    all_info=results_up.merge(results_intron, how='outer',on="ref_gene_name",suffixes=('_upstream', '_intron'))
    #write all_info file to table for testing best algorithm to use to generate cumulative score based on training data
    #all_info.to_csv("%s/%s_all_info.txt" %(output_dir,jobID),sep='\t',index=False)
    return all_info

def writeSpecies_hit_files(s_dfs,fimo_dfs,ortho_dfs):
    
    #Make all hits file for each species 
    def sort_by_ref_gene(df):
        try:df=df.sort_values(["ref_gene_name","region"],ascending=[True,False])
        except:df=df.sort_values(["gene_name","region"],ascending=[True,False])
        return df
    #merge on fimo and ortho u
    s_m_all_hits_dfs=[]
    s1=s_dfs[0]
    s1_m_all_hits=s1.merge(fimo_dfs[0],how="outer",right_index=True, left_on="hit_id")
    #print(s1_m_all_hits.columns.values,'s1_m_all_hits')
    s1_m_all_hits=sort_by_ref_gene(s1_m_all_hits)
    
    s_m_all_hits_dfs.append(s1_m_all_hits)
    
    for i in range(1,len(speciesList)):
        s_m_all_hits=mergeOrtho_fimo(s_dfs[i],ortho_dfs[i-1],fimo_dfs[i],"inner")
        s_m_all_hits=sort_by_ref_gene(s_m_all_hits)
        s_m_all_hits.columns=["hit id","gene name","reference genome ortholog","homology type","region","distance from gene", "coordinate", "strand", "score", "fimo p-value","fimo q-value", "matched_sequence"]

        s_m_all_hits.to_csv("%s/motif_match_data_per_species/%s_%s_genome_motif_match_results.csv" %(output_dir,jobID,speciesList[i]),sep=',',index=False)
    
        s_m_all_hits_dfs.append(s_m_all_hits)

    #write each species table to csv
    s1_m_all_hits.columns=['dna', 'start_x', 'end' ,'hit_id' ,'ref_gene_name', 'gene_strand' ,'distance from gene',
     'coordinate', 'region', 'gene_dna' ,'gene_start' ,'gene_end', 'WBgene', '# motif_id',
     'motif_alt_id' ,'sequence_name' ,'start_y', 'stop', 'strand', 'score' ,'fimo p-value',
     'fimo q-value' ,'matched_sequence']

    s1_m_all_hits=s1_m_all_hits[["hit_id","ref_gene_name","region","distance from gene", "coordinate", "strand", "score", "fimo p-value","fimo q-value", "matched_sequence"]]

    s1_m_all_hits.to_csv("%s/motif_match_data_per_species/%s_%s_genome_motif_match_results.csv" %(output_dir,jobID,speciesList[0]),sep=',',index=False)    
    return s_m_all_hits_dfs

def fill_NA(l,x):
    n=[]
    for i in l:
        if i=='NA':
            n.append(x)
        elif i=='nan':
            n.append(x)
        else:
            n.append(i)
    return n

def min_max_species_data(species_dfs,df2):
    all_scores=[]
    all_dist=[]
    all_counts=[]
    for df in species_dfs:
        all_scores=all_scores+list(df["score"].values)
        all_dist=all_dist+list(df["distance from gene"].values)
    
    for h in list(df2.columns.values):
        if 'count' in h:
            all_counts=all_counts+list(df2[h].values)
            
    all_counts=fill_NA(all_counts,0)
    
    m1=min(all_scores)
    mx1=max(all_scores)
    m2=min(all_dist)
    mx2=max(all_dist)
    max_count=max(all_counts)
    return m1,mx1,m2,mx2,max_count

def norm_data(df,min_scores,max_scores,min_dist,max_dist,max_count):
    #norm each column of data using min and max
    headers=getHeaders(df)
    #fill NA values with with min score, or max dist, or min conservation...
    
    def norm_list(minl,maxl,l):
        l_normed=[(i-minl)/(maxl-minl) for i in l]
        return l_normed
    
    for h in headers:
        l=list(df[h].values)
        if 'score' in h:
            l=fill_NA(l,min_scores)
            l_normed=norm_list(min_scores,max_scores,l)
            df[h]=l_normed
        if 'count' in h:
            l=fill_NA(l,0)
            l_normed=norm_list(0,max_count,l)
            df[h]=l_normed
        if 'dist' in h:
            l=fill_NA(l,max_dist)
            l_normed=norm_list(min_dist,max_dist,l)
            l_normed=[(1-i) for i in l_normed]
            df[h]=l_normed
        if 'conservation' in h:
            l=fill_NA(l,0)
            l_normed=norm_list(0,5,l)
            df[h]=l_normed
    df=df.fillna(value=0)
    #df.to_csv("%s/%s_all_info_normed.txt" %(output_dir,jobID),sep='\t',index=False)
    
    return df

def getHeaders(df):
    headers=[]
    for i in list(df.columns.values):
        if 'avg_species_' not in i:
            if 'dist' not in i:
                if 'score' in i:
                    headers.append(i)
                if 'count' in i:
                    headers.append(i)
                if 'conservation' in i:
                    headers.append(i)
    return headers


def import_orthologs(speciesList):
    print("import_orthologs")
    ortho_dfs=[]
    for i in range(1,len(speciesList)):
        o=pd.read_csv("%s/ortholog_files/%s_%s_ortho.txt" %(TargetOrtho_path,speciesList[0],speciesList[i]),sep="\t",header=0)
        o.columns=["ref_gene_name","gene_name","ht"]
        ortho_dfs.append(o)
    return ortho_dfs

def getDF_by_region(df,region):
    print("getDF_by_region")
    df_reg=df.loc[df["region"]==region]
    return df_reg

def getMatches(g,o,df,result_column):
    print("getDF_by_region")
    orthos=list(set(list(o.loc[o["g1"]==g]["g2"].values)))
    matches=list(set(list(df.loc[df["gene_name"].isin(orthos)][result_column].values)))
    return matches

def mergeSpecies(m1,m2):
    print("mergeSpecies")
    m=m1.merge(m2,how="inner",on=["ref_gene_name","region"])
    return m

def compute_gene_info(region,sm_dfs):
    print("compute_gene_info")
    s1_m=sm_dfs[0]
    ref_genes=list(set(list(set(list(s1_m.loc[s1_m["region"]==region]["ref_gene_name"].values)))))
    df=writeResults(ref_genes,region,s_m_dfs)
    return df



def normList(df,header): 
    print("normList")
    l=list(df[header].values)
    maxl=float(max(l))
    minl=float(min(l))
    l_normed=[(i-minl)/(maxl-minl) for i in l]
    return l_normed


def adjust_probs(probs):
    print("adjust_probs")
    """takes a list of classifier probabilities of form (class 0 probability, class 1 probability) 
    and transforms class 0 probabilites by multiplying by -1. This results is high probabilty of class 1
    having the largest values and high probability of class 0 having the lowerst values"""
    plist=[]
    for i in range(len(probs)):
        ps= str(probs[i]).split()
        B=float(ps[-1][:-1])
        A=float(ps[-2])
        if A>B:plist.append(float(A*-1))
        else:plist.append(float(B))
    return plist

def rankData(data):
    print("rankData")
    """takes a list of numbers of returns the list in place with elements replaced by rank order"""
    uni=list(set(data))
    uni.sort()
    uni.reverse()
    uniDic={}
    for i in range(len(uni)):
        uniDic[uni[i]]=i+1
    rankedList=[uniDic[n] for n in data]
    return rankedList


def predict_targets(df_normed):
    #train model on all unc-3 data, classify whole genome results, assess GO terms from ranked results 

    #import training data
    if len(speciesList)==5:
        train=pd.read_csv("%s/data/unc3_ASE_AIY_training_data.txt" %(TargetOrtho_path), sep='\t')
    #if all 8 species included use this training set
    if len(speciesList)==8:
        train=pd.read_csv("%s/data/unc3_ASE_AIY_training_dataj_unc3_p0007478j_ASE_p000424403j_AIY_p000257719_8species.txt" %(TargetOrtho_path),sep='\t')
    headers=getHeaders(train)
    SVM=svm.SVC(kernel='rbf',probability=True,class_weight='balanced')
    model=SVM
    feature_set=headers
    X_train=train[feature_set]
    y_train=train["class_label"]

    #predict dataset class probability based on training data (unc-3, AIY, ASE combined data)
    X_test=df_normed[feature_set]
    probs =model.fit(X_train, y_train).predict_proba(X_test)
    y_preds=model.fit(X_train, y_train).predict(X_test)
    df_normed["class_pred"]=y_preds

    #adjust class probabilities for ranking data 
    #(multiply class 0 probs by -1 so that highest values are most probable class 1, and lowest values are more probable class 0)
    adj_probs=adjust_probs(probs)
    df_normed.loc[:,"class_prob"]=adj_probs
    #rank data by class predictions (TF target gene or unknown )
    ranks=rankData(adj_probs)
    return adj_probs,ranks

def updateDF(df,adj_probs,ranks):
    #add classifier prediction probability based ranks to data and write to file.
    df.loc[:,"Rank"]=ranks
    df.loc[:,"Rank_Pct"]=df.Rank.rank(pct=True)*100
    df.loc[:,"class_prob"]=adj_probs
    df=df.sort_values(by=["Rank"],ascending=[1])

    headers=list(df.columns.values)
    new_headers=["ref_gene_name","Rank","Rank_Pct","class_prob","conservation_upstream","conservation_intron"]
    
    for c in range(1,len(speciesList)+1):
        for h in headers:
            if "species%s" %c in h:
                if "upstream" in h:
                    new_headers.append(h)
    for c in range(1,len(speciesList)+1):
        for h in headers:
            if "species%s" %c in h:
                if "intron" in h:
                    new_headers.append(h)
    for h in headers:
        if 'avg_species' in h:
            if 'upstream' in h:
                new_headers.append(h)
    for h in headers:
        if 'avg_species' in h:
            if 'intron' in h:
                new_headers.append(h)
                
    df=df[new_headers]

    #write the results to a file
    df.to_csv("%s/%s_TargetOrtho2_ranked_genes_summary.csv" %(output_dir,jobID),sep=',',index=False)


def main():
    try:
        print(jobID)
        mkDirs(jobID,output_dir)
        writeParams()
        runFIMO_threaded(speciesList,jobID)
        runBEDOPS_all(speciesList)
        s_dfs=set_up_gene_associations(speciesList)
        clear_temp()
        ortho_dfs=import_orthologs(speciesList)
        s_m_dfs,fimo_dfs=set_up_merged_gene_fimo_data(speciesList,s_dfs,ortho_dfs)

        upstream_df,intron_df=computeStats(s_m_dfs)

        all_info=merge_upstream_intron(speciesList,upstream_df,intron_df)
        ##all_info=compute_combo_features(all_info,speciesList)
        species_dfs=writeSpecies_hit_files(s_dfs,fimo_dfs,ortho_dfs)

        #norm data for scoring
        min_scores,max_scores,min_dist,max_dist,max_count=min_max_species_data(species_dfs,all_info)
        all_info_normed=norm_data(all_info,min_scores,max_scores,min_dist,max_dist,max_count)
        adj_probs,ranks=predict_targets(all_info_normed)
        all_info=merge_upstream_intron(speciesList,upstream_df,intron_df)
        updateDF(all_info,adj_probs,ranks)
        print("output directory: \'%s\'" %(output_dir))
    except Exception as e: 
        print(e)
        #clear entire output directory and temp file.
        #clear_error()
main()





