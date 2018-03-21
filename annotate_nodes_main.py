#!/usr/bin/python

#abigailc@Leviathan-Whale February 14 2018 - began
#abigailc@Leviathan-Shark March 19th 2018 - updated

#overall idea

#I have a large number of phylogenetic trees. 
#for each subclade (as defined by a taxonomic search algorithm) of a massive tree of SOD I have
#   1 species tree (ribosomal concat) - bestTree made in RAxML, rooted at a specific outgroup sequence
#   1 gene tree (SOD) - with 100 bootstraps made in RAxML, unrooted or midpoint rooted
#   1 (theoretical) reconciliation file produced by RANGER DTL (these need to be run, it will take forever manually)

#ideal output: a species tree with annotated nodes. a gene tree with annotated nodes. a summary file.
#annotations = #bs support, # times IDed as transfer, #times IDed as Duplication. etc

#what about where topology doesn't match? -> look at the BS file itself.
#                                         -> i am defining clades by all of their downstream tips, so exact topology is irrelevent



###################
#testing log
"""
get sample of my data




"""
############################

#imports real packages
import glob
import re
import os
#imports my modules (?)
import annotate_classes
import annotate_ranger
import annotate_parseranger

###############
#functions

#main
#march2018
def overall_annotation_function(path_to_species_trees, path_to_gene_trees, project_name = "myproject", path_to_ranger_outputs = "", ):
    """
    this function should take in folders (full of gene and species trees built from several datasets)
    it should then correlate them to each other by name of the dataset
    compare the data using RangerDTL (running all 100 bootstraps of the gene tree individually against the species tree consensus for each dataset)
    parse the output of RangerDTL
    map the output onto each clade of the species tree
    map the output onto each clade of each bootstrap tree
    write each to file so we can visually see what's going on
    write the data from each to tab-seperated-file so we can do stats on what's going on
    --> visualization of HGT events + congruence of gene and species tree
    
    args:
        path_to_species_trees : path to the folder you have your species trees saved in
        path_to_gene_trees : path to the folder you have your gene trees saved in
    optional_args:
        project_name : the name of this project (was used to generate the files i currently have)
    
    output:
        one output file per dataset identified -> should be a mapping of reticulations and bootstraps onto the best species tree.
    """

    #initially gather the names of the datasets from the species_trees folder
    dataset_names = gather_dataset_names(path_to_species_trees, project_name, "_CC")
    #create an object of class PhyloData for each unique dataset
    phylodata_objects = []
    for name in dataset_names:
        phylodata_objects.append(annotate_classes.PhyloData(name))
    #for each object, have it try and assign itself the correct files
    print("populating phylodata objects")
    populate_objects(phylodata_objects, project_name, path_to_species_trees, path_to_gene_trees, path_to_ranger_outputs)
    #run the visualizer for each object
    parse_and_visualize(phylodata_objects, project_name)

##############3
#subfunctions

#march2018
def parse_and_visualize(phylodata_objects, projectname):
    """
    this function runs the parse_ranger_clade function (as seen in the module annotate_parseranger)
    once for each dataset/phylodata object

    args: 
        phylodata_objects : list of phylodata objects, each must be populated with rangerDTL output files (100 bs) and the species tree
        projectname : the name of this project, for naming conventions
    returns:
        nothing
    outputs:
        four files per dataset:
        -annotated specties tree
        -annotated bootstraps file
        -tab-sep species retics
        -tab-sep bootstrap retics
    """
    print("beginning to parse ranger outputs and annotate trees")
    os.system("mkdir results")
    for dataset in phylodata_objects:
        datasetname = dataset.name
        with open (dataset.species_tree) as sp:
            str_species_tree = sp.read().strip()
        success = annotate_parseranger.parse_ranger_clade(dataset.ranger_dtl_outputs, "results/"+projectname+datasetname, "results/"+projectname+datasetname, str_species_tree)
        print("finished parsing and annotating the dataset "+datasetname)
    print("finished parsing and annotating all data. look in ./results to see results.")
    
#march2018
def populate_objects(phylodata_objects, project_name, path_to_species_trees, path_to_gene_trees, path_to_ranger_outputs):
    """
    this function will try and associate each phylodata object with the correct
    species_besttree
    gene_bootstrap_trees
    and rangerDTL output files (if they exist)
    args:
        list of phylodata objects
        name of project
        paths (to species trees, to bootstrap gene trees, to rangerDTL
    returns
        True if everything was associated
        False if something has gone horribly awry

    """


    #try and populate the species and gene files. should work.
    for obj in phylodata_objects:
        #print("Populating species trees")
        obj.populate_species_tree(path_to_species_trees)
        #print("Populating gene trees")
        obj.populate_gene_boots(path_to_gene_trees)


    #now try and populate ranger output, if not make directory and run run_rangerDTL
    for obj in phylodata_objects:
        #print("Checking for rangerDTL outputs")
        exists = obj.populate_ranger_dtl_outputs(path_to_ranger_outputs)
        if exists is False:
            #run the program.
            print("Running RangerDTL")
            path_to_ranger_outputs, list_of_ranger_outputs = annotate_ranger.run_rangerDTL(obj, project_name)
            #print("Checking for new rangerDTL outputs")
            exists = obj.populate_ranger_dtl_outputs(path_to_ranger_outputs)
            if exists is False:
                print ("error in rangerdtl_output assignation")
                raise SystemExit
    return True

#helper functions
#march2018
def gather_dataset_names(path_to_files, project_name, suffix):
    """this function will populate a list of the files in /species_trees/
    and then it will use REGEX to extract the name of the dataset associated with each
    for use in associating all of the files together. (#TODO maybe i can use sqlite for this?)
    
    please ensure that all files in the given folder are species_tree bestTrees from RAxML
    for this specific project (this is default output of DTL_PROJECT.py)
    """
    species_trees_list = list_files_in_given_folder(path_to_files)
    dataset_names = extract_dataset_names(species_trees_list, project_name, suffix)
    return dataset_names

#march2018
def extract_dataset_names(list_of_files, prefix = "", suffix = ""):
    """ this function takes in a list of files, and uses regex to extract the name of each
    dataset (this was not saved for some reason in initial processing)
    
    args:
        list_of_files : a list of filenames to extract dataset_names from
    optional inputs include:
        prefix : (which was used in the generation of these files, eg SOD1 or squalinesynth)
        suffix : (was occasionally used in generation of these files, eg _gene or _species)
    outputs:
        dataset_names: a list of the extracted dataset names
    #this is a docstring test (see if it works)
    >>> extract_dataset_names(["RAxML_bestTree.SOD_test_species"], "SOD_", "_species")
    ["test"]

    """
    dataset_names = []
    for filename in list_of_files:
        dataname = re.sub("(.*?)(\.)("+prefix+")(.*)("+suffix+")", "\\4", filename)
        dataset_names.append(dataname)
    assert dataset_names != [], "dataset_names did not populate"
    return dataset_names

#march2018
def list_files_in_given_folder(path_to_folder):
    """ this takes in a path, and should return a list of all of the files
    found in that folder, where each file name is a string.

    args
        path_to_folder : a path from the location you are running from to the folder you care about.
    returns
        a list of strings. strings are the files within that folder.
    """
    file_names_list = []
    for file_name in glob.glob(path_to_folder+"/*"):
        file_names_list.append(file_name)
    assert file_names_list != [], "failed to populate folder"+path_to_folder
    return file_names_list




#argument parser so i can run from command line?


if __name__ == "__main__":

    print("Running in terminal")
    #imports real packages
    import glob
    import re
    import os
    import sys
    import argparse
    #imports my modules (?) #they have to be in the same folder right now.
    import annotate_classes
    import annotate_ranger
    import annotate_parseranger

    parser = argparse.ArgumentParser(description="All")

    #all arguments are 
    parser.add_argument("directory", nargs='?', default=os.getcwd(), type=str, help="type name of directory to run in eg SOD_project")
    parser.add_argument("-p", "--projectname", action = "store", default = "myproject", help="type projectname eg SOD_ver2")
    parser.add_argument("-b", "--bootstrap_path", action = "store", default= os.getcwd(), help="type path to your bootstrap files")
    parser.add_argument("-s", "--species_path", action = "store", default= os.getcwd(), help="type path to your speciestree bestTree files")
    parser.add_argument("-r", "--rangerout_path", action = "store", default="", help="type path to your rangerDTL files if you have them")

    args = parser.parse_args()

    os.chdir(args.directory)
    print("moved to "+args.directory)

    overall_annotation_function(path_to_species_trees=args.species_path, path_to_gene_trees=args.bootstrap_path, project_name=args.projectname, path_to_ranger_outputs=args.rangerout_path)
   