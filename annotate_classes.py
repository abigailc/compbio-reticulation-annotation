#module for annotate_nodes_main.py
#update abigailc@leviathan april 2 2018

import glob

class PhyloData(object):
    """this class should keep track of the associated files for a
    given instance

    optional settings
        name_of_dataset : should be alphanumeric identifier of this specific dataset
        species_tree : should be the raxml species tree output (best consensus tree)    
        gene_tree : should be the raxml 100-bootstraps output for the gene tree associated with this dataset
        bipartitions_table : should be the output file from my shitty code - redo this
        ranger_dtl_output : should be the output file for this dataset from RangerDTL
        output_tree: where the final output will be written to
    """
    def __init__(self, name_of_dataset):
        """ this just initialized the parameters I expect to need
            and gives the instance a name """
        self.name = name_of_dataset
        self.species_tree = None
        self.gene_boots = None
        self.bipartitions_table = None
        self.ranger_dtl_outputs = None
        self.who_are_you()
    #QUESTION #TODO: is it possible to have written all populate_X functions as one single
    #function with a line like self.VARIABLE = blah, where variable is an arg to this method like "species_tree"
    def populate_species_tree(self, path="", require = "bestTree",):
        """
        given a path, look in that folder for a file that contains your name
        verify that it is the type of file you want using require (require = "bestTree", "bipartitions", etc)
        args: 
            path : where to look
        optionalargs:
            require : added verification step that this string is in the found file name
            path : path to go to when looking for files.
        return:
            the name of the file you found (which is also saved as an instance variable)
        """
        file_list = []
        file_list = list_files_in_given_folder(path)
        for file in file_list:
            if self.name in file:
                if require in file:
                    #self.variable
                    self.species_tree = file
                    return file
        print ("error, could not assign a species tree to dataset: "+self.name)
        raise SystemExit
    
    def populate_gene_boots(self, path="", require = "bootstrap"):
        """
        given a path, look in that folder for a file that contains your name
        verify that it is the type of file you want using require (require = "bestTree", "bipartitions", etc)
        args: 
            path : where to look
        optionalargs:
            require : added verification step that this string is in the found file name
            path : path to go to when looking for files.
        return:
            the name of the file you found (which is also saved as an instance variable)
        """
        file_list = []
        file_list = list_files_in_given_folder(path)
        for file in file_list:
            if self.name in file:
                if require in file:
                    self.gene_boots = file
                    return file_list
        print ("error, could not assign a gene_bootstrap file to dataset: "+self.name)
        raise SystemExit

    def populate_ranger_dtl_outputs(self, path = "", require = "RangerOut"):
        """
        given a path, look in that folder for a file that contains your name
        verify that it is the type of file you want using require (require = "bestTree", "bipartitions", etc)
        args: 
            path : where to look
        optionalargs:
            require : added verification step that this string is in the found file name
            path : path to go to when looking for files.
        return:
            the list of the files you found (which is also saved as an instance variable)
        """
        file_list = []
        #this should contain 100 output files
        output_list = []
        file_list = list_files_in_given_folder(path)
        for file in file_list:
            if self.name in file:
                if require in file:
                    output_list.append(file)
        if len(output_list) != 100:
            print(str(len(output_list))+" ranger_dtl_outputs were found for "+self.name)
            return False
        self.ranger_dtl_outputs = output_list
        return output_list
        print ("could not assign a rangerDTL output file to: "+self.name)
        return False
        #TODO: implement running the dtl on anything this errors

    def who_are_you(self):
        #print your name :)
        print("I am "+self.name)

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