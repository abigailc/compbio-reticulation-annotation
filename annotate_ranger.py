
#update abigailc@leviathan march 20

import os
import re

def run_rangerDTL(obj, project_name):
	"""
	this is the main function called in this module
	it should take a phylodata object and produce rangerDTL outputs
	by seperating the bootstraps, writing them each to a ranger input file,
	running all the optroots, running all the rangers, saving and returning the outputs.

	args:
		obj : a phylodata object with attributes name, species_tree, and gene_boots assigned
		project_name : whatever you want to call this project. 
	"""
	list_of_inputs = []
	bootfile = obj.gene_boots
	speciesfile = obj.species_tree

	current_obj_name = obj.name
	#this will be a list of strings, where each string is one bootstrap tree in .newick form
	sep_boots_list = separate_bootstrap_trees(bootfile)
	os.system("mkdir ranger_inputs")
	list_of_inputs, species_string = write_inputs(sep_boots_list, speciesfile, project_name+"_"+current_obj_name)
	ranger_outputs = run_programs(list_of_inputs, species_string, project_name+"_"+current_obj_name)
	undo_replace_underscores(ranger_outputs)
	print("generated rangerDTL output files at ranger_outputs/")
	#note that /ranger_inputs/ and /ranger_outputs/ is NEVER part of the filenames here 
	return "ranger_outputs/", ranger_outputs

#march2018
def undo_replace_underscores(ranger_outputs):
	"""
	opens all the ranger output files, replaces the XX0XX in each tipname 
	with an underscore like it was originally

	args: list of file names
	output: modifies the contents of those files to have correct underscore placement.
	"""
	for file in ranger_outputs:
		with open("ranger_outputs/"+file) as old:
			text = old.read()
		newtext = text.replace("XX0XX", "_")
		with open("ranger_outputs/"+file,"w") as new:
			new.write(newtext)
	return ranger_outputs

#march2018
def separate_bootstrap_trees(bootfile):
	""" 
	this takes in a bootstrap file, and outputs each individual tree as one string in a list
	
	args: 
		bootfile : contains the 100 bootstraps newline separated (as output from raxml)
	output:
		sep_boots_list : a list containing each bootstrap tree as a sting.
	"""
	list_of_boot_trees = []
	with open(bootfile) as openbootfile:
		for line in openbootfile:
			list_of_boot_trees.append(line.strip())
	assert list_of_boot_trees != [], "failed to populate bootstrap file: "+bootfile
	return list_of_boot_trees

#march2018
def ensure_comparable_tips(sep_boots_list, species_string):
	str_to_tips, overall_tiplist = String_to_Tips(sep_boots_list[0])
	newtip_info_dict = {}
	tip_newtip_tups = []
	#keep gi number in the gene_bootstraps to avoid identicle-tip-errors.

	for tip in overall_tiplist:
		newtip = re.sub("(.*?\|)([A-z_0-9]*)(\|gi.*)", "\\2", tip)
		taxinfo = re.sub("(.*?\|)([A-z_0-9]*)(\|gi.*)", "\\1", tip)
		ginum = re.sub("(.*?\|)([A-z_0-9]*)(\|gi#)(.*)", "\\4", tip)
		newtip_info_dict[newtip] = taxinfo
		#newtip should be just Genus and species, with XX0XX separating the two
		#followed by an underscore and the gi number.
		tip_newtip_tups.append([tip,replace_underscores(newtip)+"_"+ginum])

	#ensure the new tips are the exact same as in the species tree.
	str_to_tips_sp, overall_tiplist_sp = String_to_Tips(species_string)
	for tip in newtip_info_dict.keys():
		if tip in overall_tiplist_sp:
			pass
		else:
			print("tip mismatch")
			print(tip)
			print(newtip_info_dict.keys())
			print(overall_tiplist_sp)
			print("why does your gene tree have a tip that your species tree does not.")
			print("this is probably going to error RIP")
			print("please re-make your gene trees without the tip:")
			print(tip)
			raise SystemExit

	#actually do the replacing
	new_boots_list = []
	for tree in sep_boots_list:
		for tippair in tip_newtip_tups:
			tree = tree.replace(tippair[0],tippair[1])
		new_boots_list.append(tree)
		

	#remove the underscores, since rangerDTL errors with them
	return new_boots_list, species_string, newtip_info_dict

#TODO
def fix_tips(bootstrap_tree, species_tree, name_mapping):
	#this will maybe add taxonomic information back in?
	pass

#march2018
def replace_underscores(string):
	newstring = string.replace("_", "XX0XX")
	return newstring

#update march2018, from old_code
def String_to_Tips(string):
	"""
	given a .newick subtree as a string, return a list of the tips

	args
		string: should be a .newick string eg "((A,B),((C,D),E))"
	returns
		overall_tiplist: a list that looks like ["A", "B", "C", "D", "E"]
		str_to_tips: a dict of all subclades to their tiplists. {"((A,B),((C,D),E))":["A", "B", "C", "D", "E"],
								(A,B):[A,B] ...etc}
	"""
	depth = 0
	depth_to_index = {}
	index = 0
	str_to_tips = {}
	list_of_strings = []
	overall_tiplist = []
	for item in string:
		if item == "(":
			depth += 1
			depth_to_index[depth] = index
		if item == ")":
			list_of_strings.append(string[depth_to_index[depth]:index+1])
			depth -= 1
		index += 1
	#print(list_of_strings)
	for item in list_of_strings:
		tiplist = []
		item2 = item.replace("(", "")
		item2 = item2.replace(")", "")
		items2 = item2.split(",")
		#print(items)
		for tip in items2:
			tipsplits = tip.split(":")
			tip = tipsplits[0]
			tiplist.append(tip)
			if tip not in overall_tiplist:
				overall_tiplist.append(tip)
		str_to_tips[item] = tiplist
	#print(str_to_tips)
	#add single tips str mapping to singletip in list.
	#for item in overall_tiplist:
	#	str_to_tips[item] = [item]
	return str_to_tips, overall_tiplist


#march2018
def write_inputs(sep_boots_list, speciesfile, current_obj_name):
	"""
	this function will take in a list of bootstrap trees, a filename of a species tree,
	and create files to be fed into rangerDTL for each bootstrap-species matchup

	args:
		sep_boots_list: a list of bootstrapped gene trees as strings
		speciesfile : a filename of the species besttree to compare to
		current_obj_name : the name of the dataset currently being compared
	"""
	i = 0
	input_list = []
	#read the species file
	with open(speciesfile, "r") as sp:
		species_string = sp.read().strip()

	#fix the incompatible naming conventions
	sep_boots_list, species_string, newtip_info_dict = ensure_comparable_tips(sep_boots_list, species_string)

	#remove the underscores (XX0XX replace)
	species_string = replace_underscores(species_string)
	#go one by one through the boots
	for boot in sep_boots_list:
		#boot = replace_underscores(boot)
		i+=1
		new_name = str(i)+current_obj_name+".RangerIn"
		with open ("ranger_inputs/"+new_name,"w") as new:
				new.write("[&R]"+species_string+"\n")
				new.write("[&R]"+boot)
		input_list.append(new_name)
	assert input_list != [], "did not create ranger outputs for "+current_obj_name
	return input_list, species_string

#march2018
def run_programs(input_list, species_string, current_obj_name):
	"""
	this should call OptRoot.linux and Ranger-DTL.linux
	if you are on a mac, please edit this
	make sure your programs are in your path, or in the folder you are running from.

	args:
		input_list : a list of ranger_inputs
		species_string : literally just the species tree in string forma
		current_obj_name : for naming purposes
	outputs:
		there should now exist 100 ranger dtl output files.
	returns:
		a list of the ranger_dtl_output files

	"""
	os.system("mkdir ranger_outputs")
	ranger_outputs = []
	for inputfile in input_list:
		optroot_out = inputfile.split(".")[0]+".OptRoot"
		ranout_file = inputfile.split(".")[0]+".RangerOut"
		os.system("./OptRoot.linux -i ranger_inputs/"+inputfile+" -o ranger_inputs/"+optroot_out)
		input_dtl = parse_optroot_output(optroot_out, species_string)
		os.system("./Ranger-DTL.linux -i ranger_inputs/"+input_dtl+" -o ./ranger_outputs/"+ranout_file) 
		ranger_outputs.append(ranout_file)
		assert os.path.isfile("ranger_outputs/"+ranout_file), "no ranoutfile generated at"+ranout_file
	return ranger_outputs

#march2018
def parse_optroot_output(optroot_out, species_string):
	"""
	this opens the output file that optroot generates,
	takes the best (or first, if there are multiple best) re-rooted bootstrap tree
	and writes it to file, followed by the species tree, as input for rangerDTL
	args:
		optroot_out: a single optroot output file
		species_string : the species tree in string form
	output:
		a newly written file containing the optimally rooted bootstrap and the species tree
		in the correct format to be fed into rangerDTL
	returns:
		name of that newly written file
	"""
	with open("ranger_inputs/"+optroot_out) as infile:
		tree = ""
		for line in infile:
			if line[0] == "(":
				if tree == "":
					tree = line
				else:
					break
		assert tree != "", "did not find an optrooted tree in"+optroot_out
	outputfile, junk = optroot_out.split(".")
	outputfile = outputfile +"RangerOut.txt"
	with open("ranger_inputs/"+outputfile, "w") as outfile:
		outfile.write("[&R]"+species_string+"\n")
		outfile.write("[&R]"+tree)
	assert os.path.isfile("ranger_inputs/"+outputfile) is True, "did not generate "+outputfile+" in parse_optroot"
	return outputfile
