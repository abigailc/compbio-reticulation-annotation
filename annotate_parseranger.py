#module for annotate_nodes_main.py

#abigailc@leviathan april 1 2018


#####imports
import re
import operator

########################
#new functions

########################
#medium functions

#march 2018
def write_the_bootstrap_tree(list_of_annotated_strings, output_file):
    """
    this function writes every bootstrap tree annotated by this program to one file, newline seperated

    args:
        list_of_annotated_strings : a list of annotated bootstrap gene trees
        output_file : a file name to write to
    output:
        a file containing all 100 bootstraps annotated 
    """
    with open(output_file, "w") as out:
        for item in list_of_annotated_strings:
            out.write(item+"\n")
    print("wrote the annotated bootstraps gene tree to "+output_file)
    return output_file

#march 2018
def write_the_species_tree(annotated_species_tree, output_file):
    """
    this function writes the species tree to file

    args:
        annotated_species_tree : a string of annotated species tree in .newick format
        output_file : a file name to write to
    output:
        a file containing the annotated species tree 
    """
    with open(output_file, "w") as out:
        out.write(annotated_species_tree)

    print("wrote the annotated species besttree to "+output_file)
    return output_file    

#march 2018
def Format_Information_Bootstrap(str_tips_dict, svd_spec,svd_dup,svd_trans,svd_loss,svd_recip, svd_donor):
    """
    this takes in a bunch of list that map subtree -> 1 if that node is a transfer/loss/etc
    suboptimal but then i don't have to re-write the code used for the species tree

    args
        str_tips_dict (contains all subclades)
        several lists of subclades (associated with eg "transfers")
    outputs:
        a single dictionary containing all of the subtrees and what the nodes should be labeled as.    
    """
    format_dict = {}
    for clade_str in str_tips_dict:
        output = ""
        if clade_str in svd_recip:
            output=output+("Recip~")
        if clade_str in svd_donor:
            output=output+("Donor~")
        if clade_str in svd_spec:
            output=output+("Sp")
        if clade_str in svd_dup:
            output=output+("Dup")
        if clade_str in svd_trans:
            output=output+("Txfr")
        if clade_str in svd_loss:
            output=output+("Loss")
        format_dict[clade_str] = output
    return format_dict

#march 2018
def Format_Information_SpeciesTree(str_tips_dict, master_speciestree):

    """
    this takes in a bunch of list that map subtree -> 1 if that node is a transfer/loss/etc
    suboptimal but then i don't have to re-write the code used for the species tree

    args
        str_tips_dict (contains all subclades)
        master_speciestree dictionary of cladestr:
        ["TransRecip", "TransDonor", "Duplication", "Loss", "Speciation", "Total"]

    outputs:
        a single dictionary containing all of the subtrees and what the nodes should be labeled as.    
    """
    format_dict = {}
    for clade in str_tips_dict:
        tips = str_tips_dict[clade]
        d = sorted(tips)
        string_cladetips  = " ".join(d)
        output = ""
        if master_speciestree[string_cladetips][0] != 0:
            output=output+("TxR"+str(master_speciestree[string_cladetips][0])+"~")
        if master_speciestree[string_cladetips][1] != 0:
            output=output+("TxD"+str(master_speciestree[string_cladetips][1])+"~")
        if master_speciestree[string_cladetips][2] != 0:
            output=output+("Dup"+str(master_speciestree[string_cladetips][2])+"~")
        if master_speciestree[string_cladetips][3] != 0:
            output=output+("Loss"+str(master_speciestree[string_cladetips][3])+"~")
        if master_speciestree[string_cladetips][4] != 0:
            output=output+("Sp"+str(master_speciestree[string_cladetips][4]))
        #this is currently always 0 because we aren't tracking overall bootstraps (oops)
        if master_speciestree[string_cladetips][5] != 0:
            output=output+("bs"+str(master_speciestree[string_cladetips][5]))

        format_dict[clade] = output
    return format_dict

#march2018
def Master_StrValDict_GeneTree(polarity_recip,polarity_donor,transfers,duplications,speciations, node_overall):
    """
    correlates a specific node (as defined by its tips) with the number of times it is IDed as various things
    args:
        str_tips_dict:dictionary of subclade:tips
        other dicts: dictionaries of tips(as a string space seperated and sorted):[value,YNSTRING]
    returns:
        one dictionary of tips(as space seperated sorted string):value_trans_recip, value_trans_donor, represents_trans, duplication, speciation, overall}
        where value is defined as "how many bootstraps this clade was identified as this thing in"
    """
    #get all potential subclade (from across all bootstrapped trees)
    overall_str_tips_dict = {}
    uber_strclade_list = list(transfers.keys())+list(speciations.keys())+list(duplications.keys())
    #for each subclade, tablulate # of transfers etc
    master_speciestree = {}
    for clade in uber_strclade_list:
        master_speciestree[clade] = [0,0,0,0,0,0]
    master_speciestree["clade"] = ["TransRecip", "TransDonor", "Transfer", "Duplication", "Speciation", "Total"]
    #populate the dict ints. order is trans_recip, trans_donor, trans_id, dup, speciation, total
    #for now only donor/recip are going to be mapped, i think.
    for str_clade in master_speciestree:
        if str_clade in polarity_recip:
            master_speciestree[str_clade][0] = master_speciestree[str_clade][0]+polarity_recip[str_clade][0]        
        if str_clade in polarity_donor:
            master_speciestree[str_clade][1] = master_speciestree[str_clade][1]+polarity_donor[str_clade][0]
        if str_clade in transfers:
            master_speciestree[str_clade][2] = master_speciestree[str_clade][2]+transfers[str_clade][0]
        if str_clade in duplications:
            master_speciestree[str_clade][3] = master_speciestree[str_clade][3]+duplications[str_clade][0]
        if str_clade in speciations:
            master_speciestree[str_clade][4] = master_speciestree[str_clade][4]+speciations[str_clade][0]
        if str_clade in node_overall:
            master_speciestree[str_clade][5] = master_speciestree[str_clade][5]+node_overall[str_clade][0]
    return master_speciestree

#march2018
def Master_StrValDict_SpeciesTree(str_species_tree, trans_recip, trans_donor, dup, loss, spec, overall):
    """
    correlates a specific node (as defined by its tips) with the number of times it is IDed as various things
    args:
        str_tips_dict:dictionary of subclade:tips
        other dicts: dictionaries of tips(as a string space seperated and sorted):[value,YNSTRING]
    returns:
        one dictionary of tips(as space seperated sorted string):value_trans_recip, value_trans_donor, value_etc...
        where value is defined as "how many bootstraps this clade was identified as this thing in"
    """
    str_tips_dict, overall_tiplist = String_to_Tips(str_species_tree)
    master_speciestree = {}
    for clade in str_tips_dict:
        tips = str_tips_dict[clade]
        d = sorted(tips)
        string_cladetips  = " ".join(d)
        master_speciestree[string_cladetips] = [0,0,0,0,0,0]

    master_speciestree["clade"] = ["TransRecip", "TransDonor", "Duplication", "Loss", "Speciation", "Total"]
    #populate the dict ints. order is trans_recip, trans_donor, dup, loss, speciation, total
    #for now only donor/recip are going to be mapped, i think.
    for str_clade in master_speciestree:
        if str_clade in trans_recip:
            master_speciestree[str_clade][0] = master_speciestree[str_clade][0]+trans_recip[str_clade][0]        
        if str_clade in trans_donor:
            master_speciestree[str_clade][1] = master_speciestree[str_clade][1]+trans_donor[str_clade][0]
        if str_clade in dup:
            master_speciestree[str_clade][2] = master_speciestree[str_clade][2]+dup[str_clade][0]
        if str_clade in loss:
            master_speciestree[str_clade][3] = master_speciestree[str_clade][3]+loss[str_clade][0]
        if str_clade in spec:
            master_speciestree[str_clade][4] = master_speciestree[str_clade][4]+spec[str_clade][0]
        if str_clade in overall:
            master_speciestree[str_clade][5] = master_speciestree[str_clade][5]+overall[str_clade][0]
    return master_speciestree

#######################

#larger functions

#march 2018
def annotate_the_bootstrap_tree(str_bs_tree, gene_tree_transfer_recipient, gene_tree_transfer_donor, speciation_res_ge, transfer_res_ge, loss_res_ge, duplicate_res_ge):
    """
    this should take in all the results of parsing the rangerDTL output for one single bootstrap
    and the string containing that bootstrap tree
    and it should return that tree annotated with what each node was found to represent when compared with the species tree
    args:
        str_bs_tree : a string containing the newick tree of this specific bootstrap
        all others : lists of the clades (in stringtips format) that contain transferrecips/donor/speciation/transferID/loss/duplicate/
    output:
        a string containing the newick tree with annotations as node-labels.
    """
    str_bs_tree = remove_node_labels(str_bs_tree)
    
    str_tips_dict, overall_tiplist = String_to_Tips(str_bs_tree)
   
    #results_of_transfer_polarity_analysis
    #gene_tree_transfer_recipient, gene_tree_transfer_donor
    #directly reading the ranger output file, correlating with subtrees via tips
    #svd is a dict {clade_as_string:value}
    #speciation_res_ge is a list of tips_in_a_clade_spacesep_sorted
    tips_to_str_dict = generate_tips_to_str(str_tips_dict)

    svd_spec = stringtips_to_cladestr(speciation_res_ge, tips_to_str_dict)
    svd_dup = stringtips_to_cladestr(duplicate_res_ge, tips_to_str_dict)
    svd_trans = stringtips_to_cladestr(transfer_res_ge, tips_to_str_dict)
    svd_loss = stringtips_to_cladestr(loss_res_ge, tips_to_str_dict)
    svd_recip = stringtips_to_cladestr(gene_tree_transfer_recipient, tips_to_str_dict)
    svd_donor = stringtips_to_cladestr(gene_tree_transfer_donor, tips_to_str_dict)

    #these are all just lists of the subclades
    dict_to_add_to_tree = Format_Information_Bootstrap(str_tips_dict, svd_spec,svd_dup,svd_trans,svd_loss,svd_recip, svd_donor)
    
    annotated_string = AddValToTree(str_bs_tree, dict_to_add_to_tree)
    
    return annotated_string
    
#march 2018
def annotate_the_species_tree(str_species_tree, master_speciestree):
    """
    This function creates an annotated string of the species tree, where each node is now
    labeled with information regarding in how many bootstraps it was found as a transfer
    recip/donor or as a loss/specition/duplication, as well as the overall bootstrap support.
    args:
        str_species_tree : a newick file of the species tree to annotate
        master_speciestree : a dictionary created by MasterStrValDict_SpeciesTree of form {clade_in_stringtips_form: [VAL VAL VAL VAL VAL VAL VAL]}
    output:
        annotated species tree string for visualization purposes
    """        
    #parse the species tree and find all tips, as well as a dictionary mapping clade-as-a-string to the-tips-the-clade-contains-as-a-list
    str_tips_dict, overall_tiplist = String_to_Tips(str_species_tree)
    #format your labels
    dict_to_add_to_tree = Format_Information_SpeciesTree(str_tips_dict, master_speciestree)
    #add the labels to the tree
    annotated_string = AddValToTree(str_species_tree, dict_to_add_to_tree)
    return annotated_string

###########################
#small functions

def remove_node_labels(string):
    """
    removes node labels as are added by rangerDTL
    """
    nolabel = re.sub("\)m\d*", ")", string)
    return nolabel


#march 2018
def generate_tips_to_str(str_tips_dict):
    """
    str tips dict is of form {clade of str: [list of tips]}, but 
    we want {space-sep-list-of-tips:clade of str}
    
    args:
        a dictionary of {clade_in_newick_form: [list_of_tips]}
    output:
        a dictionary of {clade_in_stringtips_form: clade_in_newick_form}
    """
    tips_to_str_dict = {}
    for item in str_tips_dict:
        sorted_list_of_tips = sorted(str_tips_dict[item])
        str_list_of_tips = " ".join(sorted_list_of_tips)
        tips_to_str_dict[str_list_of_tips] = item
    return tips_to_str_dict

#march2018
def stringtips_to_cladestr(list_stringtips, stringtips_to_cladestr_dict):
    """
    all this does is take a list of space-seperated-sorted-tips and a dictionary mapping tips:clade
    and return a list of the clades that match the given list of space-sep-tips
    args:
        list_stringtips: a list of clades in stringtip format
        stringtips to cladestr dictionary
    outputs:
        list_cladestr: a list of clades in cladestr format
    """
    cladestr_list = []
    for item in list_stringtips:
        if item in stringtips_to_cladestr_dict:
            cladestr_list.append(stringtips_to_cladestr_dict[item])
    return cladestr_list

#march 2018
def AddValToTree(str_tree, str_val_dict):
    """
    this function should take an input .newick style string and a dictionary of subtree:label
    it will then iteratively replace the nested subtrees with labeled subtrees.

    args:
        str_tree is a string of a newick
        str val dict is {substring of newick:value to insert}
    output:
        a string of a labeled newick
    """
    #solved the error with the ordering of the replacement strings --- needs to go insiude out.
    #dict is inharentlly unordered. lets make it ordered then. i can count commas
    list_of_items = []
    for item in str_val_dict:
        list_of_items.append(item)
    #needed to replace in correct inside-out order to preserve string recognition
    comma_list = Order_List_By_Commas(list_of_items)
    #skip anything with zero commas from the list --- to avoid mapping NodeLabels onto terminal tips
    for item in comma_list:
        if "," in item:
            str_tree = str_tree.replace(item, item+str_val_dict[item])
    return str_tree

#march 2018
def Order_List_By_Commas(list_input):
    """
    this takes in a list of strings, counts how many commas are in each
    and then reverse-sorts them by # of commas
    
    args:
        list_input: a list of strings to be ordered by # commas
    output:
        an ordered list of strings, where [0] has the most commas
    """
    comma_dict = {}
    for item in list_input:
        commanum = 0
        for character in item:
            if character == ",":
                commanum +=1
        comma_dict[item] = commanum
    sorted_x = sorted(comma_dict.items(), key=operator.itemgetter(1), reverse=True)
    #sorted_x will be a list of sorted tuples, those with MOST commas in front.
    #now turn it into a list
    comma_list = []
    for item in sorted_x:
        comma_list.append(item[0])
    assert comma_list != [], "did not populate comma list"
    return comma_list

#march 2018
def print_master_dict(input_master, output_file):
    """
    this nicely writes the master dict to a tab-seperated file

    input: a dictionary of clades:[alue_trans_recip, value_trans_donor, v_dup, v_loss, v_spec, v_all]
    output: a tab-seperated file
    """
    with open(output_file, "w") as new:
        new.write("clade\t"+str(input_master["clade"][0])+"\t"+str(input_master["clade"][1])+"\t"+str(input_master["clade"][2])+"\t"+str(input_master["clade"][3])+"\t"+str(input_master["clade"][4])+"\t"+str(input_master["clade"][5])+"\n")
        for k, v in input_master.items():
            new.write(k+"\t"+str(v[0])+"\t"+str(v[1])+"\t"+str(v[2])+"\t"+str(v[3])+"\t"+str(v[4])+"\t"+str(v[5])+"\n")
    print("wrote master dict to "+output_file)
    #writes the nicely formatted dictionary

#edit march 19 2018
def determine_polarity_on_the_gene_tree(results_polarity):
    """
    this function takes in a list of tuples.
    first entry is a string of a newick subclade that was IDed as a transfer
    second entry is a string of the tips identified as the donor clade (from the species tree) 
    
    it compares each side of the deepest split of the gene_tree_subclade
    to the species_tree_donor_tips
    to identify which side of the deepest split is the donor, and which the recipient.

    it outputs a list of donor clades and recipient clades.
    these can then be annotated (on the branch) as D or R

    """

    gene_tree_transfer_recipient = []
    gene_tree_transfer_donor = []
    errors_polarity = 0
    #iterate through each transfer that was found
    for transfer_found in results_polarity:
        transferclade_gene_string = transfer_found[0]
        #these tips are from the RECIPIENT clade
        donor_species_tips = transfer_found[1]
        #l stands for list (is a list of tips)
        #s stands for string (is a full clade)
        #1 is the first side of the split, 2 is the second side of the split
        l1,l2,s1,s2 = DeepestSplit(transferclade_gene_string, True, True)    
        l1_val = 0
        l2_val = 0
        donor_species_tips_list = donor_species_tips.split(" ")
        for item in donor_species_tips_list:
            for l1_item in l1:
                if l1_item == item:
                   l1_val +=1
            for l2_item in l2:
                if l2_item == item:
                    l2_val +=1
            #larger valued side contains more tips of what was called the recipient in the species tree
            #whichever side is largest should be saved as the recipient lineage
            #whichever side is smaller as the donor (smaller should = 0)
        #print(l1_val)
        #print(l2_val)
        if l1_val > l2_val:
            #add the string-version of tips space seperated to a list, just like the other lists eg transfer_recip+res+sp
            #sort the list first so that it always matches preexisting clades
            strtips, alltips = String_to_Tips(s1)
            strtips2, alltips2 = String_to_Tips(s2)
            tips_list_l1 = sorted(alltips)
            str_of_l1 = ' '.join(tips_list_l1)
            tips_list_l2 = sorted(alltips2)
            str_of_l2 = ' '.join(tips_list_l2)
            gene_tree_transfer_recipient.append(str_of_l1)
            gene_tree_transfer_donor.append(str_of_l2)
        elif l2_val > l1_val:
            strtips, alltips = String_to_Tips(s1)
            strtips2, alltips2 = String_to_Tips(s2)
            tips_list_l1 = sorted(alltips)
            str_of_l1 = ' '.join(tips_list_l1)
            tips_list_l2 = sorted(alltips2)
            str_of_l2 = ' '.join(tips_list_l2)
            gene_tree_transfer_recipient.append(str_of_l2)
            gene_tree_transfer_donor.append(str_of_l1)
        else:
            errors_polarity += 1
            #print("current equal errors:"+str(errors_polarity))
            continue
    #we need to append the list WITH the gi numbers associated.
    #print(gene_tree_transfer_recipient)
    return gene_tree_transfer_recipient, gene_tree_transfer_donor

#update of old code March 2018
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
    #    str_to_tips[item] = [item]
    return str_to_tips, overall_tiplist

#################################
#main function

#update march 2018
def parse_ranger_clade(rangerfiles, where_to_write_annotated_bootstraps, where_to_write_annotated_speciestrees, str_species_tree):
    """
    this takes in a list of rangerDTL output files, reads them,
    and returns tab-seperated files listing occurances, a visual of the species tree, and a visual of bs file

    args:
        rangerfiles: a list of all the files to run on (with path appended to each please)
        where_to_write_annotated_bootstraps (this should be a path)
        where_to_write_annotated_speciestrees (this should be a path)
        str_species_trees: a string of the species besttree.
            #TODO see if this causes error opening in FigTree. may need to use ranger-out-species tree.
    output:
        writes 4 files:
        -clade: what it was IDed as (species tree)
        -clade: what is was IDed as (gene tree)
        -annotated bootstraps newick file
        -annotated species_tree newick file
    """
    #all this section does is make a bunch of empty dictionaries that will be used by parse_ranger_clades
    
    #this tracks actual gene tree node that received the transfer, instead of which represent the transfer
    node_polarity = {}
    node_polarity_recip = {}
    node_polarity_donor = {}
    
    #this counts number of times node represents a transfer (in the case of species tree, we count donor and recip seperately)
    node_transfers_sp_recip = {}
    node_transfers_sp_donor = {}
    node_transfers_ge = {}
    #this counts number of times node represents a duplication
    node_duplication_sp = {}
    node_duplication_ge = {}
    #this counts number of times node represents a speciation
    node_speciation_sp = {}
    node_speciation_ge = {}
    #NOTWORKING - cannot infer from current datasets
    #this counts number of times node represents a loss 
    node_loss_sp = {}
    node_loss_ge = {}
    #this counts how many times node is ever found
    node_overall_sp = {}
    node_overall_ge = {}

    """
    this will one-by-one open a ranger output file, parse it, save an annotated optrooted bootstrap tree, and update the dictionaries.
    """
    #bootstrap number
    bs_num = 1
    #list of annotated bootstrap_gene_tree strings in order
    annotated_bootstrap_list = []
    #will iterate through all available bootstraps
    print("one_by_one")
    for singlebs in rangerfiles:
        #first, we parse a single bootstrap file.
        results_sp, results_ge, results_polarity, str_bs_tree = Parse_Single_Ranger_Output_BOTHWAYS(singlebs)
        #parse results with respect to species tree
        speciation_res_sp = results_sp[0]
        transfer_recip_res_sp = results_sp[1]
        transfer_donor_res_sp = results_sp[2]
        loss_res_sp = results_sp[3]
        duplicate_res_sp = results_sp[4]
        overall_res_sp = results_sp[5]

        #parse results with respect to gene tree
        speciation_res_ge = results_ge[0]
        transfer_res_ge = results_ge[1]
        loss_res_ge = results_ge[2]
        duplicate_res_ge = results_ge[3]
        overall_res_ge = results_ge[4]

        #deal with polarity in the gene tree
        gene_tree_transfer_recipient, gene_tree_transfer_donor = determine_polarity_on_the_gene_tree(results_polarity)
        #these are just lists of clade-strs

        #write the results of this bootstrap onto it in newick format and save
        #each node should have SP TX or DU
        #branches below a TX should have R or D if it is determinable
        annbs = annotate_the_bootstrap_tree(str_bs_tree, gene_tree_transfer_recipient, gene_tree_transfer_donor, speciation_res_ge, transfer_res_ge, loss_res_ge, duplicate_res_ge)
        annotated_bootstrap_list.append(annbs)

        #now, update the overall species_tree dictionaries
        node_transfers_sp_recip = dictionary_check_bipart_BOTHWAYS(node_transfers_sp_recip, transfer_recip_res_sp)
        node_transfers_sp_donor  = dictionary_check_bipart_BOTHWAYS(node_transfers_sp_donor, transfer_donor_res_sp)
        node_loss_sp  = dictionary_check_bipart_BOTHWAYS(node_loss_sp, loss_res_sp)
        node_duplication_sp  = dictionary_check_bipart_BOTHWAYS(node_duplication_sp, duplicate_res_sp)
        node_speciation_sp  = dictionary_check_bipart_BOTHWAYS(node_speciation_sp, speciation_res_sp)
        node_overall_sp  = dictionary_check_bipart_BOTHWAYS(node_overall_sp, overall_res_sp)

        #and the overall gene tree dicitonaries
        node_transfers_ge  = dictionary_check_bipart_BOTHWAYS(node_transfers_ge, transfer_res_ge)
        node_loss_ge  = dictionary_check_bipart_BOTHWAYS(node_loss_ge, loss_res_ge)
        node_duplication_ge  = dictionary_check_bipart_BOTHWAYS(node_duplication_ge, duplicate_res_ge)
        node_speciation_ge  = dictionary_check_bipart_BOTHWAYS(node_speciation_ge, speciation_res_ge)
        node_overall_ge  = dictionary_check_bipart_BOTHWAYS(node_overall_ge, overall_res_ge)

        #and the polarity determination, maybe not relevent for species tree
        node_polarity_donor = dictionary_check_bipart_BOTHWAYS(node_polarity_donor, gene_tree_transfer_donor)
        node_polarity_recip = dictionary_check_bipart_BOTHWAYS(node_polarity_recip, gene_tree_transfer_recipient)
    #now we have finished looping, so write to file all the bootstraps
    write_the_bootstrap_tree(annotated_bootstrap_list, where_to_write_annotated_bootstraps+"_ANNboots.newick")
    #and all the clade_specific_info 
    master_genetree = Master_StrValDict_GeneTree(node_polarity_recip,node_polarity_donor,node_transfers_ge,node_duplication_ge,node_speciation_ge,node_overall_ge)
    print_master_dict(master_genetree, where_to_write_annotated_bootstraps+"_master_gene_spreadsheet")
    #and all the species_tree_info
    master_speciestree = Master_StrValDict_SpeciesTree(str_species_tree, node_transfers_sp_recip, node_transfers_sp_donor, node_duplication_sp, node_loss_sp, node_speciation_sp, node_overall_sp)
    print_master_dict(master_speciestree, where_to_write_annotated_speciestrees+"_master_species_spreadsheet")
    annotated_species_tree = annotate_the_species_tree(str_species_tree,master_speciestree)
    write_the_species_tree(annotated_species_tree, where_to_write_annotated_speciestrees+"_ANNspecies.newick")
    #now write the bipartitions file! haha its still a mess

    #the old way of doing this
    bipartitions_table = where_to_write_annotated_bootstraps+"_bipartitions.txt"
    with open(bipartitions_table, "w") as bipart:
        bipart.write("BASED ON SPECIES TREE\n")
        bipart.write("Transfer_Recipients_Sp\n")
        sort_write_the_dict(bipart, node_transfers_sp_recip)
        bipart.write("Transfer_Donors_Sp\n")
        sort_write_the_dict(bipart, node_transfers_sp_donor)
        bipart.write("Loss_Sp\n")
        sort_write_the_dict(bipart, node_loss_sp)
        bipart.write("Duplicate_Sp\n")
        sort_write_the_dict(bipart, node_duplication_sp)
        bipart.write("Speciation_Sp\n")
        sort_write_the_dict(bipart, node_speciation_sp)
        bipart.write("Overall_Sp\n")
        sort_write_the_dict(bipart, node_overall_sp)
        bipart.write("BASED ON GENE TREE\n")
        bipart.write("Transfers_by_recip_clade\n")
        sort_write_the_dict(bipart, node_polarity_recip)
        bipart.write("Transfers_by_donor_clade\n")
        sort_write_the_dict(bipart, node_polarity_donor)
        bipart.write("Transfer_Ge\n")
        sort_write_the_dict(bipart, node_transfers_ge)
        bipart.write("Loss_Ge\n")
        sort_write_the_dict(bipart, node_loss_ge)
        bipart.write("Duplicate_Ge\n")
        sort_write_the_dict(bipart, node_duplication_ge)
        bipart.write("Speciation_Ge\n")
        sort_write_the_dict(bipart, node_speciation_ge)
        bipart.write("Overall_Ge\n")
        sort_write_the_dict(bipart, node_overall_ge)

    
    #print("completed parsing and wrote to: "+bipartitions_table)

    return True


##################################
#OLD FUNCTIONS
#GENERALLY UNTOUCHED
#here be dragons

#update june 29
def dictionary_check_bipart_BOTHWAYS(node_dict, list_res):
    #what does node1 look like though?
    #space-seperated tips sorted alphabetically!
    #node_dict is something like node_loss_dict = [node1:1, node2:4]
    #list_res is something like loss_list = [node1, node2, node4]
    #tracks each bootstrap y/n for each node's presence. also have an overall value (#of times node was yes)
    #0 = total
    #1 = YNNNNYYYYYNNN
    bs_num = 0
    for each_node in node_dict:
        #track how many N already should exist for use later.
        bs_num = len(node_dict[each_node][1])
        err = "yes"
        #look to see if we should add a Y or N to each member of the clade dictionary.
        for result_node in list_res:
            #if transfer exists in the dictionary and in current bootstrap, add a Y and a number.
            if result_node == each_node:
                node_dict[result_node][0] += 1
                node_dict[result_node][1].append("Y")
                err = "no"
        #if transfer in dictionary but not in current bootstrap, add a N and no number.
        if err == "yes":
            node_dict[each_node][1].append("N")        
    #if the transfer is NOT in the dictionary, add it with the appropriate number of NOs.
    for result_node in list_res:
        if result_node in node_dict:
            pass
        else:
            inner_list = []
            for item in range(bs_num):
                inner_list.append("N")
            inner_list.append("Y")
            node_dict[result_node] = [1,inner_list]
    return node_dict

#copied from do_subsampling jul 28
def DeepestSplit(tree, mfix = False, gifix = False):
    #provided subtree will be text. like
    #(((((((((((Elephantulus_edwardii:0.02656660969,Orycteropus_afer:0.006886559777):0.002186555754,Trichechus_manatus:0.0207867919)89:0.0042892984,(Erinaceus_europaeus:0.0150776695,Sus_scrofa:0.01871033399)44:0.002137367876)98:0.001261277757,Equus_asinus:0.006995816459)11:0.0003479426504,((Pteropus_vampyrus:0.01260656735,(Pan_paniscus:0.007271418717,Manis_javanica:0.01036292681)13:0.001322734885)9:0.001552612205,(Dasypus_novemcinctus:0.009902206268,Tupaia_chinensis:0.02113855241)30:0.001364556134)12:0.001840714523)4:0.001979025907,Galeopterus_variegatus:0.00671926141)25:0.004692110764,Ochotona_princeps:0.02065394194)71:0.005193088649,Mus_musculus:0.01442576505)60:0.02036347729,(Monodelphis_domestica:0.007867166076,Sarcophilus_harrisii:0.01059171724)100:0.02083374525)100:0.007086087413,Ornithorhynchus_anatinus:0.03185159824)65:0,Gekko_japonicus:0.04818786247)Root;
    #returns: ["tipa1","tipa2", "tipa3"],["tipb1","tipb2","tipb3"]
    # where the first list represents strings in the first clade made by deepest split, and second list represents strings in the other clade.
    
    #if mfix is True:
    #    if tree[-1] == ",":
    #        tree = tree[:-1]+")"
    p = False
    indent = 0
    first = ""
    switch = False
    last = ""
    if "Root;" in tree:
        tree.replace("Root;", ";")
    for character in tree:
        if character == "(":
            indent +=1
        if character == ",":
            if indent == 1:
                switch = True
            else:
                indent = indent - 1
        if switch is False:
            first = first+character
        else:
            last = last+character
    first = first[1:] #cuts the initial paren
    last = last[1:-2] #cuts seperating comma and final paren +;
    #now we have the two deepest subtrees.
    #lets extract the tip names.
    #CONSIDER that you need to exclude any NODE LABELS.
    first_edit = re.sub("(:[^,]*)", "", first)
    first_edit = re.sub("[\(\)]", "", first_edit)
    first_edit_list = first_edit.split(",")
    last_edit = re.sub("(:[^,]*)", "", last)
    last_edit = re.sub("[\(\)]", "", last_edit)
    last_edit_list = last_edit.split(",")
    #returns list of tip_names in side1 of deepest split, then side2.
    # rets [a,b,c], [d,e,f], "(a,(b,c)", "(d,e),f)"

    #short term fix for the m - node values
    if mfix is True:
        new_first_edit_list = []
        new_last_edit_list = []
        for item in first_edit_list:
            new_first_edit = re.sub("([^_]*)(_)([^m]*)(m)(.*)", "\\1\\2\\3", item)
            new_first_edit_list.append(new_first_edit)
        for item in last_edit_list:
            new_last_edit = re.sub("([^_]*)(_)([^m]*)(m)(.*)", "\\1\\2\\3", item)
            new_last_edit_list.append(new_last_edit)
        first_edit_list = new_first_edit_list
        last_edit_list = new_last_edit_list
    if gifix is True:
        new_first_edit_list = []
        new_last_edit_list = []
        for item in first_edit_list:
            new_first_edit = re.sub("(.*)(_)(\d*)", "\\1", item)
            new_first_edit_list.append(new_first_edit)
        for item in last_edit_list:
            new_last_edit = re.sub("(.*)(_)(\d*)", "\\1", item)
            new_last_edit_list.append(new_last_edit)
        first_edit_list = new_first_edit_list
        last_edit_list = new_last_edit_list
    first = first+";"
    last = last+";"
    return first_edit_list, last_edit_list, first, last

#copied from do_subsampling jul 28
def Node_Newick(tree):
    #provided subtree will be text. like ((A,B)x1,(C,D)x2,E)x3;
    node_to_newick_dict = {}
    indent = 0
    first = ""
    name = False
    track = False
    tracking = {}
    tracking[0]=0
    last = ""
    if "Root;" in tree:
        tree.replace("Root;", ";")
    for i in range(len(tree)):
        character = tree[i]
        if character == "(":
            indent +=1
            tracking[indent] = i
            name = False
            track = True
        elif character == ")":
            if name is True:
                node_to_newick_dict[currentname] = currentstring
                currentstring = tree[tracking[indent]:i+1]
                currentname = ""
                track = False
                name = True
            else:
                currentstring = tree[tracking[indent]:i+1]
                track = True
                name = True
                currentname = ""
            indent = indent - 1
        elif character == ",":
            if name is True:
                node_to_newick_dict[currentname] = currentstring
                track = False
                name = False
            else:
                track = False
                name = False
        elif character == ";":
            if name is True:
                node_to_newick_dict[currentname] = currentstring
                currentstring = tree[tracking[indent]:i+1]
                currentname = ""
                track = False
        else:
            if name is True:
                currentname = currentname+character
            last = last+character
    #return should be like : {x1:(A,B), x2:(C,D), x3:((A,B)x1,(C,D)x2,E)}
    return node_to_newick_dict

#update feb 2018
def sort_write_the_dict(openfile, diction):
    """
    not currently using
    """
    for i, value in sorted(diction.items(), key=lambda v: v[0], reverse=True):
        i_string = str(i)
        str1 = ''.join(value[1])
        openfile.write(i_string+"\t"+str(diction[i][0])+"\t"+str1+"\n")

#update june29 in do_subsampling
def Parse_Single_Ranger_Output_BOTHWAYS(rangerfile):

    #initial july
    res_JULY = []

    #initial species lists
    speciation_res_sp = []
    transfer_recip_res_sp = []
    transfer_donor_res_sp = []
    loss_res_sp = []
    duplicate_res_sp = []
    overall_res_sp = []

    #initial gene lists
    speciation_res_ge = []
    transfer_res_ge = []
    loss_res_ge = []
    duplicate_res_ge = []
    overall_res_ge = []

     #initialize lists
    results_sp = [speciation_res_sp, transfer_recip_res_sp, transfer_donor_res_sp, loss_res_sp, duplicate_res_sp, overall_res_sp]
    results_ge = [speciation_res_ge, transfer_res_ge, loss_res_ge, duplicate_res_ge, overall_res_ge]

    #set to recognize gene and species trees
    prep = False
    sp_prep = False
    with open (rangerfile) as ranger:
        for line in ranger:
            if "Gene Tree:" in line:
                prep = True
            elif "Species Tree:" in line:
                sp_prep = True
            elif "minimum reconciliation cost" in line:
                continue
            elif prep is True:
                #make the gene tree dict
                prep = False
                str_bs_tree = line.strip()
                gene_tree_dict = Newick_Tree_Nodes_To_Dict(line)
                node_newick_gene = Node_Newick(line)
            elif sp_prep is True:
                #make the species tree dict
                sp_prep = False
                str_sp_tree = line.strip()
                species_tree_dict = Newick_Tree_Nodes_To_Dict(line)
                node_newick_species = Node_Newick(line)

            elif "Speciation" in line:
                #store info one speciation
                #m3 = LCA[MeiothermusXX0XXcerbereus_654400501, MeiothermusXX0XXtimidus_517278376]: Speciation, Mapping --> n7
                node_finder = line.split(" ")
                node_ge = node_finder[0]
                node_sp = node_finder[-1]
                node_ge = node_ge.strip()
                node_sp = node_sp.strip()
                #convert and add gene version
                list_of_tips_ge = gene_tree_dict[node_ge]
                str_of_tips_ge = ' '.join(list_of_tips_ge)
                speciation_res_ge.append(str_of_tips_ge)
                #convert and add species version
                list_of_tips_sp = species_tree_dict[node_sp]
                str_of_tips_sp = ' '.join(list_of_tips_sp)
                speciation_res_sp.append(str_of_tips_sp)
            #
             # m8 = LCA[PyrobaculumXX0XXneutrophilum_501318754, PyrobaculumXX0XXoguniense_504114203]: Transfer, Mapping --> n11, Recipient --> n12
            elif "Transfer" in line:
                node_finder = line.split(" ")
                node_ge = node_finder[0]
                node_sp_recip = node_finder[-1]
                node_sp_donor = node_finder[-4]
                node_sp_donor = node_sp_donor.strip(",")
                node_ge = node_ge.strip()
                node_sp_donor = node_sp_donor.strip()
                node_sp_recip = node_sp_recip.strip()

                #convert and add gene version
                list_of_tips_ge = gene_tree_dict[node_ge]
                str_of_tips_ge = ' '.join(list_of_tips_ge)
                transfer_res_ge.append(str_of_tips_ge)
                #convert and add species version
                list_of_tips_sp = species_tree_dict[node_sp_donor]
                str_of_tips_sp = ' '.join(list_of_tips_sp)
                transfer_donor_res_sp.append(str_of_tips_sp)
                list_of_tips_sp = species_tree_dict[node_sp_recip]
                str_of_tips_sp = ' '.join(list_of_tips_sp)
                transfer_recip_res_sp.append(str_of_tips_sp)
        
                pair_JULY = (node_newick_gene[node_ge], str_of_tips_sp)
                #pair july is the newick of the specific bootstrap , the tips sorted and reattached
                #this should be ("(a,b),(c,d)", [c, d])
                res_JULY.append(pair_JULY)
        
            elif "Loss" in line:
                node_finder = line.split(" ")
                node_ge = node_finder[0]
                node_sp = node_finder[-1]
                node_ge = node_ge.strip()
                node_sp = node_sp.strip()
                #convert and add gene version
                list_of_tips_ge = gene_tree_dict[node_ge]
                str_of_tips_ge = ' '.join(list_of_tips_ge)
                loss_res_ge.append(str_of_tips_ge)
                #convert and add species version
                list_of_tips_sp = species_tree_dict[node_sp]
                str_of_tips_sp = ' '.join(list_of_tips_sp)
                loss_res_sp.append(str_of_tips_sp)

    #             do the third thing
            elif "Duplication" in line:
                node_finder = line.split(" ")
                node_ge = node_finder[0]
                node_sp = node_finder[-1]
                node_ge = node_ge.strip()
                node_sp = node_sp.strip()
                #convert and add gene version
                list_of_tips_ge = gene_tree_dict[node_ge]
                str_of_tips_ge = ' '.join(list_of_tips_ge)
                duplicate_res_ge.append(str_of_tips_ge)
                #convert and add species version
                list_of_tips_sp = species_tree_dict[node_sp]
                str_of_tips_sp = ' '.join(list_of_tips_sp)
                duplicate_res_sp.append(str_of_tips_sp)
    #            
            elif "Leaf Node" in line:
                continue
            else:
                #this should ignore the uh initial SPecies Tree": and the species tree itseld and any \n lines
                continue
    return results_sp, results_ge, res_JULY, str_bs_tree


#OLD
def Newick_Tree_Nodes_To_Dict(genetree_string):
    #genetree_string = (((A,B)n1,(C,D)n2)n3),E)n4
    node_to_tips_dict = {}
    finished = False
    overall_tips_list = []

    #get node name: tips
    while finished is False:   
        #get the name of the node, eg N1
        node_name = re.sub ("(.*)(\()([^\(|\)]*)(\))([a-z][0-9]*)([^\)]*)(.*)", "\\5", genetree_string)
        #get stuff inside of outermost parentheses
        node_tips_string = re.sub ("(.*)(\()([^\(|\)]*)(\))([a-z][0-9]*)([^\)]*)(.*)", "\\3", genetree_string)
        #should give the first inner clade it sees -- in this case, (A,B)
        #the tips will be seperated by commas
        genetree_string = re.sub ("(.*)(\()([^\(|\)]*)(\))([a-z][0-9]*)([^\)]*)(.*)", "\\1\\3\\6\\7", genetree_string)
        node_tips_string = node_tips_string.strip()
        node_name = node_name.strip()
        tips_list = node_tips_string.split(",")
        #sort so that m1 = a,b is the same as m1 = b,a
        #m1 = [a,b]
        tips_list = sorted(tips_list)
        for item in tips_list:
            if item not in overall_tips_list:
                overall_tips_list.append(item)
        node_to_tips_dict[node_name] = tips_list
        it = 0
        for item in genetree_string:
            if item == "(":
                it += 1
        if it == 0:
            finished = True
    #print(node_to_tips_dict)
    #get node name: string of that clade
    #I think I also want to add add the single tips to the dictionary with node_name = tip and tips_list = tip.
    for tip in overall_tips_list:
        node_to_tips_dict[tip] = [tip]

    #print("node to dict")
    #print(node_to_tips_dict)
    return node_to_tips_dict

    #this dict will be used to determine which tips received the gene via a transfer even if the transfer is deeper in the tree.
