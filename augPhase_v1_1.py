'''
2013, 2014
by Hyun Uk Kim, Jae Yong Ryu and Kyu-Sang Hwang
'''

# import Model, Reaction, Metabolite classes in COBRA tool 
from cobra import Model, Reaction, Metabolite
from Bio import SeqIO
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.io.sbml import write_cobra_model_to_sbml_file

import copy
import os
import pickle
import re
import urllib2


def make_locusTag_geneID_nonBBH(gbkFile, fileType, inputFile, outputFile):
    fp2 = open(outputFile, "w")
    locusTag_geneID_dict = {}
    geneID_locusTag_dict = {}

#Reads GenBank file
    record = SeqIO.read(gbkFile, fileType)

    for feature in record.features:
        fp1 = open(inputFile,"r")
        if feature.type == 'CDS':
            if 'db_xref' in feature.qualifiers:

#Reads first non-BBH file    
		targetLocusTag = fp1.readline()
		while targetLocusTag:
	            targetLocusTag = targetLocusTag.strip()

#Looks for identical ORF from non-BBH list
		    if feature.qualifiers['locus_tag'][0] == targetLocusTag:
			print "feature.qualifiers['locus_tag']:", feature.qualifiers['locus_tag']

#Standard .gbk has "GI" for the first db_xref and "GeneID" for the second db_xref.
#Following lines take whichever comes first.
		 	geneID = feature.qualifiers.get('db_xref')
		 	print "feature.qualifiers.get('db_xref'):", geneID
			geneID = geneID[0].split(':')
			geneID_type = geneID[0].strip()
			geneID = geneID[1].strip()
			print "geneID:", geneID 

#Writing output files
			fp2.write(str(feature.qualifiers['locus_tag'][0])+"\t")
			fp2.write(geneID+'\n')

#Saves data in Dictionary
			locusTag_geneID_dict[feature.qualifiers['locus_tag'][0]] = geneID
			geneID_locusTag_dict[geneID] = feature.qualifiers['locus_tag'][0] 
 		    targetLocusTag = fp1.readline()
    fp1.close()
    fp2.close()
    return locusTag_geneID_dict, geneID_locusTag_dict



#Converts ncbi_geneid(e.g. 944762) to species_geneid (e.g. b0031) in KEGG database by using KEGGAPI
#Input: ncbi_gene_id (string) (e.g., 1217641)
#Return: Locus_tag ID in KEGG database (e.g., sma:SAV_7535)
def get_species_locusTag(ncbi_geneid):
    try:
	url = "http://rest.kegg.jp/conv/genes/ncbi-gi:%s"%(ncbi_geneid)
    except:
	url = "http://rest.kegg.jp/conv/genes/ncbi-geneid:%s"%(ncbi_geneid)
#Open and read data for the results of query in url address
    data = urllib2.urlopen(url).read()
    sptlist = data.strip().split()
    print "sptlist:", sptlist
    species_locusTag = sptlist[1]
    print "species_locusTag:", species_locusTag
    return species_locusTag


#Retrieves EC numbers relevant to their  species_locusTag from KEGG
#Input: Species_locusTag in KEGG database (e.g., sma:SAV_7535)
#Output: A set of EC numbers in list form (e.g., ['4.1.3.6 ', '4.1.3.34]'])
def get_ECNumberList_from_locusTag(species_locusTag):
    url = "http://rest.kegg.jp/get/%s"%(species_locusTag)
    gene_info_text = urllib2.urlopen(url).read()
    split_text = gene_info_text.strip().split('\n')
    for line in split_text:
        sptlist = line.split()
        if sptlist[0] == 'ORTHOLOGY':
            ECNumbers = re.findall(r'[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+[\s\]]', line)
            return ECNumbers
    return []


#Retrives a list of reaction IDs using their EC numbers from KEGG
#Input: E.C number in string form (e.g., 4.1.3.6)
#Output: reactionID in list form (e.g., ['R00362'])
def get_rxnid_from_ECNumber(enzymeEC):
    url = "http://rest.kegg.jp/get/enzyme:%s"%(enzymeEC)
    ecinfo_text = urllib2.urlopen(url).read()

#Original line also extracted genes in other organisms: R50912; R50345 (NOT rxnid)
#The HTTP error was solved by putting "\\b" only at the end (not at the front) in order to also include reaction ID followed by "other" in KEGG
    rxnid_set = re.findall(r'\s+R[0-9]+[0-9]+[0-9]+[0-9]+[0-9]'+'\\b', ecinfo_text)
    rxnid_list = []    
    for each_set in rxnid_set:
        rxnid = re.findall('R[0-9]+[0-9]+[0-9]+[0-9]+[0-9]+', each_set)
        rxnid_list+=rxnid

#Removes redundancy
        rxnid_list = list(set(rxnid_list))
    return rxnid_list


#Get reaction information using its ID from KEGG
#Input: KEGG rxnid in string form (e.g., R00362)
#Output: Reaction information for 'Name', 'Definition', and 'Equation' as dictionary form {'NAME': 'citrate oxaloacetate-lyase', 'DEFINITION': Citrate <=> Acetate + Oxaloacetate, 'EQUATION': C00158 <=> C00033 + C00036}
def get_rxnInfo_from_rxnid(rxnid):  
    url = "http://rest.kegg.jp/get/rn:%s"%(rxnid)
    reaction_info_text = urllib2.urlopen(url).read()
    split_text = reaction_info_text.strip().split('\n')
    NAME = ''
    DEFINITION = ''
    EQUATION = ''
    ENZYME = ''
    PATHWAY = ''

    for line in split_text:
        sptlist = line.split()
        if sptlist[0].strip() == 'NAME':
            NAME = ' '.join(sptlist[1:])
        if sptlist[0].strip() == 'DEFINITION':
            DEFINITION = ' '.join(sptlist[1:])
        if sptlist[0].strip() == 'EQUATION':
            EQUATION = ' '.join(sptlist[1:])
        if sptlist[0].strip() == 'ENZYME':
            ENZYME = ' '.join(sptlist[1:])

#Considers only reactions mapped in pathways
        if sptlist[0].strip() == 'PATHWAY':        
            PATHWAY = ' '.join(sptlist[1:])
	    return {'NAME':NAME, 'DEFINITION':DEFINITION, 'EQUATION':EQUATION, 'ENZYME':ENZYME, 'PATHWAY':PATHWAY}

#Four nested function calling four functions above
#Input: locusTag_geneID_dict (e.g., { '1216669' : 'SAV_6272'})
#Output: GPR_dict (e.g. '1211478': {'2.1.1.132': {'R05149': {'DEFINITION': '2 S-Adenosyl-L-methionine + Precorrin 6Y <=> 2 S-Adenosyl-L-homocysteine + Precorrin 8X + CO2', 'EQUATION': '2 C00019 + C06319 <=> 2 C00021 + C06408 + C00011', 'NAME': 'S-Adenosyl-L-methionine:1-precorrin-6Y'}}})
def make_all_rxnInfo(locusTag_geneID_dict, outputFile1, outputFile2, outputFile3, outputFile4):
    fp1 = open(outputFile1,'w')
    fp2 = open(outputFile2,'w')
    fp3 = open(outputFile3,'w')
    fp4 = open(outputFile4,'w')
    rxnid_info_dict = {}
    rxnid_geneid_dict = {}
    rxnid_locusTag_dict = {}
 
    for locusTag in locusTag_geneID_dict.keys():
        ncbi_geneid = locusTag_geneID_dict[locusTag]
	species_locusTag = get_species_locusTag(ncbi_geneid) 
        ECnumbers = get_ECNumberList_from_locusTag(species_locusTag)
        
	if ECnumbers:
            for enzymeEC in ECnumbers:
                enzymeEC = enzymeEC.replace(']','')
                enzymeEC = enzymeEC.replace(' ','')
                
                rxnid_list = get_rxnid_from_ECNumber(enzymeEC)
                for rxnid in rxnid_list:
                    rxnid_info_dict[rxnid] = get_rxnInfo_from_rxnid(rxnid)

		    if rxnid not in rxnid_geneid_dict.keys():
		        rxnid_geneid_dict[rxnid] = [(ncbi_geneid)]
		        rxnid_locusTag_dict[rxnid] = [(locusTag)]

#Appends additional different genes to the same reaction ID
		    elif rxnid in rxnid_geneid_dict.keys():
			rxnid_geneid_dict[rxnid].append((ncbi_geneid))
			rxnid_locusTag_dict[rxnid].append((locusTag))
	
                    print locusTag, ncbi_geneid, rxnid, rxnid_info_dict[rxnid]
		    print "\n"
	else:
	    print locusTag, ": KEGG info NOT available"
	    print "\n"
	    print >>fp1, '%s' %(locusTag)
    print >>fp2, '%s' %(rxnid_info_dict)
    print >>fp3, '%s' %(rxnid_geneid_dict)
    print >>fp4, '%s' %(rxnid_locusTag_dict)

    fp1.close()
    fp2.close()
    fp3.close()
    fp4.close()
    return rxnid_info_dict, rxnid_geneid_dict, rxnid_locusTag_dict


def check_existing_rxns(kegg_mnxr_dict, templateModel_bigg_mnxr_dict, rxnid_info_dict, outputFile1, outputFile2):
    fp1 = open(outputFile1,'w')
    fp2 = open(outputFile2,'w')
    rxnid_to_add_list =[]

    for rxnid in rxnid_info_dict.keys():
#Considers only reactions mapped in pathways
	if rxnid in kegg_mnxr_dict.keys():
            kegg_mnxr = kegg_mnxr_dict[rxnid]

#Checks with reactions in the template model through MNXref
            if kegg_mnxr in templateModel_bigg_mnxr_dict.values():
                print >>fp1, "%s\t%s" % (rxnid, kegg_mnxr)
            elif rxnid not in rxnid_to_add_list:
                rxnid_to_add_list.append(rxnid)
#            	print 'ID of reaction-to-add: ', rxnid
              	print >>fp2, "%s" % (rxnid)
    fp1.close()
    fp2.close()
    return rxnid_to_add_list

#Output: MNXR for the reactions to add, converted from KEGG rxnid
def get_mnxr_using_kegg(rxnid_to_add_list, kegg_mnxr_dict, outputFile1):
    mnxr_to_add_list = []
    fp1 = open(outputFile1, "w")
    for rxnid in rxnid_to_add_list:
	if rxnid in kegg_mnxr_dict.keys():
	    mnxr_to_add_list.append(kegg_mnxr_dict[rxnid])
	    print >>fp1, "%s\t%s" %(rxnid, kegg_mnxr_dict[rxnid])
    fp1.close()
    return mnxr_to_add_list

#Creating: e.g., {'R03232': {'f1p': -1.0, 'C04261': 1.0, 'fru': 1.0, 'C00615': -1.0}}
#Metabolites are presented primarily with bigg, otherwise with KEGG
def extract_rxn_mnxm_coeff(mnxr_to_add_list, mnxr_rxn_dict, mnxm_bigg_compound_dict, mnxm_kegg_compound_dict, mnxr_kegg_dict, outputFile1, outputFile2):
    fp1 = open(outputFile1, "w")
    fp2 = open(outputFile2, "w")
    rxnid_mnxm_coeff_dict = {}
    mnxm_coeff_dict = {}

    for mnxr in mnxr_to_add_list:
	unparsed_equation = mnxr_rxn_dict[mnxr]
	print unparsed_equation

#"substrates" and "products" contain stoichiometric coeff of each compound
	sptReaction = unparsed_equation.split('=')
	substrates = sptReaction[0].strip()
	products = sptReaction[1].strip()

#Discards polymerization reactions with undefinite coeff's
#e.g., 1 MNXM9 + (n+2) MNXM90033 = 1 MNXM5617 + (n) MNXM90033
	if '(' not in substrates and '(' not in products:
#Creating: e.g., {bigg compoundID:(-1)coeff}, {kegg compoundID:(-1)coeff} or {mnxm:(-1)coeff}
	    substrates = substrates.split(' + ')
	    mnxm_coeff_dict = {}
	    for substrate in substrates:
	        substrate = substrate.split()
	        if substrate[1] in mnxm_bigg_compound_dict.keys():
		    mnxm_coeff_dict[mnxm_bigg_compound_dict[substrate[1]]] = float(substrate[0])*-1
	        elif substrate[1] in mnxm_kegg_compound_dict.keys():
		    mnxm_coeff_dict[mnxm_kegg_compound_dict[substrate[1]]] = float(substrate[0])*-1
	        else:
		    mnxm_coeff_dict[substrate[1]] = float(substrate[0])*-1

#Creating: e.g., {bigg compoundID:coeff}, {kegg compoundID:coeff} or {mnxm:coeff}
	    products = products.split(' + ')
	    for product in products:
	        product = product.split()
	        if product[1] in mnxm_bigg_compound_dict.keys():
		    mnxm_coeff_dict[mnxm_bigg_compound_dict[product[1]]] = float(product[0])
	        elif product[1] in mnxm_kegg_compound_dict.keys():
		    mnxm_coeff_dict[mnxm_kegg_compound_dict[product[1]]] = float(product[0])
	        else:
		    mnxm_coeff_dict[product[1]] = float(product[0])

#Creating: e.g., {'R03232': {'f1p': -1.0, 'C04261': 1.0, 'fru': 1.0, 'C00615': -1.0}}
	    rxnid_mnxm_coeff_dict[mnxr_kegg_dict[mnxr]] = mnxm_coeff_dict

#Discards polymerization reactions with undefinite coeff's 
#e.g., 1 MNXM9 + (n+2) MNXM90033 = 1 MNXM5617 + (n) MNXM90033
	else:
	    print >>fp2, '%s\t%s' %(mnxr, unparsed_equation)
    fp1.write(str(rxnid_mnxm_coeff_dict))
    for key in rxnid_mnxm_coeff_dict.keys():
	print >>fp1, '%s\t%s' %(key, rxnid_mnxm_coeff_dict[key])

    fp1.close()
    fp2.close()
    return rxnid_mnxm_coeff_dict


#Extracts compound information for name, entry, and formula using compound ID from KEGG
def get_compoundInfo(compoundID):
    url = "http://rest.kegg.jp/get/%s"%(compoundID)
    data = urllib2.urlopen(url).read()
    split_text = data.strip().split('\n')
    NAME = ''
    FORMULA = ''
    for line in split_text:
        sptlist = line.split()
        if sptlist[0].strip() == 'NAME':
            NAME = sptlist[1].strip()
            NAME = NAME.replace(';','')
        elif sptlist[0].strip() == 'FORMULA':
            FORMULA = sptlist[1].strip()
    return {'NAME':NAME, 'FORMULA':FORMULA}


def add_nonBBH_rxn(modelPrunedGPR, rxnid_info_dict, rxnid_mnxm_coeff_dict, rxnid_locusTag_dict, bigg_mnxm_compound_dict, kegg_mnxm_compound_dict, mnxm_compoundInfo_dict, outputFile1):

    fp1 = open(outputFile1, "w")
    fp1.write("Reaction ID"+"\t"+"Reaction name"+"\t"+"Lower bound"+"\t"+"Upper bound"+"\t"+"Reaction equation"+"\t"+"GPR"+"\t"+"Pathway"+"\n")
    for rxnid in rxnid_mnxm_coeff_dict.keys():
	print rxnid
#ID
	rxn = Reaction(rxnid)
	fp1.write(str(rxn.id)+"\t")
#Name
#Soem reaction IDs do not have NAME despite the presence of PATHWAY
	rxn.name = rxnid_info_dict[rxnid]['NAME']
	fp1.write(str(rxn.name)+"\t")
  
#Reversibility / Lower and upper bounds
	rxn.lower_bound = -1000
	rxn.uppwer_bound = 1000
	fp1.write(str(rxn.lower_bound)+"\t")
	fp1.write(str(rxn.upper_bound)+"\t")

#Metabolites and their stoichiometric coeff's
	for metab in rxnid_mnxm_coeff_dict[rxnid]:
	    metab_compt = '_'.join([metab,'c'])

#Adding metabolites already in the model
	    if metab_compt in modelPrunedGPR.metabolites:
		rxn.add_metabolites({modelPrunedGPR.metabolites.get_by_id(metab_compt):rxnid_mnxm_coeff_dict[rxnid][metab]})

#Adding metabolites with bigg compoundID, but not in the model
	    elif metab in bigg_mnxm_compound_dict.keys():
		mnxm = bigg_mnxm_compound_dict[metab]
		metab_compt = Metabolite(metab, formula = mnxm_compoundInfo_dict[mnxm][1], name = mnxm_compoundInfo_dict[mnxm][0], compartment='c')
		rxn.add_metabolites({metab_compt:rxnid_mnxm_coeff_dict[rxnid][metab]})

#Adding metabolites with KEGG compoundID and not in the model
	    else:
		keggID = get_compoundInfo(metab)
		metab_compt = Metabolite(metab, formula = keggID['FORMULA'], name = keggID['NAME'], compartment='c')
		rxn.add_metabolites({metab_compt:rxnid_mnxm_coeff_dict[rxnid][metab]})
	fp1.write(rxn.reaction+"\t")

#GPR association
	if len(rxnid_locusTag_dict[rxnid]) == 1:
	    gpr = '( %s )' %(rxnid_locusTag_dict[rxnid][0])
	else:
	    gpr = ' or '.join(rxnid_locusTag_dict[rxnid])
	    gpr = '( %s )' %(gpr)
	rxn.add_gene_reaction_rule(gpr)
	fp1.write(rxn.gene_reaction_rule+"\t")

#Subsystem
	rxn.subsystem = rxnid_info_dict[rxnid]['PATHWAY']
	fp1.write(str(rxn.subsystem)+"\n")

#E.C. number: not available feature in COBRApy
#Objective coeff: default
	rxn.objective_coefficient = 0

#Addition of a reaction to the model
	modelPrunedGPR.add_reaction(rxn)
    fp1.close()

    target_model = copy.deepcopy(modelPrunedGPR)
    return target_model

