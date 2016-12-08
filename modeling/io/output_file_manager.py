
#Copyright 2014-2016 BioInformatics Research Center, KAIST
#Copyright 2014-2016 Novo Nordisk Foundation Center for Biosustainability, DTU

import logging
import pickle
from cobra.io.sbml import write_cobra_model_to_sbml_file, create_cobra_model_from_sbml_file


def generate_outputs_primary_model(model, modelPrunedGPR, target_model, options):
    #Output files
    #Model reloading and overwrtting are necessary for model consistency:
    #e.g., metabolite IDs with correct compartment suffices & accurate model stats
    #This can also mask the effects of model error (e.g., undeclared metabolite ID)
    #Cobrapy IO module seems to have an error for adding new reactions
    write_cobra_model_to_sbml_file(target_model,
            './%s/2_primary_metabolic_model/target_model_%s.xml'
            %(options.outputfolder, options.orgName))
    target_model = create_cobra_model_from_sbml_file(
            './%s/2_primary_metabolic_model/target_model_%s.xml'
            %(options.outputfolder, options.orgName))
    write_cobra_model_to_sbml_file(target_model,
            './%s/2_primary_metabolic_model/target_model_%s.xml'
            %(options.outputfolder, options.orgName))

    #Output on screen
    model = pickle.load(open('%s/model.p' %(options.input1),'rb'))
    logging.info("Stats of 'primary' metabolic model")
    logging.info("Number of genes: %s; %s; %s" %(len(model.genes),
                len(modelPrunedGPR.genes), len(target_model.genes)))
    logging.info("Number of reactions: %s; %s; %s" %(len(model.reactions),
                len(modelPrunedGPR.reactions), len(target_model.reactions)))
    logging.info("Number of metabolites: %s; %s; %s" %(len(model.metabolites),
                len(modelPrunedGPR.metabolites), len(target_model.metabolites)))

    fp1 = open('./%s/2_primary_metabolic_model/target_model_reactions.txt'
                %options.outputfolder, "w")
    fp2 = open('./%s/2_primary_metabolic_model/target_model_metabolites.txt'
                %options.outputfolder, "w")
    fp1.write("Reaction ID"+"\t"+"Reaction name"+"\t"+"Reaction equation"+"\t"+"GPR"+"\t"+"Pathway"+"\n")
    fp2.write("Metabolite ID"+"\t"+"Metabolite name"+"\t"+"Formula"+"\t"+"Compartment"+"\n")

    for j in range(len(target_model.reactions)):
        rxn = target_model.reactions[j]
        print >>fp1, '%s\t%s\t%s\t%s\t%s\t%s' %(rxn.id, rxn.name,
                rxn.reaction, rxn.gene_reaction_rule, rxn.subsystem)

    for i in range(len(target_model.metabolites)):
        metab = target_model.metabolites[i]
        print >>fp2, '%s\t%s\t%s\t%s' %(metab.id, metab.name, metab.formula,
                metab.compartment)

    fp1.close()
    fp2.close()


def generate_outputs_secondary_model(target_model_complete, options):
    #Output files
    #Model reloading and overwrtting are necessary for model consistency:
    #e.g., metabolite IDs with correct compartment suffices & accurate model stats
    #This can also mask the effects of model error (e.g., undeclared metabolite ID)
    #Cobrapy IO module seems to have an error for adding new reactions
    write_cobra_model_to_sbml_file(target_model_complete,
            './%s/4_complete_model/target_model_complete.xml'
            %options.outputfolder)

    #Output on screen
    model = pickle.load(open('%s/model.p' %(options.input1),'rb'))
    logging.info("Stats of 'secondary' metabolic model")
    logging.info("Number of genes: %s" %len(target_model_complete.genes))
    logging.info("Number of reactions: %s" %len(target_model_complete.reactions))
    logging.info("Number of metabolites: %s" %len(target_model_complete.metabolites))

    fp1 = open('./%s/4_complete_model/target_model_complete_reactions.txt'
            %options.outputfolder, "w")
    fp2 = open('./%s/4_complete_model/target_model_complete_metabolites.txt'
            %options.outputfolder, "w")
    fp1.write("Reaction ID"+"\t"+"Reaction name"+"\t"+"Reaction equation"+"\t"+"GPR"+"\t"+"Pathway"+"\n")
    fp2.write("Metabolite ID"+"\t"+"Metabolite name"+"\t"+"Formula"+"\t"+"Compartment"+"\n")

    for j in range(len(target_model_complete.reactions)):
        rxn = target_model_complete.reactions[j]
        print >>fp1, '%s\t%s\t%s\t%s\t%s\t%s' %(rxn.id, rxn.name,
                rxn.reaction, rxn.gene_reaction_rule, rxn.subsystem)

    for i in range(len(target_model_complete.metabolites)):
        metab = target_model_complete.metabolites[i]
        print >>fp2, '%s\t%s\t%s\t%s' %(metab.id, metab.name, metab.formula,
                metab.compartment)

    fp1.close()
    fp2.close()

