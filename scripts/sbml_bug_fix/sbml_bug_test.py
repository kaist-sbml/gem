import cobra.test
from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import write_cobra_model_to_sbml_file, create_cobra_model_from_sbml_file


def add_new_rxn():
    rxn_dict = {'R00618':{'thmtp': 1.0, 'h': -1.0, 'h2o': 1.0, 'pi': -1.0, 'thmpp': -1.0}, 'R00125':{'h': -2.0, 'ap4a': 1.0, 'h2o': 1.0, 'adp': -2.0}, 'R01968':{'dgmp': -1.0, 'h2o': -1.0, 'pi': 1.0, 'dgsn': 1.0}}

    model = create_cobra_model_from_sbml_file('iMK1208.xml')
    #model = cobra.test.create_test_model("ecoli")
    #model = cobra.test.create_test_model("salmonella")

    for rxnid in rxn_dict.keys():

        rxn = Reaction(rxnid)
        rxn.name = rxnid
        rxn.lower_bound = -1000
        rxn.uppwer_bound = 1000

        #Metabolites and their stoichiometric coeff's
        for metab in rxn_dict[rxnid].keys():
            metab_compt = '_'.join([metab,'c'])

            #Adding metabolites already in the model
            if metab_compt in model.metabolites:
                print "Metabolite %s: Already present in model" %metab_compt
                rxn.add_metabolites({model.metabolites.get_by_id(
                    metab_compt):rxn_dict[rxnid][metab]})

            #Adding metabolites not in the model
            else:
                print "Metabolite %s: To be added" %metab
                metab_compt = Metabolite(metab_compt, compartment='c')
                rxn.add_metabolites({metab_compt:rxn_dict[rxnid][metab]})

            rxn.gene_reaction_rule = '( dummy_gene )'
            rxn.objective_coefficient = 0

        model.add_reaction(rxn)

    write_cobra_model_to_sbml_file(model, './test_model.xml', use_fbc_package=False)

if __name__ == '__main__':
    add_new_rxn()
