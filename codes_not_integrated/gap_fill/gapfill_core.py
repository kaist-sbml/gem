'''
2015 Jae-Yong Ryu
2015 Kyu-Sang Hwang
2015 Hyun Uk Kim
'''

#gurobipy should be imported in this manner
from gurobipy import *
from cobra import Model, Reaction, Metabolite
from cobra.io.sbml import create_cobra_model_from_sbml_file

class gapfilling_precursor():
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''

    def load_cobra_model(self, cobra_model):
        import time

        self.cobra_model = cobra_model
        model = cobra_model
        model.optimize()
        print('\nSimulated desired product is %1.3f' % model.solution.f)
        model_metabolites = []
        model_reactions = []
        Lower_Boundary_Constraint = {}
        Upper_Boundary_Constraint = {}
        objective_reaction = ''
        for each_metabolite in model.metabolites:
            model_metabolites.append(str(each_metabolite.id))

        Smatrix = {}

        for each_reaction in model.reactions:
            if each_reaction.objective_coefficient == 1.0:
                objective_reaction = str(each_reaction.id)

            reactant_list = each_reaction.get_reactants()
            reactant_coff_list = each_reaction.get_coefficients(reactant_list)
            product_list = each_reaction.get_products()
            product_coff_list = each_reaction.get_coefficients(product_list)

            for i in range(len(reactant_list)):
                Smatrix[(str(reactant_list[i].id), str(each_reaction.id))] = reactant_coff_list[i]

            for i in range(len(product_list)):
                Smatrix[(str(product_list[i].id), str(each_reaction.id))] = product_coff_list[i]

            model_reactions.append(str(each_reaction.id))
            lb = each_reaction.lower_bound
            ub = each_reaction.upper_bound
            if lb < -1000.0:
                lb = float('-inf')
            if ub > 1000.0:
                ub = float('inf')
            Lower_Boundary_Constraint[str(each_reaction.id)] = lb
            Upper_Boundary_Constraint[str(each_reaction.id)] = ub

        self.model_metabolites = model_metabolites

        self.model_reactions = model_reactions

        self.Smatrix = Smatrix

        self.Lower_Boundary_Constraint = Lower_Boundary_Constraint

        self.Upper_Boundary_Constraint = Upper_Boundary_Constraint

        self.objective = objective_reaction

        self.update_model_info(model)

        return (model_metabolites, model_reactions, Smatrix, Lower_Boundary_Constraint, Upper_Boundary_Constraint, objective_reaction)


    def update_model_info(self, model):
        self.model_reaction_info={}
        self.model_genes=[]
        self.model_gene_info={}
        self.model_metabolite_target_Reaction_info={}

        self.subsystem = []
        self.subsystem_reaction_info = {}

        self.model_genes =  model.genes
        self.metabolite_list = []

        for reaction in model.reactions:
            GPR_info = reaction.gene_reaction_rule
            self.model_reaction_info[str(reaction.id)]={}
            self.model_reaction_info[str(reaction.id)]['GPR_str'] = reaction.gene_reaction_rule
            self.model_reaction_info[str(reaction.id)]['genes'] = reaction.get_gene()
            self.model_reaction_info[str(reaction.id)]['products'] = reaction.get_products()
            self.model_reaction_info[str(reaction.id)]['reactants'] = reaction.get_reactants() 

            for each_reactant in self.model_reaction_info[str(reaction.id)]['reactants']:
                if each_reactant not in self.metabolite_list:
                    self.metabolite_list.append(str(each_reactant))

            for each_product in self.model_reaction_info[str(reaction.id)]['products']:
                if each_product not in self.metabolite_list:
                    self.metabolite_list.append(str(each_product))

            for each_reactant in self.model_reaction_info[str(reaction.id)]['reactants']:
                if each_reactant not in self.model_metabolite_target_Reaction_info.keys():
                    self.model_metabolite_target_Reaction_info[str(each_reactant)] = [str(reaction.id)]
                else:
                    self.model_metabolite_target_Reaction_info[str(each_reactant)].append(str(reaction.id))

            self.model_reaction_info[str(reaction.id)]['subsystem'] = str(reaction.subsystem)

            if str(reaction.subsystem) not in self.subsystem_reaction_info.keys():
                self.subsystem_reaction_info[str(reaction.subsystem)] = [str(reaction.id)]
            else:
                self.subsystem_reaction_info[str(reaction.subsystem)].append(str(reaction.id))

            self.subsystem.append(str(reaction.subsystem))

            for gene in reaction.get_gene():
                str_gene = str(gene)
                if str_gene not in self.model_gene_info.keys():
                    self.model_gene_info[str_gene]=[str(reaction.id)]
                else:
                    self.model_gene_info[str_gene].append(str(reaction.id))

        self.metabolite_list = list(set(self.metabolite_list))

        self.model = model.to_array_based_model()
        self.subsystem = list(set(self.subsystem))
        return model


    def change_reversibility(self,target_reaction, cobra_model):
        IrreversibleReaction  = []
        for each_reaction in cobra_model.reactions:
            if each_reaction.reversibility == 0:
                IrreversibleReaction.append(each_reaction.id)
            each_reaction.reversibility = 1
            each_reaction.lower_bound = -1000.0
            each_reaction.upper_bound = 1000.0

        status, ObjVal, rev_changed_reaction  =  self.run_Reversibility(target_reaction=target_reaction, IrreversibleReaction=IrreversibleReaction)
        print status
        print ObjVal
        print rev_changed_reaction
        return rev_changed_reaction


    def run_Reversibility(self, target_reaction, flux_constraints={}, IrreversibleReaction=[], inf_flag = False):

        print 'Change reaction reversibility ... '
        rev_changed_reactions = []
        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions

        Smatrix = self.Smatrix

        LowerFlux=self.Lower_Boundary_Constraint
        UpperFlux=self.Upper_Boundary_Constraint

        if inf_flag == False:
            for key in LowerFlux.keys():
                if LowerFlux[key] == float("-inf"):
                    LowerFlux[key] = -1000.0

            for key in UpperFlux.keys():
                if UpperFlux[key] == float("inf"):
                    UpperFlux[key] = 1000.0

        pairs, coffvalue = multidict(Smatrix)
        pairs = tuplelist(pairs)

        m = gurobipy.Model('Gap Filling')
        m.setParam('OutputFlag', 0)
        m.reset()

        #Create variables
        v = {}
        b_bool = {}

        for each_reaction in model_reactions:
            each_reaction = each_reaction.strip()
            if each_reaction in flux_constraints.keys():
                v[each_reaction] = m.addVar(lb = flux_constraints[each_reaction][0], ub = flux_constraints[each_reaction][1] , name=each_reaction)
            else:
                v[each_reaction] = m.addVar(lb = LowerFlux[each_reaction], ub = UpperFlux[each_reaction] , name=each_reaction)

            b_bool[each_reaction] = m.addVar(vtype = GRB.BINARY, name=each_reaction)

        m.update()
        print IrreversibleReaction
        for each_reaction in IrreversibleReaction:
            m.addConstr(v[each_reaction] >= -1000.0*b_bool[each_reaction])
            m.addConstr(v[each_reaction] <= 1000.0)

        m.addConstr(v[target_reaction]>= 0.01)

        m.update()

        #Add constraints
        for each_metabolite in model_metabolites:
            if len(pairs.select(each_metabolite, '*')) == 0:
                continue
            m.addConstr(quicksum(v[reaction]*coffvalue[metabolite,reaction] for metabolite, reaction in pairs.select(each_metabolite, '*')) == 0)

        m.update()

        m.setObjective(quicksum((b_bool[each_reaction]) for each_reaction in IrreversibleReaction), GRB.MINIMIZE)
        m.optimize()

        print 'Model Status : ', m.status
        if m.status == 2:
            print 'Objective value : ', m.ObjVal
            ReactionFlux = {}
            for reaction in IrreversibleReaction:
                if b_bool[reaction].x == 1.0:
                    rev_changed_reactions.append(reaction)

        return m.status, m.ObjVal, rev_changed_reactions


    def fill_gap(self, target_reaction, cobra_model, universal_reaction_model):
        universal_reactions = [each_reaction.id for each_reaction in universal_reaction_model.reactions]
        cobra_model.add_reactions(universal_reaction_model.reactions)
        cobra_model.optimize()
        print cobra_model.solution.f
        status, ObjVal, added_reaction  =  self.run_GapFill(target_reaction=target_reaction, UniversalReactions=universal_reactions)
        print status
        print ObjVal
        print added_reaction
        return

    #Based on MILP optimization, this function is followed by the equations as below.
    #Satish_Kumar_et_al,_BMC_Bioinformatics,_2007 -> GapFill algorithm
    def run_GapFill(self, target_reaction, flux_constraints={}, UniversalReactions=[]):

        print 'Start Gap filling  ... '
        added_reactions = []

        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions

        Smatrix = self.Smatrix

        LowerFlux=self.Lower_Boundary_Constraint
        UpperFlux=self.Upper_Boundary_Constraint

        for key in LowerFlux.keys():
            if LowerFlux[key] == float("-inf"):
               LowerFlux[key] = -1000.0

        for key in UpperFlux.keys():
            if UpperFlux[key] == float("inf"):
               UpperFlux[key] = 1000.0

        pairs, coeffvalue = multidict(Smatrix)
        pairs = tuplelist(pairs)

        m = gurobipy.Model('Gap Filling')
        m.setParam('OutputFlag', 0)
        m.reset()

        #Create variables
        v = {}
        b_bool = {}

        for each_reaction in model_reactions:
            each_reaction = each_reaction.strip()
            if each_reaction in flux_constraints.keys():
                v[each_reaction] = m.addVar(lb = flux_constraints[each_reaction][0], ub = flux_constraints[each_reaction][1] , name=each_reaction)
            else:
                v[each_reaction] = m.addVar(lb = LowerFlux[each_reaction], ub = UpperFlux[each_reaction] , name=each_reaction)

            b_bool[each_reaction] = m.addVar(vtype = GRB.BINARY, name=each_reaction)

        m.update()

        for each_reaction in UniversalReactions:
            m.addConstr(v[each_reaction]-1000.0*b_bool[each_reaction] <= 0)
            m.addConstr(v[each_reaction]+1000.0*b_bool[each_reaction] >= 0)

        #Check tolerance: too small number may cause infeasibility
        #addConstr modification : 0.01 -> 1.00
        m.addConstr(v[target_reaction]>= 1.00)

        m.update()

        #Add constraints
        for each_metabolite in model_metabolites:
            if len(pairs.select(each_metabolite, '*')) == 0:
                continue
            m.addConstr(quicksum(v[reaction]*coeffvalue[metabolite,reaction] for metabolite, reaction in pairs.select(each_metabolite, '*')) == 0)

        m.update()

        m.setObjective(quicksum((b_bool[each_reaction]) for each_reaction in UniversalReactions), GRB.MINIMIZE)
        m.optimize()

        print 'Model Status : ', m.status
        if m.status == 2:
            print 'Objective value : ', m.ObjVal
            ReactionFlux = {}
            for reaction in UniversalReactions:
                if b_bool[reaction].x == 1.0:
                    added_reactions.append(reaction)

        return m.status, m.ObjVal, added_reactions

