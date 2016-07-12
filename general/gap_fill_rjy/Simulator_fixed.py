'''
Created on 2014. 6. 4.

@author: user
'''

import os
import copy
import cobra.io as io
from gurobipy import *
#from MBEL.Model.GPR_Manipulation import convertGPRstringToListFormat
import numpy as np
from pandas import DataFrame

class Simulator(object):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
    def get_simulation_condition(self):
        
        return
    
    def run_flux_tuning(self, target_info_dic, wild_opt_flux, flux_constraints={}, manipulation_level = 0.2, method='moma', inf_flag = False):
        import itertools  
        print 'Start Manipulate_flux ... '
        targets = target_info_dic.keys()
        
        each_target_info = {}
        manipulation_combination = []
        all_target_reactions = []
        index_cnt = 0
        for each_target in targets:
            target_reactions = target_info_dic[ each_target ]['reaction']
            all_target_reactions+=target_reactions
            mode = target_info_dic[ each_target ]['mode']
            each_target_info[ each_target ] = [ target_reactions ]
            if mode == 'a':
                manipulation_ratio_list = list(np.arange( manipulation_level, 1, manipulation_level))  
                manipulation_ratio_list = np.around( manipulation_ratio_list, 3)              
            elif mode == 'i':
                manipulation_ratio_list = list(np.arange( 1, 2+manipulation_level, manipulation_level))
                manipulation_ratio_list = np.around( manipulation_ratio_list, 3)    
            print manipulation_ratio_list
            manipulation_combination.append( manipulation_ratio_list )
            each_target_info[ each_target ] = [ target_reactions, index_cnt ]
            index_cnt+=1
            
        all_target_reactions = list(set(all_target_reactions))           
        all_manipulation_candidates = list(itertools.product(*manipulation_combination))        

        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions
        
        objective = self.objective
            
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
                 
        pairs, coffvalue = multidict( Smatrix )
        pairs = tuplelist(pairs)         
 
        m = Model('Manipulate_flux_sum')
        m.setParam('OutputFlag', 0)
        m.reset()

        # create variables
        v = {}         
        for each_reaction in model_reactions:   
            if each_reaction in flux_constraints.keys():
                v[each_reaction] = m.addVar(lb = flux_constraints[each_reaction][0], ub = flux_constraints[each_reaction][1] , name=each_reaction)             
            else:
                v[each_reaction] = m.addVar(lb = LowerFlux[each_reaction], ub = UpperFlux[each_reaction] , name=each_reaction)
        m.update()     
    
        # Mass balance        
        for each_metabolite in model_metabolites:            
            if len( pairs.select(each_metabolite, '*') ) == 0:
                continue
            m.addConstr( quicksum( v[reaction]*coffvalue[metabolite,reaction] for metabolite, reaction in pairs.select(each_metabolite, '*') ) == 0 )  
        m.update()
       
        Result_Dic = {}
        
        for manipulation_candidate in all_manipulation_candidates:
            reaction_flux_ratio_info = {} 
            for i in range(len(each_target_info)):
                target = targets[i]
                index = each_target_info[target][1]
                k_info = manipulation_candidate[index]
                target_reactions = each_target_info[target][0] 
                for each_reaction in target_reactions:
                    reaction_flux_ratio_info[ each_reaction ] = k_info
            
            const_reaction_info_dic = {}     
            addr_list = []                
            for each_target_reaction in all_target_reactions:                  
                const_reaction_info_dic[ each_target_reaction ] = reaction_flux_ratio_info[each_target_reaction] * wild_opt_flux[each_target_reaction]
                                                         
            m.update() 
            #m.setObjective(quicksum( (v[each_target_reaction] - const_reaction_info_dic[each_target_reaction])*(v[each_target_reaction] - const_reaction_info_dic[each_target_reaction]) for each_target_reaction in all_target_reactions ), GRB.MINIMIZE)
            #m.optimize()
            #moma_status = m.status            
            #print m.status
            
            moma_const_dic = {}
            moma_status = 2
            if moma_status == 2:
                for each_reaction in all_target_reactions:
                    moma_const_dic[ each_reaction ] = const_reaction_info_dic[each_reaction]                
                                 
                for each_addr in addr_list:
                    removeConstraintIndex = m.getConstrs().index(each_addr)
                    m.remove(m.getConstrs()[removeConstraintIndex])
                m.update()    
    
                for each_target_reaction in all_target_reactions:                    
                    addr = m.addConstr( v[ each_target_reaction ] == (moma_const_dic[ each_target_reaction ]) )
                    addr_list.append( addr )     
                m.update()         
                       
                if method == 'fba':
                    m.setObjective(v[objective], GRB.MAXIMIZE)
                elif method == 'moma':
                    m.setObjective(quicksum( (v[each_reaction] - wild_opt_flux[each_reaction])*(v[each_reaction] - wild_opt_flux[each_reaction]) for each_reaction in model_reactions ), GRB.MINIMIZE)
                else:
                    m.setObjective(v[objective], GRB.MAXIMIZE)                    
                
                m.update()
                m.optimize()    
                if m.status == 2:
                    simulation_obj = m.ObjVal
                    print 'Objective : ',simulation_obj
                    each_flux_dist = {}
                    for each_reaction in model_reactions:
                        each_flux_dist[each_reaction] = v[each_reaction].x  
 
                    # make key
                    Item_Key = ''
                    for each_key in reaction_flux_ratio_info:
                        each_value =  reaction_flux_ratio_info[ each_key ]
                        Item_Key = Item_Key+'[%s_%.2f]#'%(each_key,each_value)
 
                    Result_Dic[ Item_Key ] = each_flux_dist
            
            for each_addr in addr_list:
                removeConstraintIndex = m.getConstrs().index(each_addr)
                m.remove(m.getConstrs()[removeConstraintIndex])
                
            m.update()      
        df = DataFrame.from_dict( Result_Dic )                
        return df.T, Result_Dic
        
    def run_MOMA(self, flux_constraints={}, wild_opt_flux={},  inf_flag = False, norm = 'Euclidean'):
        
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
                 
        pairs, coffvalue = multidict( Smatrix )
        pairs = tuplelist(pairs)         
 
        m = Model('MOMA')
        m.setParam('OutputFlag', 0)
        m.reset()
        
        # create variables
        v = {}         
            
        for each_reaction in model_reactions:   
            if each_reaction in flux_constraints.keys():
                v[each_reaction] = m.addVar(lb = flux_constraints[each_reaction][0], ub = flux_constraints[each_reaction][1] , name=each_reaction)             
            else:
                v[each_reaction] = m.addVar(lb = LowerFlux[each_reaction], ub = UpperFlux[each_reaction] , name=each_reaction)
                            
        m.update()
        
        # Add constraints
        for each_metabolite in model_metabolites:            
            if len( pairs.select(each_metabolite, '*') ) == 0:
                continue
            m.addConstr( quicksum( v[reaction]*coffvalue[metabolite,reaction] for metabolite, reaction in pairs.select(each_metabolite, '*') ) == 0 )            
           
        m.update()
        
        if norm == 'Euclidean':
            m.setObjective(quicksum( (v[each_reaction] - wild_opt_flux[each_reaction])*(v[each_reaction] - wild_opt_flux[each_reaction]) for each_reaction in model_reactions ), GRB.MINIMIZE)        
            m.optimize()
            
                
        if m.status == 2:        
            ReactionFlux = {}
            for reaction in model_reactions:
                ReactionFlux[reaction] = float(v[reaction].x)
                if abs(float(v[reaction].x)) <= 1e-6:
                    ReactionFlux[reaction] = 0.0
                #print reaction, float(v[reaction].x)
                
            return m.status, m.ObjVal, ReactionFlux       
        else:
            return m.status, False, False
            
    def run_FBA(self, new_objective='', flux_constraints={}, inf_flag = False, internal_min_flag = False, mode='max'):
        #print 'Start simple FBA simulation ... '
        
        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions
        
        if new_objective == '':
            objective = self.objective
        else:
            objective = new_objective
            
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
                 
        pairs, coffvalue = multidict( Smatrix )
        pairs = tuplelist(pairs)         
 
        m = Model('FBA')
        m.setParam('OutputFlag', 0)
        m.reset()
         
        # create variables
        v = {}         
        for each_reaction in model_reactions:   
            if each_reaction in flux_constraints.keys():
                v[each_reaction] = m.addVar(lb = flux_constraints[each_reaction][0], ub = flux_constraints[each_reaction][1] , name=each_reaction)             
            else:
                v[each_reaction] = m.addVar(lb = LowerFlux[each_reaction], ub = UpperFlux[each_reaction] , name=each_reaction)
            
 
        m.update() 
         
        # Add constraints
        for each_metabolite in model_metabolites:            
            if len( pairs.select(each_metabolite, '*') ) == 0:
                continue
            m.addConstr( quicksum( v[reaction]*coffvalue[metabolite,reaction] for metabolite, reaction in pairs.select(each_metabolite, '*') ) == 0 )            
             
        m.update()

        if mode == 'max':
            m.setObjective(v[objective], GRB.MAXIMIZE)        
        elif mode == 'min':
            m.setObjective(v[objective], GRB.MINIMIZE)
        m.optimize()
      
        if m.status == 2:
            obj_value = m.ObjVal
            #print 'Objective value : ', m.ObjVal
            if internal_min_flag == True:
                m.addConstr( v[objective] == obj_value )
                m.update()
                m.setObjective(quicksum( (v[each_reaction] * v[each_reaction] ) for each_reaction in model_reactions ), GRB.MINIMIZE)
                m.optimize()
                if m.status == 2:
                    obj_value = m.ObjVal
                    print 'Flux minimization objective value : ', m.ObjVal
                    ReactionFlux = {}
                    for reaction in model_reactions:
                        ReactionFlux[reaction] = float(v[reaction].x)
                        if abs(float(v[reaction].x)) <= 1e-6:
                            ReactionFlux[reaction] = 0.0                    
            else:                   
                ReactionFlux = {}
                for reaction in model_reactions:
                    ReactionFlux[reaction] = float(v[reaction].x)
                    if abs(float(v[reaction].x)) <= 1e-6:
                        ReactionFlux[reaction] = 0.0
                
            return m.status, obj_value, ReactionFlux 
        return m.status, False, False 
        
    def run_simple_FBA(self, new_objective='', inf_flag = False):
        print 'Start simple FBA simulation ... '
        
        model_metabolites = self.model_metabolites
        model_reactions = self.model_reactions
        
        if new_objective == '':
            objective = self.objective
        else:
            objective = new_objective
            
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
                 
        pairs, coffvalue = multidict( Smatrix )
        pairs = tuplelist(pairs)         
 
        m = Model('FBA')
        m.setParam('OutputFlag', 0)
        m.reset()
         
        # create variables
        v = {}         
        for each_reaction in model_reactions:                
            v[each_reaction] = m.addVar(lb = LowerFlux[each_reaction], ub = UpperFlux[each_reaction] , name=each_reaction)
 
        m.update()
         
        # Add constraints
        for each_metabolite in model_metabolites:            
            if len( pairs.select(each_metabolite, '*') ) == 0:
                continue
            m.addConstr( quicksum( v[reaction]*coffvalue[metabolite,reaction] for metabolite, reaction in pairs.select(each_metabolite, '*') ) == 0 )            
             
        m.update()
 
        m.setObjective(v[objective], GRB.MAXIMIZE)        
        m.optimize()
        
        print 'Model Status : ', m.status        
        if m.status == 2:
            #print 'Objective value : ', m.ObjVal
            ReactionFlux = {}
            for reaction in model_reactions:
                ReactionFlux[reaction] = float(v[reaction].x)
                if abs(float(v[reaction].x)) <= 1e-6:
                    ReactionFlux[reaction] = 0.0
                print reaction, ReactionFlux[reaction]
        
        return m.status, m.ObjVal, ReactionFlux 
        
    def readModel(self, filename):
        ## save model information
        model = io.sbml.create_cobra_model_from_sbml_file(filename) 
        return self.loadCobraModel( model )
    
    def loadCobraModel(self, cobra_model):       
        self.cobra_model = cobra_model 
        model = cobra_model        
        model.optimize()
        print('\nSimulated growth rate is %1.3f' % model.solution.f)
        model_metabolites = []
        model_reactions = []
        Lower_Boundary_Constraint = {}
        Upper_Boundary_Constraint = {}
        objective_reaction = ''
        for each_metabolite in model.metabolites:
            model_metabolites.append( str(each_metabolite.id) ) 
 
        Smatrix = {}
        
        for each_reaction in model.reactions:            
            if each_reaction.objective_coefficient == 1.0:
                objective_reaction = str(each_reaction.id)
                
            reactant_list = each_reaction.get_reactants()
            reactant_coff_list = each_reaction.get_coefficients( reactant_list)
            product_list = each_reaction.get_products()
            product_coff_list = each_reaction.get_coefficients( product_list )
            
            for i in range( len(reactant_list)):
                Smatrix[ ( str(reactant_list[i].id), str(each_reaction.id) ) ] = reactant_coff_list[i]      
                
            for i in range( len(product_list)):
                Smatrix[ ( str(product_list[i].id), str(each_reaction.id) ) ] = product_coff_list[i]
                 
            model_reactions.append( str(each_reaction.id) )
            lb = each_reaction.lower_bound
            ub = each_reaction.upper_bound
            if lb < -1000.0:
                lb = float('-inf')
            if ub > 1000.0:
                ub = float('inf')
            Lower_Boundary_Constraint[ str(each_reaction.id) ] = lb
            Upper_Boundary_Constraint[ str(each_reaction.id) ] = ub
        
        self.model_metabolites = model_metabolites
        self.model_reactions = model_reactions
        self.Smatrix = Smatrix
        self.Lower_Boundary_Constraint = Lower_Boundary_Constraint
        self.Upper_Boundary_Constraint = Upper_Boundary_Constraint
        self.objective = objective_reaction
        
        self.update_model_info( model )
        
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
            #print reaction           
            self.model_reaction_info[str(reaction.id)]={}
            self.model_reaction_info[str(reaction.id)]['GPR_str'] = reaction.gene_reaction_rule
            #self.model_reaction_info[str(reaction.id)]['GPR_list'] = convertGPRstringToListFormat( reaction.gene_reaction_rule )
            self.model_reaction_info[str(reaction.id)]['genes'] = reaction.get_gene()
            self.model_reaction_info[str(reaction.id)]['products'] = reaction.get_products()
            self.model_reaction_info[str(reaction.id)]['reactants'] = reaction.get_reactants()       
            
            for each_reactant in self.model_reaction_info[str(reaction.id)]['reactants']:
                if each_reactant not in self.metabolite_list:
                    self.metabolite_list.append( str(each_reactant) )
                    
            for each_product in self.model_reaction_info[str(reaction.id)]['products']:
                if each_product not in self.metabolite_list:
                    self.metabolite_list.append( str(each_product) )
            
            for each_reactant in self.model_reaction_info[str(reaction.id)]['reactants']:
                if each_reactant not in self.model_metabolite_target_Reaction_info.keys():
                    self.model_metabolite_target_Reaction_info[ str(each_reactant) ] = [str(reaction.id)]
                else:
                    self.model_metabolite_target_Reaction_info[ str(each_reactant) ].append( str(reaction.id) )                    
            
            self.model_reaction_info[str(reaction.id)]['subsystem'] = str(reaction.subsystem)
            
            if str(reaction.subsystem) not in self.subsystem_reaction_info.keys():
                self.subsystem_reaction_info[str(reaction.subsystem)] = [str(reaction.id)]
            else:
                self.subsystem_reaction_info[str(reaction.subsystem)].append( str(reaction.id) )
             
            self.subsystem.append( str(reaction.subsystem) )
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
        
    def get_reactions_from_multiple_gene_deletion(self, gene_list):        
        reaction_list = []
        reaction_candidates = []
        for each_gene in gene_list:
            gene_name = str(each_gene)        
            tmp_reaction_candidates = self.model_gene_info[gene_name]
            reaction_candidates+=tmp_reaction_candidates
            
        reaction_candidates = list(set(reaction_candidates))
        
        for reaction in reaction_candidates:
            boolean_list = []
            GPR_list = self.model_reaction_info[reaction]['GPR_list']
            for each_item in GPR_list:
                if type(each_item) == list:
                    tmp_boolean_list = []
                    for each_gene in each_item:
                        if each_gene in gene_list:
                            tmp_boolean_list.append(0.0)
                        else:
                            tmp_boolean_list.append(1.0)
                    value=1        
                    for j in range(len(tmp_boolean_list)):
                        value = value * tmp_boolean_list[j]
                    boolean_list.append(value)                           
                else:
                    if each_item in gene_list:
                        boolean_list.append(0.0)
                    else:
                        boolean_list.append(1.0)
            value=0        
            for i in range(len(boolean_list)):
                value = value + boolean_list[i]
            if value == 0:
                reaction_list.append(reaction)                
        return reaction_list
    
    def get_metabolic_models_from_active_reactions(self, FluxInfo):
        cobra_model = copy.deepcopy( self.cobra_model )
        remove_reaction_candidate_list = []
        
        for each_reaction in cobra_model.reactions:
            if each_reaction.id in FluxInfo:
                if FluxInfo[each_reaction.id] == 0.0:
                    remove_reaction_candidate_list.append( each_reaction )
        cobra_model.remove_reactions( remove_reaction_candidate_list ) 
        cobra_model.optimize()
        print 'Simulated Growth : ', cobra_model.solution.f           
        return cobra_model
    
    def get_cobra_model(self):        
        return copy.deepcopy(self.cobra_model) 
                    
if __name__ == '__main__':    
    import time    
    obj = Simulator()
    obj.readModel( './Ecoli_MFA.xml')    
 
    
