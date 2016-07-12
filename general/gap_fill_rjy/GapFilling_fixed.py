'''
Created on 2014. 6. 4.

@author: user
'''

import copy
import os

from gurobipy import *

import Simulator
from cobra.io.sbml import create_cobra_model_from_sbml_file
from cobra.io.sbml import write_cobra_model_to_sbml_file    

class GapFill(Simulator.Simulator):
    '''
    classdocs
    '''

    def __init__(self):
        '''
        Constructor
        '''
        
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
                 
        pairs, coffvalue = multidict( Smatrix )
        pairs = tuplelist(pairs)         
 
        m = Model('Gap Filling')
        m.setParam('OutputFlag', 0)
        m.reset()
         
        # create variables
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
         
        # Add constraints
        for each_metabolite in model_metabolites:            
            if len( pairs.select(each_metabolite, '*') ) == 0:
                continue
            m.addConstr( quicksum( v[reaction]*coffvalue[metabolite,reaction] for metabolite, reaction in pairs.select(each_metabolite, '*') ) == 0 )            
             
        m.update()
        
        m.setObjective( quicksum((b_bool[each_reaction]) for each_reaction in IrreversibleReaction ), GRB.MINIMIZE)        
        m.optimize()
        
        print 'Model Status : ', m.status        
        if m.status == 2:
            print 'Objective value : ', m.ObjVal
            ReactionFlux = {}
            for reaction in IrreversibleReaction:
                if b_bool[reaction].x == 1.0:
                    rev_changed_reactions.append( reaction )
                    
        return m.status, m.ObjVal, rev_changed_reactions      
    
        
    def run_GapFill(self, target_reaction, flux_constraints={}, UniversalReactions=[], inf_flag = False):
        
        print 'Start Gap filling  ... '
        added_reactions = []
        
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
 
        m = Model('Gap Filling')
        m.setParam('OutputFlag', 0)
        m.reset()
         
        # create variables
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
        
        m.addConstr(v[target_reaction]>= 0.01)
                
        m.update()
         
        # Add constraints
        for each_metabolite in model_metabolites:            
            if len( pairs.select(each_metabolite, '*') ) == 0:
                continue
            m.addConstr( quicksum( v[reaction]*coffvalue[metabolite,reaction] for metabolite, reaction in pairs.select(each_metabolite, '*') ) == 0 )            
             
        m.update()
        
        m.setObjective( quicksum((b_bool[each_reaction]) for each_reaction in UniversalReactions ), GRB.MINIMIZE)        
        m.optimize()
        
        print 'Model Status : ', m.status        
        if m.status == 2:
            print 'Objective value : ', m.ObjVal
            ReactionFlux = {}
            for reaction in UniversalReactions:
                if b_bool[reaction].x == 1.0:                    
                    added_reactions.append( reaction )
                    
        return m.status, m.ObjVal, added_reactions      
    
    def change_reversibility(self,target_reaction, IrreversibleReaction ):        
        status, ObjVal, rev_changed_reaction  =  self.run_Reversibility(target_reaction=target_reaction, IrreversibleReaction=IrreversibleReaction)
        print status
        print ObjVal
        print rev_changed_reaction        
        return

    def fill_gap(self, target_reaction, universal_reaction_model):        
        universal_reactions = [ each_reaction.id for each_reaction in universal_reaction_model.reactions ]
        cobra_model.add_reactions( universal_reaction_model.reactions )
        cobra_model.optimize()
        print cobra_model.solution.f    
        status, ObjVal, added_reaction  =  self.run_GapFill(target_reaction=target_reaction, UniversalReactions=universal_reactions)
        print status
        print ObjVal
        print added_reaction
        return        
    
if __name__ == '__main__':    
    
    
    
    # For Gap Filling
    # Ecoli_MFA.xml (original version)
    # Block_Ecoli_MFA.xml # PtsG KO
    # Universal_Ecoli_MFA.xml
    
    #cobra_model = create_cobra_model_from_sbml_file('./Gapped_Ecoli_MFA.xml')
    #universal_model = create_cobra_model_from_sbml_file('./Universal_Ecoli_MFA.xml')
    #cobra_model.add_reactions( universal_model.reactions )
    #obj = GapFill()
    #obj.loadCobraModel(cobra_model) # load merged model
    #obj.fill_gap( 'biomass', universal_model )
    
    ######################################################
    
    # For Change Reversibility
    # For Gap Filling
    # Ecoli_MFA.xml (original version)
    # Block_Ecoli_MFA.xml # PtsG KO
    # Universal_Ecoli_MFA.xml
    
    cobra_model = create_cobra_model_from_sbml_file('./Gapped_Ecoli_MFA.xml')
    
    IrreversibleReaction  = []
    for each_reaction in cobra_model.reactions:
        if each_reaction.reversibility == 0:
            IrreversibleReaction.append( each_reaction.id )
        each_reaction.reversibility = 1
        each_reaction.lower_bound = -1000.0
        each_reaction.upper_bound = 1000.0

    cobra_model.optimize()
    print cobra_model.solution.f    
    
    obj = GapFill()
    obj.loadCobraModel(cobra_model)
    obj.change_reversibility('biomass',IrreversibleReaction)
