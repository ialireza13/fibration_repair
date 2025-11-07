param_data = {
    ## this changes the balancing from hard to soft constraints
    ## leave this as true in general
    'hard_flag': True,
    ## David refers to the original constraints which are not necessarily
    ## minimal balanced. Bryant are some constraints which cut out feasible 
    ## solutions but 
    'model_type': 'v1',
    ## both, add_only, or rm_only as strings
    'rm_add_flag': "both",
    'InDegOneFlag': True,
    'save_output': True,
    ##key parameter to change
    ##this is the optimality guarantee which 
    'mip_gap': .05
}