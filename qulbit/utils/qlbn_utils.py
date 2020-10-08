
# Given an experiment string name, generate the corresponding title for plots
def get_title( experiment ):
	"""
	experiment: a string containing an experiment name
	"""
	
	experiment_title = experiment
	
	if experiment == "pd_shafir_tversky":
		experiment_title = "Shafir & Tversky (1992) Prisoners Dilemma Experiment"

	if experiment == "pd_li_taplin_game1":
		experiment_title = "Li & Taplin (2002) Prisoners Dilemma Experiment #1"

	if experiment == "pd_li_taplin_game2":
		experiment_title = "Li & Taplin (2002) Prisoners Dilemma Experiment #2"

	if experiment == "pd_li_taplin_game3":
		experiment_title = "Li & Taplin (2002) Prisoners Dilemma Experiment #3"

	if experiment == "pd_li_taplin_game4":
		experiment_title = "Li & Taplin (2002) Prisoners Dilemma Experiment #4"

	if experiment == "pd_li_taplin_game5":
		experiment_title = "Li & Taplin (2002) Prisoners Dilemma Experiment #5"

	if experiment == "pd_li_taplin_game6":
		experiment_title = "Li & Taplin (2002) Prisoners Dilemma Experiment #6"

	if experiment == "pd_li_taplin_game7":
		experiment_title = "Li & Taplin (2002) Prisoners Dilemma Experiment #7"
		
	return experiment_title 
	
