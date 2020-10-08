
import pandas as pd
import numpy as np


# Bayesian networks only support discrete values, they are not able to deal with continuous variables
# given that our breast cancer dataset has continous variables, we need to discretise them
def discretize_dataframe( df, num_bins, class_var ):
	r=np.array(range(num_bins+1))/(1.0*num_bins)
	
	# quantiles are building using pandas.qcut
	# The "class" column is just copied.
	l=[]
	
	for col in df.columns.values:
		if col!=class_var:
			l.append(  pd.DataFrame( pd.qcut( df[col],r, duplicates='drop',precision=2),columns=[col]))
		else:
			l.append( pd.DataFrame( np.round(df[col],2),columns=[col]))
			
	treated = pd.concat(l, join='outer', axis=1)
	return treated
	
