from specification_datasets import SpecDatasets
import lib.utils as utils
import numpy as np

class CalcAverageFRAPValencyLength(SpecDatasets):
		
	def __init__( self ):
		#SpecDatasetsFRAP.__init__(self, target_dataset )
		self.center = None
		self.suffix_unified = 'FRAP_'
		self.suffixes = ['FRAP{}'.format(i) for i in range(8)]
		
	def convert_save( self ):
		for valency in self.valencies:
			for length in self.lengths:
				prefix = self.filename_edited_matrix(valency, length)
				print( prefix )
				molecular_concentrations = []
				for suffix_ in self.suffixes:
					d = utils.load(self.dir_edited_data, prefix, suffix_)
					molecular_concentrations.append( d['molecular_concentration_in_target_area'] )
				
				d_ = {}
				d_['time_steps'] = d['time_steps']
				d_['molecular_concentration_in_target_area'] = np.mean( molecular_concentrations ,  axis=0 )
				d_['std_molecular_concentration_in_target_area'] = np.std(  molecular_concentrations ,  axis=0 )
				
				utils.save(self.dir_edited_data, prefix, self.suffix_unified, d_)

if __name__ == '__main__':
	
	obj = CalcAverageFRAPValencyLength()
	obj.C_valency_length_FRAP_Control_fine_sampling(frap = True)
	obj.convert_save()
	
	

