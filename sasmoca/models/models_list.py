#!/usr/bin/python

import sys

from models import (pLUV_POPC_OmpLA_RecBuf, 
					pLUV_DLPC_OmpLA_RecBuf)

from models.virtual_barrel import (Stacked_disks_singleSLD,
								   Stacked_disks_singleSLD_dimer)


###########################################################################################
###########################################################################################
def ChooseFunction (function ):
	" Choose Function "
	
	pltOptions = {  }

######################

	if function=="pLUV_POPC_OmpLA_RecBuf":
		intensity	= pLUV_POPC_OmpLA_RecBuf
	elif function=="pLUV_DLPC_OmpLA_RecBuf":
		intensity	= pLUV_DLPC_OmpLA_RecBuf

	elif function=="Stacked_disks_singleSLD":
		intensity	= Stacked_disks_singleSLD		
	elif function=="Stacked_disks_singleSLD_dimer":
		intensity	= Stacked_disks_singleSLD_dimer	
	else:
		sys.exit("--- This function name does not exist")
	
	return ( intensity, pltOptions )
	

###########################################################################################
###########################################################################################

