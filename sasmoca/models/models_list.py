#!/usr/bin/python

import sys

from models import (LUV_DMPC,
					LUV_POPC,
					LUV_POPE, 
					pLUV_DLPC_OmpLA_RecBuf,
					pLUV_POPC_OmpLA_RecBuf)

from models.virtual_barrel import (Stacked_disks_singleSLD,
								   Stacked_disks_singleSLD_dimer)


###########################################################################################
###########################################################################################
def ChooseFunction (function ):
	" Choose Function "
	
	proteo = False

######################

	if function=="LUV_DMPC":
		intensity	= LUV_DMPC

	elif function=="LUV_POPC":
		intensity	= LUV_POPC

	elif function=="LUV_POPE":
		intensity	= LUV_POPE

	elif function=="pLUV_POPC_OmpLA_RecBuf":
		intensity	= pLUV_POPC_OmpLA_RecBuf
		proteo = True

	elif function=="pLUV_DLPC_OmpLA_RecBuf":
		intensity	= pLUV_DLPC_OmpLA_RecBuf
		proteo = True

	elif function=="Stacked_disks_singleSLD":
		intensity	= Stacked_disks_singleSLD		
	elif function=="Stacked_disks_singleSLD_dimer":
		intensity	= Stacked_disks_singleSLD_dimer	
	else:
		sys.exit("--- This function name does not exist")
	
	return ( intensity, proteo )
	

###########################################################################################
###########################################################################################

