#!/usr/bin/python

import sys

from models import (IPLUV_Stk_SDP_POPCbase_RecBuf, 
					IPLUV_Stk_SDP_DLPCbase_RecBuf)

from models.Dummy_barrel import (Stacked_9disks, 
                                 Stacked_disks_singleSLD,
                                 Stacked_disks_singleSLD_dimer)


###########################################################################################
###########################################################################################
def ChooseFunction (function ):
	" Choose Function "
	
	pltOptions = {  }

######################

	if function=="IPLUV_Stk_SDP_POPCbase_RecBuf":
		intensity	= IPLUV_Stk_SDP_POPCbase_RecBuf
	elif function=="IPLUV_Stk_SDP_DLPCbase_RecBuf":
		intensity	= IPLUV_Stk_SDP_DLPCbase_RecBuf

	elif function=="Stacked_9disks":
		intensity	= Stacked_9disks
	elif function=="Stacked_disks_singleSLD":
		intensity	= Stacked_disks_singleSLD		
	elif function=="Stacked_disks_singleSLD_dimer":
		intensity	= Stacked_disks_singleSLD_dimer	
	else:
		sys.exit("--- This function name does not exist")
	
	return ( intensity, pltOptions )
	

###########################################################################################
###########################################################################################

