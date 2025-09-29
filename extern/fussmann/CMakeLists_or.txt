add_library(fabm_models_fussmann OBJECT
            fussmann_model_library.F90
            P1.F90
	    C1.F90
	    X.F90
	    Y.F90
	    Z.F90
	    )


add_dependencies(fabm_models_fussmann fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2014 - Karsten Bolding and Jorn Bruggeman (BB)         !
#-----------------------------------------------------------------------
