add_library(fabm_models_ihamocc OBJECT
            ihamocc_model_library.F90
            shared.F90
            oxygen.F90
            carbon.F90
            bromo.F90
            cfc.F90
            dms.F90
            alkalinization.F90
            iron.F90
            nitrogen.F90
            preformed_tracer.F90
            phytoplankton.F90
            zooplankton.F90
            detritus.F90
            cisonew.F90
            natdic.F90
            mixed_layer.F90
            sediment_bypass.F90
            light.F90
            tracer.F90
           )

target_link_libraries(fabm_models_ihamocc PRIVATE fabm_base)

#-----------------------------------------------------------------------
# Copyright (C) 2023 - BB
#-----------------------------------------------------------------------
