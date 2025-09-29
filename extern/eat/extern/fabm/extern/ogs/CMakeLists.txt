add_library(fabm_models_ogs OBJECT
            ogs_model_library.F90
            bfm_pelagic_base.F90
            Phyto.F90
            PelBac.F90
            PelChem.F90
            PelOxygen.F90
            MicroZoo.F90
            MesoZoo.F90
            shared.F90
            light.F90
            light_spectral.F90
            lidata.F90
            getrmud.F90 
            zenith_angle.F90
            CalciteDissolution.F90
            PelagicCSYS.F90
	    BenthicLayer.F90	    
           )

target_link_libraries(fabm_models_ogs PRIVATE fabm_base)
add_subdirectory(Forward_Adjoint EXCLUDE_FROM_ALL)
install(TARGETS adj EXPORT fabmConfig DESTINATION ${CMAKE_INSTALL_LIBDIR})
target_link_libraries(fabm_models_ogs PUBLIC adj)


