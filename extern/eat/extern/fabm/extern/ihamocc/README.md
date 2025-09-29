# FABM-iHAMOCC

This is a [FABM](https://fabm.net) port of the iHAMOCC biogeochemical model which forms part of the [Norwegian Earth System Model version 2 (NorESM2)](https://doi.org/10.5194/gmd-13-2393-2020). It is based on the iHAMOCC code that comes with the master version of [BLOM](https://noresm-docs.readthedocs.io/en/noresm2/model-description/ocn_model.html) as of June 6th 2023.

This FABM port was done as Part of the [Ocean ICU](https://ocean-icu.eu/) project
Model code has been kept as faithful to the original code as possible, but has been changed in a few places (see below).

In addition, FABM-iHAMOCC currently only includes the sediment_bypass option, meaning that there are no explicit representation of sediment dynamics, beyond a redistribution of tracers into the water column.

The code has been split up into submodules of which several are optional and can be in- or excluded by adjusting the runtime configuration (`fabm.yaml`), no code change or recompilation needed. 
The following core modules are interdependent, and should be considered mandatory, for the time being:
* alkalinization.F90    | Surface alkalinity flux
* carbon.F90            | Abiotic carbon chemistry
* detritus.F90          | Remineralization and export of biological production
* iron.F90              | Surface iron deposition
* light.F90             | Light absorption and flux-at-depth
* nitrogen.F90          | Nitrogen and laughing gas surface flux, fixation and nitrification
* oxygen.F90            | Oxygen surface exchange and solubility
* phytoplankton.F90     | Phytoplankton primary production
* tracer.F90  | Infrastructure for setting up tracers that do not have their own module, such as doc and phosphate
* sediment_bypass.F90   | Redistributes bottom fluxes to the water column
* zooplankton.F90       | Zooplankton production and grazing

The following modules are can be considered optional:
* bromo.F90             | Bromoform surface exchange, production and breakdown
* cfc.F90               | cfc gas surface exchange
* cisonew.F90           | carbon 13 and 14 isotope ratios
* dms.F90               | dimethyl sulfide surface exchange, production and breakdown
* natdic.F90            | natural carbon tracers

Lastly, the following modules provide direct or diagnostic support for other modules:
* mixed_layer           | Calculates mixed layer depth
* preformed_tracer.F90  | provides diagnostics of oxygen, alkalinity, co2 and phosphate before mixing
* shared.F90            | provides fixed parameter values and functions that are common across modules

Currently, the code is not modularized to the extent that e.g. several phytoplankton and zooplankton instances can be run simultaneously. This is mainly due to the complex way that processes in detritus, phytoplankton and zooplankton are interlinked. While a more comprehensive modularization is probably possible, it would require creating many aggregate variables to hold combined contributions of e.g. POC production from each phytoplankton etc.

## How to build

This code must be compiled together with FABM. To do this, provide the following additional arguments to cmake [when you build FABM](https://github.com/fabm-model/fabm/wiki/Building-and-installing): `-DFABM_INSTITUTES=ihamocc -DFABM_IHAMOCC_BASE=<IHAMOCCDIR>`

Here, `<IHAMOCCDIR>` is the directory with the FABM-IHAMM code (the same directory that contains this readme file). Note that `-DFABM_INSTITUTES=ihamocc` will make FABM compile IHAMOCC as the *only* available biogeochemical model. If you additionally want to have access to other biogeochemical models included with FABM, you can set `FABM_INSTITUTES` to a semi-colon separated list, e.g., `-DFABM_INSTITUTES="ihamocc;ersem"` (to prevent the shell from interpreting the semi-colons, you typically have to enclose this list with quotes).

For instance, to use IHAMOCC with the latest stable release of the [General Ocean Turbulence Model (GOTM)](https://gotm.net/), do the following:

```
git clone --recurse-submodules -b v6.0 https://github.com/gotm-model/code.git gotm
git clone https://github.com/fabm-model/fabm.git
git clone https://github.com/BoldingBruggeman/fabm-ihamocc.git
mkdir build
cd build
cmake ../gotm -DFABM_BASE=../fabm -DFABM_INSTITUTES=ihamocc -DFABM_IHAMOCC_BASE=../fabm-ihamocc
make install
```

This will install the GOTM executable with support for IHAMOCC at `~/local/gotm/bin/gotm`.

## How to run a FABM-IHAMOCC simulation

A `fabm.yaml` file with the iHAMOCC configuration is provided under `<IHAMOCCDIR>/testcases`. You can drop this file in the working directory of a FABM-compatible model such as GOTM to use it during simulation. Note that for GOTM, you will also need to ensure that `fabm/use` is set to `true` in `gotm.yaml`. Otherwise GOTM would run with physics only.

## To do

* Create a sediment module for FABM-iHAMOCC, based on the BLOM-iHAMOCC sediment module.
* There is a known bug in bromo.F90, line 119, where there is a unit mismatch in the abiotic losses equation for bromoform.
* If possible, recode the model to be fully modularized, to take full advantage of FABM

## Differences from the BLOM version of iHAMOCC

* Zooplankton feeding formulation has been changed to a standard Michaelis-Menten dynamics equation. The original version was in a non-standard form, and included references to primary production rate and minimum phytoplankton concentration. As this expression was unclear in its derivation, required additional unfortunate dependencies and was difficult to convert from discrete to continuous-time form, it was changed to the version of zooplankton grazing from [Six and Meier-Reimer 1996](https://doi.org/10.1029/96GB02561).
* detritus.F90 and nitrogen.F90 both requires mixed layer depth as a dependency. As iHAMOCC does not include a mixed layer depth model, and as not all physical models provide mixed layer depth, we have adopted the PISCES mixed layer depth model for use with FABM-IHAMOCC, and included it in its distribution.
* currently, the sediment module from iHAMOCC has not been converted to FABM. All setup therefore must use the sediment_bypass option for now.
* ...