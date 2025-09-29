#! /bin/bash

pickle_file=fussmann.pickle
echo picklefile $pickle_file
sbatch --export=ALL,pickle_file=$pickle_file run_sensitivity_template.sbatch
