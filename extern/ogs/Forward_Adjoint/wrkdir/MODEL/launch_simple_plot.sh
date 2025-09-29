#! /bin/bash

#for method in 'L-BFGS-B' 'TNC' 'SLSQP' 'Powell' 'trust-constr'; do
for method in 'L-BFGS-B'; do
	echo $method
        python plot_bio_results.py  -m ${method}
        python plot_rrs_results.py  -m ${method}
done
