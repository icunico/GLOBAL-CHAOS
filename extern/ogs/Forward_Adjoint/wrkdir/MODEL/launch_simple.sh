#! /bin/bash

#for method in 'L-BFGS-B' 'TNC' 'SLSQP' 'Powell' 'trust-constr'; do
for method in 'L-BFGS-B'; do
	echo $method
	bio=inversion_result_bio_${method}.txt
	rrs=inversion_result_rrs_${method}.txt
	kd=inversion_result_kd_${method}.txt
        python analytic_inversion.py -fb ${bio} -fr ${rrs} -fk ${kd} -m ${method}
done
