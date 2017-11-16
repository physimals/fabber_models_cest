#!/bin/sh

/Volumes/HDD/Research/Oxford/Fabber/fabber_models_cest/build/fabber_cest \
 --data=Data_SNR_B103.nii.gz --mask=Data_SNR_B103.nii.gz  \
 --method=vb --noise=white --model=cest --data-order=singlefile \
 --max-iterations=20 --t12prior --output=SNR_B103_qMT_t12prior  \
 --spec=dataspec_B103.txt --pools=poolmat_ThreePool.txt --ptrain=ptrain.txt  \
 --save-model-fit --satspoil --lineshape=superlorentzian --TR=4 --EXFA=10