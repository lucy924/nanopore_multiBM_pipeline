# Getting CIBERSORT docker to work

Following /home/dejlu879/ProjectProtocol/cibersort/README_CIBERSORTxFractions.txt

In this directory:  
```
docker pull cibersortx/fractions
```

Now go to the file directories  
```
cd /external/analyses/lucy/nanopore_multiBM_pipeline/results/test3/test_docker
dir_path="/external/analyses/lucy/nanopore_multiBM_pipeline/results/test3/test_docker"
docker run -v ${dir_path}:/src/data -v ${dir_path}:/src/outdir cibersortx/fractions --username lucy.picard@postgrad.otago.ac.nz --token 62c6ce7d617349ca44a7ee8efb7df65c --mixture test3.CS_mix_matrix.txt --sigmatrix test3.CS_bladder_ref.txt --label test3 --perm 1 --QN FALSE --verbose TRUE
```

## Docker result (commandline)
(base) dejlu879@rtis-path-nanopore:/external/analyses/lucy/nanopore_multiBM_pipeline/results/test3/test_docker$ docker run -v ${dir_path}:/src/data -v ${dir_path}:/src/outdir cibersortx/fractions --username lucy.picard@postgrad.otago.ac.nz --token 62c6ce7d617349ca44a7ee8efb7df65c --mixture test3.CS_mix_matrix.txt --sigmatrix test3.CS_bladder_ref.txt --perm 1 --QN FALSE --verbose TRUE
>Running CIBERSORTxFractions...
>[Options] username: lucy.picard@postgrad.otago.ac.nz
>[Options] token: 62c6ce7d617349ca44a7ee8efb7df65c
>[Options] mixture: test3.CS_mix_matrix.txt
>[Options] sigmatrix: test3.CS_bladder_ref.txt
>[Options] perm: 1
>[Options] QN: FALSE
>[Options] verbose: TRUE
>=============CIBERSORTx Settings===============
>Mixture file: /src/data/test3.CS_mix_matrix.txt
>Signature matrix file: /src/data/test3.CS_bladder_ref.txt
>Number of permutations set to: 1
>Enable verbose output
>==================CIBERSORTx===================
P       100     %
column Cancer CD14 CD19 CD4_Eff CD56 CD8 Endothelial Eos Fibroblast Neu Treg P-value Correlation RMSE %Completed
1 0.714953 0 0.0176704 0 0.0248744 0.0134504 0.026305 100

## Website result
>Running CIBERSORTxFractions...
>[Options] mixture: files/lucy.picard@postgrad.otago.ac.nz/test3-CS_mix_matrix.txt
>[Options] sigmatrix: files/lucy.picard@postgrad.otago.ac.nz/test3-CS_bladder_ref.txt
>[Options] perm: 1
>[Options] verbose: TRUE
>[Options] QN: FALSE
>[Options] outdir: files/lucy.picard@postgrad.otago.ac.nz/results/
>[Options] label: Job21
>=============CIBERSORTx Settings===============
>Mixture file: files/lucy.picard@postgrad.otago.ac.nz/test3-CS_mix_matrix.txt 
>Signature matrix file: files/lucy.picard@postgrad.otago.ac.nz/test3-CS_bladder_ref.txt 
>Number of permutations set to: 1 
>Enable verbose output
>==================CIBERSORTx===================
P	100	%	
column Cancer CD14 CD19 CD4_Eff CD56 CD8 Endothelial Eos Fibroblast Neu Treg P-value Correlation RMSE %Completed 
1 0.71552773682421 0 0.0177456741447506 0 0.0256672410471266 0.0121458982222691 0.0263005364546374 100 
        Cancer CD14       CD19 CD4_Eff       CD56       CD8 Endothelial
Beta 0.7155277    0 0.01774567       0 0.02566724 0.0121459  0.02630054
            Eos Fibroblast         Neu     Treg P-value Correlation      RMSE
Beta 0.03113752  0.0234229 0.003669513 0.144383       0   0.8813419 0.6724026
Suffix:  Results
All done.
