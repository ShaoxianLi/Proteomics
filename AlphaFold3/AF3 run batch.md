### STEP 1: Generate JSON files in batch and file List
Run get json files.ipynb on HPC to get all the JSON files
Make sure this script is in the same folder with the CSV file download from TMT viewer. The UniprotID should be in column 2 and GeneSymbol should be in column 3 
### STEP 2: submit the job
- Before submitting this script, run the following command once:
```
ls $HOME/af3_input_files/citrate_hits/*.json > $HOME/af3_input_files/citrate_hits_list.txt
```

- Then confirm the number of JSONs (for updating the #BSUB -J "af3_json[1-1293]" line):
```
wc -l $HOME/af3_input_files/citrate_hits_list.txt
```

- Save this script as something like ```submit_af3_json_jobs.bash```