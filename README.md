# protein-data-processing

## Argument parameters to tune
### For PDB data
```
filter_date = datetime.datetime(2020, 5, 1)  
resolution_threshold = 9.0  
require_xray = True  # Set to False if X-ray is not required
```

### For AlphaFold Database
```
# Threshold for pLDDT  
plddt_threshold = 0.7  
# Distance threshold for long-range contacts (in Å)  
distance_threshold = 8.0  
# Minimum sequence separation for long-range contacts  
min_seq_distance = 12
```
