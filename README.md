### Installation
- ```pymatgen```, ```spinney``` package is needed

- POTCAR directory should be correctly configured for ```pymatgen```

### Sample input file
```json
{
    "system_name": "Gallium_Oxide",
    "structure":{
        "path": "init_structure",
        "file": "POSCAR"
    },
    "submit_queue": "std40",
    "defects":[
        {
            "type": "remove",
            "init_atom": "Ga",
            "charge": [-1, 0]
        }
    ],
    "post_process": {
        "correction": true,
        "gap_range": [0, 3]
    } 
}
```
Please save it in json file format

- The initial structure is in ./init_structure/POSCAR
- submit_queue: ["std40","std2","f96"]
- defects type: ["remove", "replace"]
- gap_range starts from 0, but can be a little bigger to ensure the results
- correction keyword should be true

### Submission
Submit the job by this command
```bash
python main.py -i input.json
# python main.py --input input.json
```