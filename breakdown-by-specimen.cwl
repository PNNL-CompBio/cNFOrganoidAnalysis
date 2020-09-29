label: breakdownfile-tool
id: breakdownfile-tool
cwlVersion: v1.0


class: CommandLineTool
baseCommand:
- python
- breakdownfiles.py


requirements:
 - class: InlineJavascriptRequirement
 - class: DockerRequirement
   dockerPull: amancevice/pandas
 - class: InitialWorkDirRequirement
   listing:
     - entryname: breakdownfiles.py
       entry: |
         #!/usr/bin/env python
         import json
         import sys
         import pandas as pd
         query_tsv=sys.argv[1]
         group_by_column=sys.argv[2]
         res = pd.read_csv(query_tsv,delimiter='\t')
         gdf = res.groupby(group_by_column)
         names = []
         allids = []
         for key,value in gdf:
            names.append(key)
            ids = [i for i in gdf.get_group(key)['id']]
            allids.append(ids)
         res={'names':names, 'ids':[a[0] for a in allids]}
         with open('cwl.json','w') as outfile:
           json.dump(res,outfile)

inputs:

- id: query_tsv
  type: File
  inputBinding:
    position: 1

- id: group_by_column
  type: string
  default: "individualID"
  inputBinding:
    position: 2

outputs:

- id: names
  type: string[]
  outputBinding:
    glob: cwl.json
    loadContents: true
    outputEval: $(JSON.parse(self[0].contents)['names'])

- id: id_array
  type: string[]
  outputBinding:
    glob: cwl.json
    loadContents: true
    outputEval: $(JSON.parse(self[0].contents)['ids'])
