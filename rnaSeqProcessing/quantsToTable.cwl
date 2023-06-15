class: Workflow
label: salmon-alignment-from-synapse
id: salmon-alignment-from-synapse
cwlVersion: v1.0

inputs:
  synapse_config:
    type: File
  idquery:
    type: string
  sample_query:
    type: string
  group_by:
    type: string
  tableparentid:
    type: string[]
  tablename:
    type: string[]
  parentid:
    type: string

requirements:
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement

outputs:
  merged:
    type: File
    outputSource: harmonize-counts/merged
steps:
    get-fv:
       run: https://raw.githubusercontent.com/Sage-Bionetworks/synapse-client-cwl-tools/v0.1/synapse-query-tool.cwl
       in:
         synapse_config: synapse_config
         query: idquery
       out: [query_result]
    get-samples-from-fv:
      run: breakdown-by-specimen.cwl
      in:
         query_tsv: get-fv/query_result
         group_by_column: group_by
      out: [names,id_array]
    get-files:
      run: https://raw.githubusercontent.com/Sage-Bionetworks-Workflows/dockstore-tool-synapseclient/master/cwl/synapse-get-tool.cwl
      scatter: synapseid
      in:
        synapse_config: synapse_config
        synapseid:
            source: get-samples-from-fv/id_array
      out: [filepath]
    get-clinical:
       run: https://raw.githubusercontent.com/Sage-Bionetworks/synapse-client-cwl-tools/v0.1/synapse-query-tool.cwl
       in:
         synapse_config: synapse_config
         query: sample_query
       out: [query_result]
    join-fileview-by-specimen:
      run: https://raw.githubusercontent.com/sgosline/synapse-workflow-cwl-tools/master/join-fileview-by-specimen-tool.cwl
      in:
        filelist: get-files/filepath
        values: get-samples-from-fv/names
        manifest_file: get-clinical/query_result
        parentid: parentid
        key: group_by
      out:
        [newmanifest]
    harmonize-counts:
      run: https://raw.githubusercontent.com/Sage-Bionetworks/rare-disease-workflows/main/rna-seq-workflow/steps/merge-counts-with-meta-tool.cwl
      in:
        synapse_config: synapse_config
        manifest:
           source: join-fileview-by-specimen/newmanifest
        files:
           source: get-files/filepath
      out:
        [merged]
    add-to-table:
      run: https://raw.githubusercontent.com/Sage-Bionetworks/rare-disease-workflows/main/synapse-table-store/synapse-table-store-tool.cwl
      in:
        synapse_config: synapse_config
        tableparentid: tableparentid
        tablename: tablename
        file: harmonize-counts/merged
      out:
        []
