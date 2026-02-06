#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

log.info """
    Mycobacterium tuberculosis Mutation Analysis Mini Pipeline (v${params.version})
    VCF and Deeplex Excel Conversion to FHIR
    Developed by SPHERES-OUCRU ID
    Documentation: https://docs.google.com/document/d/1loZhheM22cWnU3taqAaef16lePK7BoNJ5lTjEuXwFYE/edit?usp=sharing
"""

include { VCF_PROCESSING }        from './workflows/vcf.nf'  
include { LINEAGE }               from './workflows/lineage.nf'
include { GENERATE_SAMPLE_REPORTS } from './workflows/report.nf'
include { FHIR }                  from './workflows/fhir.nf'
include { VALIDATE }              from './workflows/validate_fhir.nf'
include { MERGE_CLINICAL_DATA }   from './workflows/merge_clinical_data.nf'
include { UPLOAD_FHIR }           from './workflows/upload_fhir.nf'
include { VERSIONS }              from './workflows/utils.nf'
include { DEEPLEX }               from './workflows/deeplex.nf'

workflow {
    // VCF input channel
    vcf_ch = Channel
        .fromPath("${params.vcf_dir}/*.vcf{,.gz}", checkIfExists: false)
        .map { file -> tuple(file.baseName.replaceFirst(/\.vcf(\.gz)?$/, ''), file) }

    // Deeplex Excel input channel
    deeplex_ch = Channel
        .fromPath("${params.deeplex_dir}/*.xlsx", checkIfExists: false)
        .filter { file -> !file.name.startsWith('~$') }

    vcf_out = VCF_PROCESSING(vcf_ch) 

    lineage_out = LINEAGE(vcf_out.filtered)

    all_annotated = vcf_out.annotated
        .map { it -> it instanceof List ? it[1] : it }
    
    lineage_files = lineage_out.lineage_results
        .map { sample_id, file_path -> file_path } 
        .collect()
    
    sample_reports = GENERATE_SAMPLE_REPORTS(all_annotated, lineage_files)
    
    fhir_out = FHIR(all_annotated, lineage_files) 
    
    clinical_metadata_ch = Channel.fromPath(params.clinical_metadata, checkIfExists: false)
        .first() 
    
    deeplex_clinical_ch = Channel.fromPath(params.clinical_metadata_deeplex, checkIfExists: false)
        .first()

    deeplex_out = DEEPLEX(deeplex_ch, deeplex_clinical_ch)

    merged_clinical_out = MERGE_CLINICAL_DATA(
        fhir_out.fhir_output, 
        clinical_metadata_ch
    )
    
    validation_out = VALIDATE(merged_clinical_out.merged_fhir)
    
    // Optional: Upload validated FHIR
    //upload_out = UPLOAD_FHIR(validation_out.validated_fhir)

    VERSIONS()
}