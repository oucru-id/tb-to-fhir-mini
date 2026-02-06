#!/usr/bin/env python3

import json
import argparse
import uuid
import sys
import os  
from datetime import datetime, timezone
from clinical_metadata_parser import load_clinical_metadata, find_matching_sample, get_clinical_value
import base64
import re

def debug_print(message):
    print(f"DEBUG: {message}", file=sys.stderr)

def create_patient_resource(sample_id, clinical_data=None):
    if not clinical_data:
        raise ValueError(f"Clinical data is required for sample {sample_id}")
    
    family_name = get_clinical_value(clinical_data, 'family_name')
    given_name = get_clinical_value(clinical_data, 'given_name')
    gender = get_clinical_value(clinical_data, 'gender', 'unknown').lower()
    birth_date = get_clinical_value(clinical_data, 'birth_date')
    nik = get_clinical_value(clinical_data, 'nik')
    
    if gender in ['laki-laki', 'pria', 'male', 'm']:
        gender = "male"
    elif gender in ['perempuan', 'wanita', 'female', 'f']:
        gender = "female"
    else:
        gender = "unknown"
    
    return {
        "resourceType": "Patient",
        "id": f"{sample_id}-patient",
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Patient"]
        },
        "active": True,
        "name": [
            {
                "use": "official",
                "family": family_name,
                "given": [given_name]
            }
        ],
        "gender": gender,
        "birthDate": birth_date,
        "identifier": [
            {
                "use": "official",
                "system": "https://fhir.kemkes.go.id/id/nik",
                "value": nik
            },
            {
                "use": "usual",
                "type": {
                    "coding": [
                        {
                            "system": "http://terminology.hl7.org/CodeSystem/v2-0203",
                            "code": "MR",
                            "display": "Medical record number"
                        }
                    ]
                },
                "system": "http://sys-ids.kemkes.go.id/mr/100007730",
                "value": sample_id
            }
        ],
        "extension": [
            {
                "url": "https://fhir.kemkes.go.id/r4/StructureDefinition/administrativeCode",
                "extension": [
                    {
                        "url": "province",
                        "valueCode": "**"  
                    },
                    {
                        "url": "city",
                        "valueCode": "**"  
                    },
                    {
                        "url": "district",
                        "valueCode": "**"  
                    },
                    {
                        "url": "village",
                        "valueCode": "**"  
                    }
                ]
            },
            {
                "url": "https://fhir.kemkes.go.id/r4/StructureDefinition/citizenshipStatus",
                "valueCode": "WNI" 
            }
        ],
        "address": [
            {
                "use": "home", 
                "type": "physical", 
                "text": get_clinical_value(clinical_data, 'address'), 
                "city": get_clinical_value(clinical_data, 'city'),  
                "state": get_clinical_value(clinical_data, 'state'), 
                "country": "ID",  
                "extension": [
                    {
                        "url": "https://fhir.kemkes.go.id/r4/StructureDefinition/administrativeCode",
                        "extension": [
                            {
                                "url": "province",
                                "valueCode": "**"
                            },
                            {
                                "url": "city",
                                "valueCode": "**"
                            },
                            {
                                "url": "district", 
                                "valueCode": "**"
                            },
                            {
                                "url": "village",
                                "valueCode": "**"
                            }
                        ]
                    },
                    {
                        "url": "http://hl7.org/fhir/StructureDefinition/geolocation",
                        "extension": [
                            {
                                "url": "latitude",
                                "valueDecimal": -1.1111 
                            },
                            {
                                "url": "longitude", 
                                "valueDecimal": 111.1111
                            }
                        ]
                    }
                ]
            }
        ]
    }

def create_specimen_resource(sample_id, clinical_data=None):
    if clinical_data:
        given_name = get_clinical_value(clinical_data, 'given_name', 'Unknown')
        family_name = get_clinical_value(clinical_data, 'family_name', 'Unknown')
        patient_display = f"{given_name} {family_name}"
    else:
        patient_display = f"Patient {sample_id}"
    
    return {
        "resourceType": "Specimen",
        "id": f"{sample_id}-specimen",
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Specimen"]
        },
        "identifier": [
            {
                "system": "http://sys-ids.kemkes.go.id/specimen/100007730",
                "value": f"SPEC-{sample_id}"
            }
        ],
        "status": "available",
        "subject": {
            "reference": f"Patient/{sample_id}-patient",
            "display": patient_display
        },
        "receivedTime": datetime.now(timezone.utc).isoformat(),
        "collection": {
            "collectedDateTime": datetime.now(timezone.utc).isoformat(),
            "collector": {
                "reference": "Practitioner/1234",
                "display": "Budi"
            },
            "method": {
                "coding": [
                    {
                        "system": "http://terminology.hl7.org/CodeSystem/v2-0488",
                        "code": "KOFFP",
                        "display": "Plate, Cough"
                    }
                ],
                "text": "Sputum collection"
            },
            "quantity": {
                "value": 2,
                "unit": "mL",
                "system": "http://unitsofmeasure.org",
                "code": "mL"
            }
        },
        "type": {
            "coding": [
                {
                    "system": "http://snomed.info/sct",
                    "code": "119334006",
                    "display": "Sputum specimen"
                }
            ],
            "text": "Sputum specimen for TB testing"
        },
        "note": [
            {
                "text": f"Collected sputum sample from {patient_display} ({sample_id}) for TB genetic testing purposes"
            }
        ]
    }

def create_organization_resource():
    return {
        "resourceType": "Organization",
        "id": "1234",
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Organization"]
        },
        "identifier": [
            {
                "use": "official",
                "system": "http://sys-ids.kemkes.go.id/organization",
                "value": "1234"
            }
        ],
        "active": True,
        "type": [
            {
                "coding": [
                    {
                        "system": "http://terminology.kemkes.go.id/CodeSystem/organization-type",
                        "code": "102",
                        "display": "Pusat Kesehatan Masyarakat"
                    }
                ],
                "text": "Puskesmas"
            }
        ],
        "name": "PUSKESMAS 1",
        "alias": ["1"],
        "telecom": [
            {
                "system": "phone",
                "value": "1234",
                "use": "work"
            },
            {
                "system": "email",
                "value": "1234@gmail.com",
                "use": "work"
            }
        ],
        "address": [
            {
                "use": "work",
                "type": "physical",
                "line": ["Jl. Raya"],
                "city": "KAB. PURBALINGGA",
                "state": "JAWA TENGAH",
                "country": "ID",
                "extension": [
                    {
                        "url": "https://fhir.kemkes.go.id/r4/StructureDefinition/administrativeCode",
                        "extension": [
                            {
                                "url": "province",
                                "valueCode": "**"
                            },
                            {
                                "url": "city",
                                "valueCode": "**"
                            },
                            {
                                "url": "district",
                                "valueCode": "**"
                            }
                        ]
                    },
                    {
                        "url": "http://hl7.org/fhir/StructureDefinition/geolocation",
                        "extension": [
                            {
                                "url": "latitude",
                                "valueDecimal": -11.1111 
                            },
                            {
                                "url": "longitude",
                                "valueDecimal": 111.1111
                            }
                        ]
                    }
                ]
            }
        ]
    }

def create_practitioner_resource():
    return {
        "resourceType": "Practitioner",
        "id": "1234",
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Practitioner"]
        },
        "identifier": [
            {
                "use": "official",
                "system": "https://fhir.kemkes.go.id/id/nik",
                "value": "1234"
            }
        ],
        "active": True,
        "name": [
            {
                "use": "official",
                "text": "Budi"
            }
        ],
        "telecom": [
            {
                "system": "phone",
                "value": "+62-1234",
                "use": "work"
            }
        ],
        "gender": "male",
        "birthDate": "YYYY-MM-DD",
        "qualification": [
            {
                "identifier": [
                    {
                        "system": "https://fhir.kemkes.go.id/id/str-kki-number",
                        "value": "1234"
                    }
                ],
                "code": {
                    "coding": [
                        {
                            "system": "https://terminology.kemkes.go.id/v1-0302",
                            "code": "STR-KKI",
                            "display": "Surat Tanda Registrasi Dokter"
                        }
                    ],
                    "text": "Surat Tanda Registrasi Dokter"
                },
                "period": {
                    "start": "YYYY-MM-DD"
                }
            }
        ]
    }

def create_practitioner_role_resource():
    return {
        "resourceType": "PractitionerRole",
        "id": "SPHERES-Nurse-Role",
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/PractitionerRole"]
        },
        "active": True,
        "practitioner": {
            "reference": "Practitioner/1234",
            "display": "Budi"
        },
        "organization": {
            "reference": "Organization/1234",
            "display": "PUSKESMAS 1"
        },
        "code": [
            {
                "coding": [
                    {
                        "system": "http://snomed.info/sct",
                        "code": "224535009",
                        "display": "Registered nurse"
                    }
                ],
                "text": "Perawat (Noisn Ners)"
            }
        ],
        "telecom": [
            {
                "system": "phone",
                "value": "+62-1234",
                "use": "work"
            }
        ]
    }

def classify_drug_resistance(observations):
    resistant_drugs_groups = set()
    detected_drugs = set()
    detected_genes = set()
    
    for obs in observations:
        codes = obs.get('code', {}).get('coding', [])
        is_panel = any(c.get('code') == '89486-5' for c in codes)
        
        if is_panel:
            components = obs.get('component', [])
            for comp in components:
                value_coding = comp.get('valueCodeableConcept', {}).get('coding', [])
                is_resistant = any(vc.get('code') == 'LA6676-6' for vc in value_coding)
                
                if is_resistant:
                    drug_display = comp.get('code', {}).get('coding', [{}])[0].get('display', '').lower()
                    
                    detected_drugs.add(drug_display.replace('[susceptibility]', '').strip())
                    
                    if 'rifampicin' in drug_display or 'rifampin' in drug_display:
                        resistant_drugs_groups.add('rifampicin')
                    elif 'isoniazid' in drug_display:
                        resistant_drugs_groups.add('isoniazid')
                    elif 'ethambutol' in drug_display:
                        resistant_drugs_groups.add('ethambutol')
                    elif 'pyrazinamide' in drug_display:
                        resistant_drugs_groups.add('pyrazinamide')
                    elif 'streptomycin' in drug_display:
                        resistant_drugs_groups.add('streptomycin')
                    elif 'ethionamide' in drug_display:
                        resistant_drugs_groups.add('ethionamide')
                    elif any(fq in drug_display for fq in ['levofloxacin', 'moxifloxacin', 'ofloxacin', 'ciprofloxacin']):
                        resistant_drugs_groups.add('fluoroquinolone')
                    elif any(sli in drug_display for sli in ['amikacin', 'kanamycin', 'capreomycin']):
                        resistant_drugs_groups.add('second_line_injectable')
                    elif any(ga in drug_display for ga in ['bedaquiline', 'linezolid']):
                        resistant_drugs_groups.add('group_a')
                    else:
                        resistant_drugs_groups.add(drug_display)

    for obs in observations:
        codes = obs.get('code', {}).get('coding', [])
        is_variant = any(c.get('code') == '69548-6' for c in codes)
        
        if is_variant:
            components = obs.get('component', [])
            
            significance = ""
            current_gene = ""
            
            for component in components:
                code_display = component.get('code', {}).get('coding', [{}])[0].get('display', '').lower()
                
                if 'gene studied' in code_display:
                    current_gene = component.get('valueCodeableConcept', {}).get('text', '')
                
                if 'genetic variation clinical significance' in code_display:
                    significance = component.get('valueCodeableConcept', {}).get('text', '')

            if current_gene and significance and "Assoc w R" in significance:
                detected_genes.add(current_gene)
    
    has_rif = 'rifampicin' in resistant_drugs_groups
    has_inh = 'isoniazid' in resistant_drugs_groups
    has_fq = 'fluoroquinolone' in resistant_drugs_groups
    has_group_a = 'group_a' in resistant_drugs_groups
    
    is_mdr = has_rif and has_inh
    is_rr = has_rif
    
    classification = "Sensitive"
    description = "No resistance mutations detected"

    if not resistant_drugs_groups:
        classification = "Sensitive"
        description = "No resistance mutations detected"
    elif (is_mdr or is_rr) and has_fq and has_group_a:
        classification = "XDR-TB"
        description = "Extensively drug-resistant tuberculosis (MDR/RR + FQ + Group A)"
    elif (is_mdr or is_rr) and has_fq:
        classification = "Pre-XDR-TB"
        description = "Pre-extensively drug-resistant tuberculosis (MDR/RR + FQ)"
    elif has_rif and has_inh:
        classification = "MDR-TB"
        description = "Multidrug-resistant tuberculosis"
    elif has_rif and not has_inh:
        classification = "RR-TB"
        description = "Rifampicin-resistant tuberculosis"
    elif has_inh and not has_rif:
        classification = "HR-TB"
        description = "Isoniazid-resistant tuberculosis"
    elif len(resistant_drugs_groups) == 1:
        if 'streptomycin' in resistant_drugs_groups:
            classification = "Streptomycin-resistant TB"
            description = "Streptomycin mono-resistant tuberculosis"
        elif 'ethionamide' in resistant_drugs_groups:
            classification = "Ethionamide-resistant TB"
            description = "Ethionamide mono-resistant tuberculosis"
        elif 'pyrazinamide' in resistant_drugs_groups:
            classification = "Pyrazinamide-resistant TB"
            description = "Pyrazinamide mono-resistant tuberculosis"
        elif 'ethambutol' in resistant_drugs_groups:
            classification = "Ethambutol-resistant TB"
            description = "Ethambutol mono-resistant tuberculosis"
        elif 'fluoroquinolone' in resistant_drugs_groups:
            has_cipro = any('ciprofloxacin' in d for d in detected_drugs)
            has_other_fq = any(fq in d for d in detected_drugs for fq in ['levofloxacin', 'moxifloxacin', 'ofloxacin'])
            
            if has_cipro and not has_other_fq:
                classification = "Ciprofloxacin-resistant TB"
                description = "Ciprofloxacin mono-resistant tuberculosis"
            else:
                classification = "Drug-resistant"
                description = f"Resistance to: {', '.join(sorted(resistant_drugs_groups))}"
        else:
            classification = "Drug-resistant"
            description = f"Resistance to: {', '.join(sorted(resistant_drugs_groups))}"
    else:
        classification = "Drug-resistant"
        description = f"Resistance to: {', '.join(sorted(resistant_drugs_groups))}"

    return classification, description, sorted(list(detected_genes)), sorted(list(detected_drugs))

def extract_lineage_info(observations):
    for obs in observations:
        codes = obs.get('code', {}).get('coding', [])
        for code in codes:
            if code.get('code') == '614-8':
                lineage_text = obs.get('valueCodeableConcept', {}).get('text', '')
                if lineage_text and lineage_text != 'unknown':
                    return lineage_text

        components = obs.get('component', [])
        for component in components:
            code_display = component.get('code', {}).get('coding', [{}])[0].get('display', '').lower()
            if 'mycobacterial strain' in code_display or 'lineage' in code_display:
                lineage_text = component.get('valueCodeableConcept', {}).get('text', '')
                if lineage_text and lineage_text != 'unknown':
                    return lineage_text
    return None

def get_resistance_conclusion_coding(resistance_class):
    coding_map = {
        "Sensitive": {
           "system": "https://terminology.kemkes.go.id/CodeSystem/episodeofcare-type",
           "code": "TB-SO",
           "display": "Tuberkulosis Sensitif Obat"
        },
        "RR-TB": {
            "system": "http://snomed.info/sct",
            "code": "415345001",
            "display": "Rifampicin resistant tuberculosis"
        },
        "HR-TB": {
            "system": "http://snomed.info/sct", 
            "code": "414546009",
            "display": "Isoniazid resistant tuberculosis"
        },
        "MDR-TB": {
            "system": "http://snomed.info/sct",
            "code": "423092005", 
            "display": "Multidrug resistant tuberculosis"
        },
        "Pre-XDR-TB": {
            "system": "http://terminology.kemkes.go.id/CodeSystem/clinical-term",
            "code": "OV000435",
            "display": "Pre-XDR"
        },
        "XDR-TB": {
            "system": "http://snomed.info/sct",
            "code": "710106005",
            "display": "Extensively drug resistant tuberculosis"
        },
        "Streptomycin-resistant TB": {
            "system": "http://snomed.info/sct",
            "code": "415622003",
            "display": "Streptomycin resistant tuberculosis"
        },
        "Ethionamide-resistant TB": {
            "system": "http://snomed.info/sct",
            "code": "414149006",
            "display": "Ethionamide resistant tuberculosis"
        },
        "Pyrazinamide-resistant TB": {
            "system": "http://snomed.info/sct",
            "code": "415222009",
            "display": "Pyrazinamide resistant tuberculosis"
        },
        "Ciprofloxacin-resistant TB": {
            "system": "http://snomed.info/sct",
            "code": "413852006",
            "display": "Ciprofloxacin resistant tuberculosis"
        },
        "Ethambutol-resistant TB": {
            "system": "http://snomed.info/sct",
            "code": "414146004",
            "display": "Ethambutol resistant tuberculosis"
        },
        "Drug-resistant": {
            "system": "http://snomed.info/sct",
            "code": "413556004",
            "display": "Antibiotic resistant tuberculosis"
        }
    }
    
    return coding_map.get(resistance_class)

def create_diagnostic_report(sample_id, observations, clinical_data=None):
    
    resistance_class, resistance_description, resistant_genes, resistant_drugs = classify_drug_resistance(observations)
    
    lineage_info = extract_lineage_info(observations)
    
    conclusion_parts = []
    
    if resistance_class == "Sensitive":
        conclusion_parts.append(resistance_class)
        conclusion_parts.append("No resistance-associated gene detected")
        conclusion_parts.append("No drug resistance detected")
    else:
        conclusion_parts.append(f"{resistance_class} ({resistance_description})")
        if resistant_genes:
            conclusion_parts.append(f"Detected resistance genes: {', '.join(resistant_genes)}")
        if resistant_drugs:
            conclusion_parts.append(f"Detected drug resistance: {', '.join(resistant_drugs)}")
    
    if lineage_info:
        conclusion_parts.append(f"TB {lineage_info} detected")
    
    conclusion_parts.append("Reference genome: NC_000962.3")
    
    conclusion = ". ".join(conclusion_parts)
    
    conclusion_codes = []
    
    if resistance_class:
        resistance_coding = get_resistance_conclusion_coding(resistance_class)
        
        display_text = resistance_class
        if resistance_class == "Sensitive":
            display_text = "Sensitive - No resistance detected"
            
        code_entry = {"text": display_text}
        
        if resistance_coding:
            code_entry["coding"] = [resistance_coding]
            
        conclusion_codes.append(code_entry)
    
    if lineage_info:
        conclusion_codes.append({
            "text": lineage_info
        })
    
    if clinical_data:
        given_name = get_clinical_value(clinical_data, 'given_name', 'Unknown')
        family_name = get_clinical_value(clinical_data, 'family_name', 'Unknown')
        patient_display = f"{given_name} {family_name}"
    else:
        patient_display = f"Patient {sample_id}"
    
    report_id = f"{sample_id}-genomic-report"
    current_time = datetime.now(timezone.utc).isoformat()
    
    html_content = f"""<div xmlns="http://www.w3.org/1999/xhtml">
<h1>TB Genomic Analysis Report</h1>
<p><strong>Patient:</strong> {patient_display}</p>
<p><strong>Sample ID:</strong> {sample_id}</p>
<p><strong>Report Date:</strong> {current_time}</p>
<p><strong>Resistance Classification:</strong> {resistance_class}</p>
<p><strong>Conclusion:</strong> {conclusion}</p>
"""
    
    if observations:
        html_content += "<h2>Detected Mutations</h2><ul>"
        for obs in observations:
            components = obs.get('component', [])
            for component in components:
                if 'valueCodeableConcept' in component:
                    text = component['valueCodeableConcept'].get('text', '')
                    if text and text != 'unknown':
                        html_content += f"<li>{text}</li>"
        html_content += "</ul>"
    
    if lineage_info:
        html_content += f"<p><strong>Mycobacterial Lineage:</strong> {lineage_info}</p>"
    
    html_content += "</div>"
    
    html_base64 = base64.b64encode(html_content.encode('utf-8')).decode('utf-8')
    
    return {
        "resourceType": "DiagnosticReport",
        "id": report_id,
        "meta": {
            "profile": ["http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genomics-report"],
            "tag": [
                {
                    "system": "http://terminology.kemkes.go.id/sp",
                    "code": "genomics",
                    "display": "Genomics"
                }
            ]
        },
        "identifier": [
            {
                "system": "http://sys-ids.kemkes.go.id/diagnostic-report/100007730",
                "value": f"TB-GEN-{sample_id}-{datetime.now().strftime('%Y%m%d')}"
            }
        ],
        "status": "final",
        "category": [
            {
                "coding": [{
                    "system": "http://terminology.hl7.org/CodeSystem/v2-0074",
                    "code": "GE",
                    "display": "Genetics"
                }]
            }
        ],
        "code": {
            "coding": [{
                "system": "http://loinc.org",
                "code": "81247-9",
                "display": "Master HL7 genetic variant reporting panel"
            }],
            "text": "TB Genomic Analysis Report"
        },
        "subject": {
            "reference": f"Patient/{sample_id}-patient",
            "display": patient_display
        },
        "encounter": {
            "reference": f"Encounter/{sample_id}-encounter",
            "display": "TB Testing Encounter"
        },
        "effectiveDateTime": current_time,
        "issued": current_time,
        "performer": [
            {
                "reference": "Organization/1234",
                "display": "PUSKESMAS 1"
            },
            {
                "reference": "Practitioner/1234",
                "display": "Budi"
            }
        ],
        "result": [{"reference": f"Observation/{obs['id']}"} for obs in observations if obs.get('id')],
        "specimen": [{
            "reference": f"Specimen/{sample_id}-specimen",
            "display": f"Sputum specimen from {patient_display}"
        }],
        "conclusion": conclusion,
        "conclusionCode": conclusion_codes,
        "presentedForm": [
            {
                "contentType": "text/html",
                "language": "en-US", 
                "title": "TB Genomic Analysis Report",
                "data": html_base64
            }
        ]
    }

def create_service_request_resource(sample_id, clinical_data=None):
    if clinical_data:
        given_name = get_clinical_value(clinical_data, 'given_name', 'Unknown')
        family_name = get_clinical_value(clinical_data, 'family_name', 'Unknown')
        patient_display = f"{given_name} {family_name}"
    else:
        patient_display = f"Patient {sample_id}"

    return {
        "resourceType": "ServiceRequest",
        "id": f"{sample_id}-service-request",
        "meta": {
            "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/ServiceRequest"]
        },
        "identifier": [
            {
                "system": "http://sys-ids.kemkes.go.id/servicerequest/100007730",
                "value": f"SR-{sample_id}"
            }
        ],
        "status": "active",
        "intent": "original-order",
        "priority": "routine",
        "category": [
            {
                "coding": [
                    {
                        "system": "http://snomed.info/sct",
                        "code": "108252007",
                        "display": "Laboratory procedure"
                    }
                ]
            }
        ],
        "code": {
            "coding": [
                {
                    "system": "http://loinc.org",
                    "code": "69548-6",
                    "display": "Genetic variant assessment"
                }
            ],
            "text": "TB Genetic Variant Assessment"
        },
        "subject": {
            "reference": f"Patient/{sample_id}-patient",
            "display": patient_display
        },
        "encounter": {
            "reference": f"Encounter/{sample_id}-encounter",
            "display": "TB Testing Encounter"
        },
        "occurrenceDateTime": datetime.now(timezone.utc).isoformat(),
        "requester": {
            "reference": "Practitioner/1234",
            "display": "Budi"
        },
        "performer": [
            {
                "reference": "PractitionerRole/SPHERES-Nurse-Role",
                "display": "Registered nurse"
            }
        ]
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Path to input FHIR bundle')
    parser.add_argument('--output', required=True, help='Path to output merged FHIR bundle')
    parser.add_argument('--clinical_metadata', help='Path to clinical metadata CSV/Excel file')
    args = parser.parse_args()

    debug_print(f"Input FHIR file: {args.input}")
    debug_print(f"Clinical metadata file: {args.clinical_metadata}")
    debug_print(f"Output file: {args.output}")

    clinical_data = {}
    if args.clinical_metadata and os.path.exists(args.clinical_metadata):
        clinical_data = load_clinical_metadata(args.clinical_metadata)
    else:
        debug_print(f"Clinical metadata file not found or not provided: {args.clinical_metadata}")

    try:
        with open(args.input, 'r') as f:
            fhir_bundle = json.load(f)

        sample_ids = set()
        all_observations = []
        
        for entry in fhir_bundle.get('entry', []):
            resource = entry.get('resource', {})
            if resource.get('resourceType') == 'Observation':
                all_observations.append(resource)
                subject_ref = resource.get('subject', {}).get('reference', '')
                
                if subject_ref.startswith('Patient/'):
                    sample_id = subject_ref.replace('Patient/', '').replace('-patient', '')
                    sample_ids.add(sample_id)

        if not sample_ids or all(sid.startswith('NC-') for sid in sample_ids):
            filename = os.path.basename(args.input)
            filename_sample_id = filename.replace('.fhir.json', '').replace('_ont', '').replace('_illumina', '')
            sample_ids.add(filename_sample_id)

        matched_samples = {}
        for sample_id in sample_ids:
            sample_clinical_data = find_matching_sample(sample_id, clinical_data)
            if sample_clinical_data:
                matched_samples[sample_id] = sample_clinical_data

        merged_bundle = {
            "resourceType": "Bundle",
            "id": str(uuid.uuid4()),
            "meta": {
                "lastUpdated": datetime.now(timezone.utc).isoformat(),
                "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Bundle"]
            },
            "type": "transaction", 
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "entry": []
        }

        org_resource = create_organization_resource()
        merged_bundle['entry'].append({
            "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
            "resource": org_resource,
            "request": {
                "method": "PUT",
                "url": f"Organization/{org_resource['id']}"
            }
        })

        practitioner_resource = create_practitioner_resource()
        merged_bundle['entry'].append({
            "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
            "resource": practitioner_resource,
            "request": {
                "method": "PUT",
                "url": f"Practitioner/{practitioner_resource['id']}"
            }
        })

        role_resource = create_practitioner_role_resource()
        merged_bundle['entry'].append({
            "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
            "resource": role_resource,
            "request": {
                "method": "PUT",
                "url": f"PractitionerRole/{role_resource['id']}"
            }
        })

        for sample_id, sample_clinical_data in matched_samples.items():
            debug_print(f"Adding patient for sample: {sample_id}")
            
            patient_resource = create_patient_resource(sample_id, sample_clinical_data)
            merged_bundle['entry'].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
                "resource": patient_resource,
                "request": {
                    "method": "PUT",
                    "url": f"Patient/{patient_resource['id']}"
                }
            })

            service_request_resource = create_service_request_resource(sample_id, sample_clinical_data)
            merged_bundle['entry'].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
                "resource": service_request_resource,
                "request": {
                    "method": "PUT",
                    "url": f"ServiceRequest/{service_request_resource['id']}"
                }
            })

        observations_by_sample = {}
        for obs in all_observations:
            subject_ref = obs.get('subject', {}).get('reference', '')
            if subject_ref.startswith('Patient/'):
                sample_id = subject_ref.replace('Patient/', '').replace('-patient', '')
                if sample_id not in observations_by_sample:
                    observations_by_sample[sample_id] = []
                observations_by_sample[sample_id].append(obs)

        for sample_id, sample_observations in observations_by_sample.items():
            sample_clinical_data = matched_samples.get(sample_id)
            
            diagnostic_report = create_diagnostic_report(
                sample_id, 
                sample_observations, 
                sample_clinical_data
            )
            
            merged_bundle['entry'].append({
                "fullUrl": f"urn:uuid:{str(uuid.uuid4())}",
                "resource": diagnostic_report,
                "request": {
                    "method": "PUT",
                    "url": f"DiagnosticReport/{diagnostic_report['id']}"
                }
            })

        for entry in fhir_bundle.get('entry', []):
            resource = entry.get('resource', {})
            resource_type = resource.get('resourceType')
            resource_id = resource.get('id')
            
            entry_with_request = {
                "fullUrl": entry.get('fullUrl', f"urn:uuid:{str(uuid.uuid4())}"),
                "resource": resource,
                "request": {
                    "method": "PUT",
                    "url": f"{resource_type}/{resource_id}" if resource_id else f"{resource_type}"
                }
            }
            merged_bundle['entry'].append(entry_with_request)

        with open(args.output, 'w') as f:
            json.dump(merged_bundle, f, indent=2)

    except Exception as e:
        debug_print(f"Error occurred: {str(e)}")
        import traceback
        debug_print(f"Traceback: {traceback.format_exc()}")
        sys.exit(1)

if __name__ == "__main__":
    main()