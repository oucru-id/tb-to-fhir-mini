#!/usr/bin/env python3

import json
import argparse
import sys
import os
import base64
import uuid
from datetime import datetime, timezone
from clinical_metadata_parser import load_clinical_metadata, find_matching_sample, get_clinical_value

def create_organization_resource():
    return {
        "resourceType": "Organization",
        "id": "1234",
        "meta": {"profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Organization"]},
        "identifier": [{"use": "official", "system": "http://sys-ids.kemkes.go.id/organization", "value": "1234"}],
        "active": True,
        "type": [{"coding": [{"system": "http://terminology.kemkes.go.id/CodeSystem/organization-type", "code": "102", "display": "Pusat Kesehatan Masyarakat"}], "text": "Puskesmas"}],
        "name": "PUSKESMAS 1",
        "alias": ["1"],
        "telecom": [
            {"system": "phone", "value": "1234", "use": "work"},
            {"system": "email", "value": "1234@gmail.com", "use": "work"}
        ],
        "address": [{
            "use": "work", "type": "physical", "line": ["Jl. Raya"],
            "city": "KAB. PURBALINGGA", "state": "JAWA TENGAH", "country": "ID",
            "extension": [
                {"url": "https://fhir.kemkes.go.id/r4/StructureDefinition/administrativeCode", "extension": [
                    {"url": "province", "valueCode": "**"},
                    {"url": "city", "valueCode": "**"},
                    {"url": "district", "valueCode": "**"}
                ]},
                {"url": "http://hl7.org/fhir/StructureDefinition/geolocation", "extension": [
                    {"url": "latitude", "valueDecimal": -11.1111},
                    {"url": "longitude", "valueDecimal": 111.1111}
                ]}
            ]
        }]
    }

def create_practitioner_resource():
    return {
        "resourceType": "Practitioner",
        "id": "1234",
        "meta": {"profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Practitioner"]},
        "identifier": [{"use": "official", "system": "https://fhir.kemkes.go.id/id/nik", "value": "1234"}],
        "active": True,
        "name": [{"use": "official", "text": "Budi"}],
        "telecom": [{"system": "phone", "value": "+62-1234", "use": "work"}],
        "gender": "male",
        "birthDate": "YYYY-MM-DD",
        "qualification": [{
            "identifier": [{"system": "https://fhir.kemkes.go.id/id/str-kki-number", "value": "1234"}],
            "code": {"coding": [{"system": "https://terminology.kemkes.go.id/v1-0302", "code": "STR-KKI", "display": "Surat Tanda Registrasi Dokter"}], "text": "Surat Tanda Registrasi Dokter"},
            "period": {"start": "YYYY-MM-DD"}
        }]
    }

def create_practitioner_role_resource():
    return {
        "resourceType": "PractitionerRole",
        "id": "SPHERES-Nurse-Role",
        "meta": {"profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/PractitionerRole"]},
        "active": True,
        "practitioner": {"reference": "Practitioner/1234", "display": "Budi"},
        "organization": {"reference": "Organization/1234", "display": "PUSKESMAS 1"},
        "code": [{"coding": [{"system": "http://snomed.info/sct", "code": "224535009", "display": "Registered nurse"}], "text": "Perawat (Noisn Ners)"}],
        "telecom": [{"system": "phone", "value": "+62-1234", "use": "work"}]
    }

def create_patient_resource(sample_id, clinical_data):
    family_name = get_clinical_value(clinical_data, 'family_name')
    given_name = get_clinical_value(clinical_data, 'given_name')
    gender = get_clinical_value(clinical_data, 'gender', 'unknown').lower()
    birth_date = get_clinical_value(clinical_data, 'birth_date')
    nik = get_clinical_value(clinical_data, 'nik')
    
    if gender in ['laki-laki', 'pria', 'male', 'm']: gender = "male"
    elif gender in ['perempuan', 'wanita', 'female', 'f']: gender = "female"
    else: gender = "unknown"
    
    return {
        "resourceType": "Patient",
        "id": f"{sample_id}-patient",
        "meta": {"profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Patient"]},
        "active": True,
        "name": [{"use": "official", "family": family_name, "given": [given_name]}],
        "gender": gender,
        "birthDate": birth_date,
        "identifier": [
            {"use": "official", "system": "https://fhir.kemkes.go.id/id/nik", "value": nik},
            {"use": "usual", "type": {"coding": [{"system": "http://terminology.hl7.org/CodeSystem/v2-0203", "code": "MR", "display": "Medical record number"}]}, "system": "http://sys-ids.kemkes.go.id/mr/100007730", "value": sample_id}
        ],
        "extension": [
            {
                "url": "https://fhir.kemkes.go.id/r4/StructureDefinition/administrativeCode",
                "extension": [
                    {"url": "province", "valueCode": "**"},
                    {"url": "city", "valueCode": "**"},
                    {"url": "district", "valueCode": "**"},
                    {"url": "village", "valueCode": "**"}
                ]
            },
            {"url": "https://fhir.kemkes.go.id/r4/StructureDefinition/citizenshipStatus", "valueCode": "WNI"}
        ],
        "address": [{
            "use": "home", "type": "physical", "text": get_clinical_value(clinical_data, 'address'),
            "city": get_clinical_value(clinical_data, 'city'), "state": get_clinical_value(clinical_data, 'state'), "country": "ID",
            "extension": [
                {"url": "https://fhir.kemkes.go.id/r4/StructureDefinition/administrativeCode", "extension": [
                    {"url": "province", "valueCode": "**"},
                    {"url": "city", "valueCode": "**"},
                    {"url": "district", "valueCode": "**"},
                    {"url": "village", "valueCode": "**"}
                ]},
                {"url": "http://hl7.org/fhir/StructureDefinition/geolocation", "extension": [
                    {"url": "latitude", "valueDecimal": -1.1111},
                    {"url": "longitude", "valueDecimal": 111.1111}
                ]}
            ]
        }]
    }

def create_service_request_resource(sample_id, clinical_data):
    given_name = get_clinical_value(clinical_data, 'given_name', 'Unknown')
    family_name = get_clinical_value(clinical_data, 'family_name', 'Unknown')
    patient_display = f"{given_name} {family_name}"

    return {
        "resourceType": "ServiceRequest",
        "id": f"{sample_id}-service-request",
        "meta": {"profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/ServiceRequest"]},
        "identifier": [{"system": "http://sys-ids.kemkes.go.id/servicerequest/100007730", "value": f"SR-{sample_id}"}],
        "status": "active",
        "intent": "original-order",
        "priority": "routine",
        "category": [{"coding": [{"system": "http://snomed.info/sct", "code": "108252007", "display": "Laboratory procedure"}]}],
        "code": {"coding": [{"system": "http://loinc.org", "code": "69548-6", "display": "Genetic variant assessment"}], "text": "TB Genetic Variant Assessment"},
        "subject": {"reference": f"Patient/{sample_id}-patient", "display": patient_display},
        "occurrenceDateTime": datetime.now(timezone.utc).isoformat(),
        "requester": {"reference": "Practitioner/1234", "display": "Budi"},
        "performer": [{"reference": "PractitionerRole/SPHERES-Nurse-Role", "display": "Registered nurse"}]
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
                
                if 'clinical significance' in code_display:
                     significance = component.get('valueCodeableConcept', {}).get('text', '')

            if current_gene and significance and ('assoc w r' in significance.lower() or 'resistance' in significance.lower()):
                 if 'not associated' not in significance.lower():
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

def get_resistance_conclusion_coding(resistance_class):
    coding_map = {
        "Sensitive": {
           "system": "https://terminology.kemkes.go.id/CodeSystem/episodeofcare-type",
           "code": "TB-SO",
           "display": "Tuberkulosis Sensitif Obat"
        },
        "RR-TB": {"system": "http://snomed.info/sct", "code": "415345001", "display": "Rifampicin resistant tuberculosis"},
        "HR-TB": {"system": "http://snomed.info/sct", "code": "414546009", "display": "Isoniazid resistant tuberculosis"},
        "MDR-TB": {"system": "http://snomed.info/sct", "code": "423092005", "display": "Multidrug resistant tuberculosis"},
        "Pre-XDR-TB": {"system": "http://terminology.kemkes.go.id/CodeSystem/clinical-term", "code": "OV000435", "display": "Pre-XDR"},
        "XDR-TB": {"system": "http://snomed.info/sct", "code": "710106005", "display": "Extensively drug resistant tuberculosis"},
        "Streptomycin-resistant TB": {"system": "http://snomed.info/sct", "code": "415622003", "display": "Streptomycin resistant tuberculosis"},
        "Ethionamide-resistant TB": {"system": "http://snomed.info/sct", "code": "414149006", "display": "Ethionamide resistant tuberculosis"},
        "Pyrazinamide-resistant TB": {"system": "http://snomed.info/sct", "code": "415222009", "display": "Pyrazinamide resistant tuberculosis"},
        "Ciprofloxacin-resistant TB": {"system": "http://snomed.info/sct", "code": "413852006", "display": "Ciprofloxacin resistant tuberculosis"},
        "Ethambutol-resistant TB": {"system": "http://snomed.info/sct", "code": "414146004", "display": "Ethambutol resistant tuberculosis"},
        "Drug-resistant": {"system": "http://snomed.info/sct", "code": "413556004", "display": "Antibiotic resistant tuberculosis"}
    }
    return coding_map.get(resistance_class)

def create_diagnostic_report(sample_id, observations, clinical_data):
    patient_name = get_clinical_value(clinical_data, 'given_name') + " " + get_clinical_value(clinical_data, 'family_name')
    
    resistance_class, resistance_description, resistant_genes, resistant_drugs = classify_drug_resistance(observations)
    
    lineage_info = None
    for obs in observations:
        codes = obs.get('code', {}).get('coding', [])
        for code in codes:
             if code.get('code') == '614-8':
                 lineage_info = obs.get('valueCodeableConcept', {}).get('text')

    conclusion_parts = []
    if resistance_class == "Sensitive":
        conclusion_parts.append(resistance_class)
        conclusion_parts.append("No resistance-associated gene detected")
    else:
        conclusion_parts.append(f"{resistance_class} ({resistance_description})")
        if resistant_genes: conclusion_parts.append(f"Detected resistance genes: {', '.join(resistant_genes)}")
        if resistant_drugs: conclusion_parts.append(f"Detected drug resistance: {', '.join(resistant_drugs)}")
    
    if lineage_info: conclusion_parts.append(f"TB {lineage_info} detected")
    conclusion = ". ".join(conclusion_parts)

    conclusion_codes = []
    resistance_coding = get_resistance_conclusion_coding(resistance_class)
    
    display_text = resistance_class
    if resistance_class == "Sensitive":
        display_text = "Sensitive - No resistance detected"
        
    code_entry = {"text": display_text}
    if resistance_coding:
        code_entry["coding"] = [resistance_coding]
    
    conclusion_codes.append(code_entry)

    current_time = datetime.now(timezone.utc).isoformat()
    
    html_content = f"""<div xmlns="http://www.w3.org/1999/xhtml">
<h1>TB Genomic Analysis Report (Deeplex)</h1>
<p><strong>Patient:</strong> {patient_name}</p>
<p><strong>Sample ID:</strong> {sample_id}</p>
<p><strong>Report Date:</strong> {current_time}</p>
<p><strong>Resistance Classification:</strong> {resistance_class}</p>
<p><strong>Conclusion:</strong> {conclusion}</p>
"""
    if observations:
        html_content += "<h2>Detected Variants</h2><ul>"
        for obs in observations:
            if obs.get('code', {}).get('coding', [{}])[0].get('code') == '85502-2': continue
            
            gene = ""
            change = ""
            for comp in obs.get('component', []):
                c_code = comp.get('code', {}).get('coding', [{}])[0].get('code')
                if c_code == '48018-6': gene = comp.get('valueCodeableConcept', {}).get('text', '')
                if c_code == '81290-9': change = comp.get('valueCodeableConcept', {}).get('text', '')
            
            if gene or change:
                html_content += f"<li>{gene}: {change}</li>"
        html_content += "</ul>"
    
    if lineage_info: html_content += f"<p><strong>Mycobacterial Lineage:</strong> {lineage_info}</p>"
    html_content += "</div>"
    html_base64 = base64.b64encode(html_content.encode('utf-8')).decode('utf-8')

    return {
        "resourceType": "DiagnosticReport",
        "id": f"{sample_id}-deeplex-report",
        "meta": {
            "profile": ["http://hl7.org/fhir/uv/genomics-reporting/StructureDefinition/genomics-report"],
            "tag": [{"system": "http://terminology.kemkes.go.id/sp", "code": "genomics", "display": "Genomics"}]
        },
        "status": "final",
        "code": {"coding": [{"system": "http://loinc.org", "code": "81247-9", "display": "Master HL7 genetic variant reporting panel"}], "text": "TB Genomic Analysis Report"},
        "subject": {"reference": f"Patient/{sample_id}-patient", "display": patient_name},
        "issued": current_time,
        "performer": [
            {"reference": "Organization/1234", "display": "PUSKESMAS 1"},
            {"reference": "Practitioner/1234", "display": "Budi"}
        ],
        "result": [{"reference": f"Observation/{obs['id']}"} for obs in observations],
        "conclusion": conclusion,
        "conclusionCode": conclusion_codes,
        "presentedForm": [{
            "contentType": "text/html",
            "language": "en-US",
            "title": "TB Genomic Analysis Report",
            "data": html_base64
        }]
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True, help='Input FHIR bundle (Observations)')
    parser.add_argument('--output', required=True, help='Output merged FHIR bundle')
    parser.add_argument('--clinical_metadata', required=True, help='Clinical metadata CSV')
    args = parser.parse_args()

    clinical_data_map = load_clinical_metadata(args.clinical_metadata)
    with open(args.input, 'r') as f:
        input_bundle = json.load(f)
    
    observations = [entry['resource'] for entry in input_bundle.get('entry', []) if entry['resource']['resourceType'] == 'Observation']
    
    if not observations:
        print("No observations found.")
        sample_id = os.path.basename(args.input).split('.')[0]
    else:
        ref = observations[0].get('subject', {}).get('reference', '')
        if ref:
            sample_id = ref.replace('Patient/', '').replace('-patient', '')
        else:
            sample_id = os.path.basename(args.input).split('.')[0]

    sample_clinical = find_matching_sample(sample_id, clinical_data_map)
    if not sample_clinical:
        print(f"Warning: No clinical data found for {sample_id}")
        sample_clinical = {'given_name': 'Unknown', 'family_name': '', 'sample_id': sample_id}

    org = create_organization_resource()
    practitioner = create_practitioner_resource()
    role = create_practitioner_role_resource()
    patient = create_patient_resource(sample_id, sample_clinical)
    service_request = create_service_request_resource(sample_id, sample_clinical)
    report = create_diagnostic_report(sample_id, observations, sample_clinical)
    
    entries = []
    entries.append({"fullUrl": f"urn:uuid:{org['id']}", "resource": org, "request": {"method": "PUT", "url": f"Organization/{org['id']}"}})
    entries.append({"fullUrl": f"urn:uuid:{practitioner['id']}", "resource": practitioner, "request": {"method": "PUT", "url": f"Practitioner/{practitioner['id']}"}})
    entries.append({"fullUrl": f"urn:uuid:{role['id']}", "resource": role, "request": {"method": "PUT", "url": f"PractitionerRole/{role['id']}"}})
    
    entries.append({"fullUrl": f"urn:uuid:{patient['id']}", "resource": patient, "request": {"method": "PUT", "url": f"Patient/{patient['id']}"}})
    entries.append({"fullUrl": f"urn:uuid:{service_request['id']}", "resource": service_request, "request": {"method": "PUT", "url": f"ServiceRequest/{service_request['id']}"}})
    entries.append({"fullUrl": f"urn:uuid:{report['id']}", "resource": report, "request": {"method": "PUT", "url": f"DiagnosticReport/{report['id']}"}})
    
    for obs in observations:
        entries.append({"fullUrl": f"urn:uuid:{obs['id']}", "resource": obs, "request": {"method": "PUT", "url": f"Observation/{obs['id']}"}})

    output_bundle = {
        "resourceType": "Bundle",
        "id": str(uuid.uuid4()),
        "meta": {"lastUpdated": datetime.now(timezone.utc).isoformat(), "profile": ["https://fhir.kemkes.go.id/r4/StructureDefinition/Bundle"]},
        "type": "transaction",
        "entry": entries
    }

    with open(args.output, 'w') as f:
        json.dump(output_bundle, f, indent=2)

if __name__ == "__main__":
    main()