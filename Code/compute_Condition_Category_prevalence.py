
from utils import *
import argparse
import glob 
import networkx as nx
import pandas as pd

mappings_dir = '../mappings/'
input_dir = '../../../BL_Work/openPrescribe/serialized/'
output_dir = '../data_prep/'


def DrugMatching(conditions, isCat=False):
    drug_association_graph = nx.read_gexf(mappings_dir + 'drug_association_graph.gexf')
    drug_cat_association_graph = nx.read_gexf(mappings_dir + 'category_association_graph.gexf')
    chem = pd.read_csv(mappings_dir + 'CHEM_MASTER_MAP.csv')

    DiseaseDrugs = {}
    for d in conditions:
        if not isCat:
            drugs = findDrugsForDisease(drug_association_graph,d ,chem)
        else:
            drugs = findDrugsForCategory(drug_cat_association_graph,d ,chem)

        for drug in drugs:
            DiseaseDrugs[drug] = {}
            DiseaseDrugs[drug]['chemName'] = drugs[drug]['name']
            DiseaseDrugs[drug]['disease'] = d
    
    disease_drug_map = {}
    for k in DiseaseDrugs:
        if DiseaseDrugs[k]['disease'] not in disease_drug_map:
            disease_drug_map[DiseaseDrugs[k]['disease']] = []
        disease_drug_map[DiseaseDrugs[k]['disease']].append(k)
    return DiseaseDrugs , disease_drug_map

def dumpDrugs(DiseaseDrugs):
    drug_map_dict = {'BNF_code':[], 'Drug_name':[] , 'Mapped_Condition': []}
    for k in DiseaseDrugs:
        drug_map_dict['BNF_code'].append(k)
        drug_map_dict['Drug_name'].append(DiseaseDrugs[k]['chemName'])
        drug_map_dict['Mapped_Condition'].append(DiseaseDrugs[k]['disease'])
    drug_map_df = pd.DataFrame.from_dict(drug_map_dict)    
    drug_map_df.to_csv(output_dir + 'Drugs.csv',index=False)

def loadLSOA_mappings():
    LSOA_dist_old = json.load(open(mappings_dir + 'GP_LSOA_PATIENTSDIST.json','rb'))
    LSOA_dist_2021 = json.load(open(mappings_dir + 'GP_LSOA_PATIENTSDIST_2021.json','rb'))
    return [LSOA_dist_old , LSOA_dist_2021]

def prepare_lsoa_GP_population():
    LSOA_patient_pop = {}
    LSOA_patients_map = json.load(open(mappings_dir + 'GPs.json','r'))
    for GP in tqdm(LSOA_patients_map):
        for lsoa in LSOA_patients_map[GP]['Patient_registry_LSOA']:
            if lsoa not in LSOA_patient_pop:
                LSOA_patient_pop[lsoa] = LSOA_patients_map[GP]['Patient_registry_LSOA'][lsoa]
            else:
                LSOA_patient_pop[lsoa] += LSOA_patients_map[GP]['Patient_registry_LSOA'][lsoa]
    return LSOA_patient_pop


def calculateTemporalMetrics_LSOA(all_presc , old = True):
    LSOA_dosage = {}
    LSOA_costs = {}
    LSOA_patient_count = {}
    fail = 0.0
    LSOA_map = {}
    [LSOA_dist_old , LSOA_dist_2021] = loadLSOA_mappings()

    #At this time we are using the same map for all files. Ideally, every year needs to have a different map.
    if old:
        dosageField = '8'
        costField = '19'
        practiceField = '2'
        LSOA_map = LSOA_dist_2021
    else:
        dosageField = 'TOTAL_QUANTITY'
        costField = '19'
        practiceField = 'PRACTICE_CODE'
        LSOA_map = LSOA_dist_2021

    for name, group in all_presc.groupby(practiceField):
        total_dosage = np.sum(group[dosageField])
        total_cost = np.sum(group[costField])
        if name in LSOA_map:        
            for k in LSOA_map[name]:
                if k not in LSOA_dosage:
                    LSOA_dosage[k] = 0.0
                    LSOA_costs[k] = 0.0
                LSOA_dosage[k]+= float(total_dosage)*float(LSOA_map[name][k])
                LSOA_costs[k]+= float(total_cost)*float(LSOA_map[name][k])
    
    return  LSOA_dosage , LSOA_costs


def writeResultFiles(monthly_borough_dosage_new , monthly_borough_costs_new , diseases):
    LSOA_patient_pop = prepare_lsoa_GP_population()
    for disease in tqdm(diseases):
        disease_dict = {'YYYYMM':[] , 'LSOA_CODE' : [] , 'Total_prescriptions' : [] ,'Dosage_ratio' :[] , 'Patient_count' : []}
        for yyyymm in monthly_borough_dosage_new:
            for LSOA_CODE in monthly_borough_dosage_new[yyyymm][disease]:
                if LSOA_CODE[0] == 'E':
                    disease_dict['YYYYMM'].append(yyyymm)
                    disease_dict['LSOA_CODE'].append(LSOA_CODE)
                    disease_dict['Total_prescriptions'].append(monthly_borough_dosage_new[yyyymm][disease][LSOA_CODE])
                    disease_dict['Dosage_ratio'].append(monthly_borough_costs_new[yyyymm][disease][LSOA_CODE])
                    disease_dict['Patient_count'].append(LSOA_patient_pop[LSOA_CODE])
        disease_df = pd.DataFrame.from_dict(disease_dict)
        filename = output_dir + disease+'_V2.csv.gz'
        disease_df.to_csv(filename,index=False,compression='gzip')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-c' , '--conditions' ,nargs="+" ,help="list of a conditions you need prescription prevalence for. Can be a list of single item ['item'] ")
    parser.add_argument('--isCat' , default=False, help=" Boolean flag to notify if the list is a set of conditions or categories. Default is False, implying conditions ")
    parser.add_argument('-s', "--start" , help="start year and month, format YYYYMM")
    parser.add_argument('-e', "--end" , help="end year and month, format YYYYMM")
    parser.add_argument('-odir', "--output_dir" , help="Directory for output files, default ../data_prep/")
    parser.add_argument('-idir', "--input_dir" , help="Directory for input nhs files, default: ../../BL_Work/openPrescribe/serialized/. Make sure that you have NHS files downloaded here. There is no sanity check here.")

    args = parser.parse_args()

    conditions = args.conditions
    start = str(args.start)
    end = str(args.end)

    idir = args.input_dir
    odir = args.output_dir

    if idir:
        print('overriding input dir')
        input_dir = idir
    if odir:
        print('overriding input dir')
        output_dir = odir


    files = glob.glob(input_dir + '/*.gz')
    files.sort()
    start_i = 0
    end_i = len(files)
    
    if start:
        for i in range(len(files)):
            year = files[i].split('/')[-1].split('.')[0].strip()
            if year == start:
                start_i = i
    if end:
        for i in range(len(files)):
            year = files[i].split('/')[-1].split('.')[0].strip()
            if year == end:
                end_i = i

    files_sub = files[start_i:end_i+1]

    print (f'Running for {conditions} between months {start} and {end}. The following files were selected \n\n\n')
    for f in files_sub:
        print(f)

    

    DiseaseDrugs, drugMap = DrugMatching(conditions , args.isCat)
    print(drugMap)

    monthly_borough_dosage_new = {}
    monthly_borough_costs_new = {}


    #Compute prevalence over selected files

    for f in tqdm(files_sub):
        month = f.split('/')[-1].split('.')[0]
        print("Working with month  " + month)
        if int(month) > 201906:
            old = False
        else:
            old = True
        
        monthly_borough_dosage_new[month] = {}
        monthly_borough_costs_new[month] = {}
        pdp = pd.read_csv(f,compression='gzip')
        for disease in tqdm(drugMap):
            print("Working with disease  " + disease)
            monthly_borough_dosage_new[month][disease] = {}
            monthly_borough_costs_new[month][disease] = {}
            drugs = drugMap[disease]
            opioids = pdp.loc[pdp['16'].isin(drugs)] #Original opioids

            monthly_borough_dosage_new[month][disease] , monthly_borough_costs_new[month][disease] = calculateTemporalMetrics_LSOA(opioids, old)


    print("Done computing LSOA level prescription prevalences, writing files")


    writeResultFiles(monthly_borough_dosage_new , monthly_borough_costs_new , conditions)

    print("Finished processing !!!!! ")
