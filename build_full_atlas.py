#!/usr/bin/env python3
"""Build the public, static, aggregate-only pAML NMF methylation atlas, v2.

Requires Python 3.9+, pandas; no patient-level data. Outputs self-contained
`data/atlas.js` and `data/build-audit.json`.

Run from the repository root:
python scripts/build_full_atlas.py \
  --annotation /path/to/Table_12C_B_ALL_Directional_CpG_Annotations_GRCh38.csv \
  --clinical /path/to/Table_13C_02_Clinical_Directional_CpG_Integrated_Atlas.csv \
  --11d /path/to/Table_11D_C_ALL_Selected_CpG_Direction_Associations.csv
"""
import argparse
import collections
import json
import math
from pathlib import Path
import pandas as pd

FLAGS = [
    'Primary_Restricted', 'Primary_Restricted_010',
    'Neighborhood_Conserved', 'Neighborhood_Conserved_010',
    'Core_Target_06', 'Core_Source_06',
    'Prior_Bridge_Source_to_Target_07a', 'Prior_Bridge_Target_to_Source_07a',
    'Prior_11C_Target_Restricted', 'Prior_11C_Target_Conserved',
    'Prior_06b_A_Higher_015', 'Prior_06b_B_Higher_015',
    'Prior_06b_Shared_Higher_015', 'Prior_06b_Shared_Lower_015',
    'Prior_06b_Shared_Higher_Formal_015', 'Prior_06b_Shared_Lower_Formal_015',
    'Prior_06b_A_Defining_Concordant_Distinct_015',
    'Prior_06b_B_Defining_Concordant_Distinct_015',
    'Prior_06b_Target_Higher_015', 'Prior_06b_Source_Higher_015',
    'Prior_06b_Target_Defining_Distinct_015',
    'Prior_06b_Source_Defining_Distinct_015',
]
R_FIELDS = [
    'CpG', 'Source', 'Target', 'Direction', 'Chromosome', 'Start_hg38',
    'End_hg38', 'Probe_Strand_hg38', 'Region_GRCh38',
    'CpG_Island_Context_GRCh38', 'Promoter_Overlap_Genes',
    'GeneBody_Overlap_Genes', 'Nearest_Gene', 'Distance_To_Nearest_Gene_bp',
    'Coordinate_Status', 'Mean_Beta_Primary', 'Mean_Beta_Incoming',
    'Mean_Beta_Source_Nonbridge', 'N_Primary', 'N_Incoming',
    'N_Source_Nonbridge', 'BH_Incoming_vs_Nonbridge',
    'BH_Primary_vs_Incoming', 'BH_Primary_vs_Nonbridge',
    'BH_TOST_Primary_Equiv_Incoming', 'P_Welch_Incoming_vs_Nonbridge',
    'P_Welch_Primary_vs_Incoming', 'P_Welch_Primary_vs_Nonbridge',
    'P_TOST_Primary_Equiv_Incoming', 'SD_Primary', 'SD_Incoming',
    'SD_Source_Nonbridge', 'Delta_Incoming_Minus_Nonbridge',
    'Delta_Primary_Minus_Incoming', 'Delta_Primary_Minus_Nonbridge',
    'W_Target_Loading_Percentile', 'W_Target_Component_Specificity',
    'W_Target_Column_Share', 'W_Target_Sensitivity_Weight', 'Conserved_Direction',
] + FLAGS
CLIN_FIELDS = ['Population', 'Source', 'Target', 'Status', 'N_Source', 'N_Bridge',
    'N_NonBridge', 'Has_11D_CpGs', 'Top_Feature', 'Top_Level',
    'Top_Feature_Type', 'Top_Effect_Bridge_Minus_NonBridge', 'Top_Effect_Unit',
    'Top_P_Primary', 'Top_FDR_Within_Population', 'Top_FDR_Overall',
    'Top_Primary_Test', 'N_Directional_CpG_Associations', 'N_Unique_CpGs',
    'N_PrimaryRestricted_CpGs', 'N_NeighborhoodConserved_CpGs',
    'N_Both_PrimaryAndConserved', 'N_PrimaryOnly', 'N_ConservedOnly',
    'N_Neither_Label', 'Has_11D_CpGs']

def n(v, digits=8):
    try:
        x = float(v)
        return round(x, digits) if math.isfinite(x) else None
    except (ValueError, TypeError):
        return None

def s(v):
    return None if pd.isna(v) or str(v).strip() == '' else str(v).strip()

def safe_int(v):
    x=n(v,0)
    return None if x is None else int(x)

def numeric_group(g):
    return int(g.removeprefix('nmf'))

def assert_columns(df,fields,name):
    missing = set(fields)-set(df.columns)
    if missing: raise ValueError(f'{name}: missing fields: {", ".join(sorted(missing))}')

def gene_names(v):
    return [x.strip() for x in str(v).split(';') if x.strip() and x.strip().lower()!='nan']

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--annotation',required=True)
    p.add_argument('--clinical',required=True)
    p.add_argument('--11d',dest='d11',help='Optional independent ordered CpG/direction key verification')
    p.add_argument('--output',default='data/atlas.js')
    p.add_argument('--public-min-cell',type=int,default=5,help='Suppress small subgroup counts; minimum 5')
    args=p.parse_args()
    if args.public_min_cell < 5: raise ValueError('Public subgroup minimum must be >=5')
    ann=pd.read_csv(args.annotation,low_memory=False)
    clin=pd.read_csv(args.clinical,low_memory=False)
    assert_columns(ann,R_FIELDS,'12C')
    assert_columns(clin,CLIN_FIELDS,'13C')
    if ann.duplicated(['CpG','Source','Target']).any(): raise ValueError('Duplicate CpG/source/target keys')
    if ann.Coordinate_Status.ne('GRCh38-mapped').any(): raise ValueError('Unexpected non-GRCh38-mapped entries')
    if ann['Region_GRCh38'].isna().any(): raise ValueError('Missing region class')
    verified=False
    if args.d11:
        d11=pd.read_csv(args.d11,usecols=['CpG','Source','Target','Direction'],low_memory=False)
        if not ann[['CpG','Source','Target','Direction']].astype(str).equals(d11[['CpG','Source','Target','Direction']].astype(str)):
            raise ValueError('11D/12C key/order mismatch; refuse unsafe merge')
        verified=True
    for c in ['CpG','Region_GRCh38','CpG_Island_Context_GRCh38','Chromosome','Start_hg38','End_hg38',
              'Probe_Strand_hg38','Promoter_Overlap_Genes','GeneBody_Overlap_Genes','Nearest_Gene','Distance_To_Nearest_Gene_bp']:
        if c!='CpG' and (ann.groupby('CpG',sort=False)[c].nunique(dropna=False)>1).any():
            raise ValueError(f'Inconsistent annotation across directions for CpG: {c}')
    if clin.duplicated(['Source','Target','Population']).any(): raise ValueError('13C duplicate population/direction')
    directions=sorted(set(zip(ann.Source,ann.Target)),key=lambda st:(numeric_group(st[0]),numeric_group(st[1])))
    all_clin=clin[clin.Population.eq('All')]
    look={(str(r.Source),str(r.Target)):r for r in all_clin.itertuples(index=False)}
    diridx={};dir_rows=[]
    for st in directions:
        r=look.get(st)
        if r is None or r.Status!='eligible': raise ValueError(f'12C direction missing eligible 13C pooled: {st}')
        diridx[st]=len(dir_rows)
        dir_rows.append([*st,int(r.N_Source),int(r.N_Bridge),int(r.N_NonBridge)])
    probes=[];probeidx={};gene_index=collections.defaultdict(list)
    rows=[];meanmap={};reconstruction_checks=0
    for rec in ann.itertuples(index=False,name=None):
        z=dict(zip(ann.columns,rec));cp=s(z['CpG'])
        pp=probeidx.get(cp)
        if pp is None:
            pp=len(probes);probeidx[cp]=pp
            probes.append([cp,s(z['Chromosome']),int(z['Start_hg38']),int(z['End_hg38']),s(z['Probe_Strand_hg38']),
                s(z['Region_GRCh38']),s(z['CpG_Island_Context_GRCh38']),s(z['Promoter_Overlap_Genes']),
                s(z['GeneBody_Overlap_Genes']),s(z['Nearest_Gene']),safe_int(z['Distance_To_Nearest_Gene_bp'])])
        src=str(z['Source']);tgt=str(z['Target']);di=diridx[(src,tgt)]
        means=[n(z[f],12) for f in ['Mean_Beta_Primary','Mean_Beta_Incoming','Mean_Beta_Source_Nonbridge']]
        obs=[safe_int(z[f]) for f in ['N_Primary','N_Incoming','N_Source_Nonbridge']]
        if any(v is None for v in obs):raise ValueError(f'NA observed counts for {cp} {src}->{tgt}')
        mask=sum(1<<i for i,f in enumerate(FLAGS) if not pd.isna(z[f]) and str(z[f]).lower() in ('true','1','t'))
        values=[pp,di,*means,*obs,
            *[n(z[k],10) for k in ['BH_Incoming_vs_Nonbridge','BH_Primary_vs_Incoming',
                                    'BH_Primary_vs_Nonbridge','BH_TOST_Primary_Equiv_Incoming']],
            mask,
            *[n(z[k],8) for k in ['W_Target_Loading_Percentile','W_Target_Component_Specificity','W_Target_Column_Share',
                                   'W_Target_Sensitivity_Weight','SD_Primary','SD_Incoming','SD_Source_Nonbridge',
                                   'P_Welch_Incoming_vs_Nonbridge','P_Welch_Primary_vs_Incoming','P_Welch_Primary_vs_Nonbridge',
                                   'P_TOST_Primary_Equiv_Incoming','Delta_Incoming_Minus_Nonbridge','Delta_Primary_Minus_Incoming',
                                   'Delta_Primary_Minus_Nonbridge']],
            s(z['Conserved_Direction'])]
        rid=len(rows);rows.append(values)
        region=probes[pp][5]
        if region=='Promoter': genes=gene_names(z['Promoter_Overlap_Genes'])
        elif region=='Gene body': genes=gene_names(z['GeneBody_Overlap_Genes'])
        elif region=='Intergenic': genes=gene_names(z['Nearest_Gene'])
        else: raise ValueError(f'Unexpected region {region}')
        for g in set(genes):gene_index[g].append(rid)
        def put(g,beta,count):
            nonlocal reconstruction_checks
            if beta is None or count is None or count<=0:return
            key=(g,pp)
            if key not in meanmap:meanmap[key]=(beta,count)
            else:
                prev,prevn=meanmap[key]
                if abs(prev-beta)>1e-6 or prevn!=count:
                    raise ValueError(f'Conflicting group beta value: {g}/{cp}: {prev,prevn} vs {beta,count}')
                reconstruction_checks+=1
        put(tgt,means[0],obs[0])
        if means[1] is not None and means[2] is not None and obs[1]+obs[2]>0:
            put(src,(means[1]*obs[1]+means[2]*obs[2])/(obs[1]+obs[2]),obs[1]+obs[2])
    groupcount={str(r.Source):int(r.N_Source) for r in all_clin.itertuples(index=False)}
    group_array=[[f'nmf{i}',(groupcount.get(f'nmf{i}') if groupcount.get(f'nmf{i}',0)>=args.public_min_cell else None)] for i in range(1,24)]
    group_means=collections.defaultdict(list)
    for (gr,pp),(beta,count) in meanmap.items():group_means[gr].append([pp,round(beta,8),count])
    for k in group_means:group_means[k].sort(key=lambda z:z[0])
    clinical=[];small_suppressed=0
    for r in clin.itertuples(index=False):
        status=str(r.Status)
        eligible=status=='eligible' and safe_int(r.N_Bridge) is not None and safe_int(r.N_Bridge)>=args.public_min_cell and safe_int(r.N_NonBridge)>=args.public_min_cell
        if not eligible:small_suppressed+=1
        clinical.append([str(r.Source),str(r.Target),str(r.Population),status,
            safe_int(r.N_Source) if eligible else None,
            safe_int(r.N_Bridge) if eligible else None,
            safe_int(r.N_NonBridge) if eligible else None,
            s(r.Top_Feature) if eligible else None,s(r.Top_Level) if eligible else None,
            s(r.Top_Feature_Type) if eligible else None,
            n(r.Top_Effect_Bridge_Minus_NonBridge,8) if eligible else None,
            s(r.Top_Effect_Unit) if eligible else None,
            n(r.Top_P_Primary,10) if eligible else None,
            n(r.Top_FDR_Within_Population,10) if eligible else None,
            n(r.Top_FDR_Overall,10) if eligible else None,
            s(r.Top_Primary_Test) if eligible else None,
            safe_int(r.N_Directional_CpG_Associations),safe_int(r.N_Unique_CpGs),
            safe_int(r.N_PrimaryRestricted_CpGs),safe_int(r.N_NeighborhoodConserved_CpGs),
            safe_int(r.N_Both_PrimaryAndConserved),safe_int(r.N_PrimaryOnly),
            safe_int(r.N_ConservedOnly),safe_int(r.N_Neither_Label),
            bool(r.Has_11D_CpGs)])
    gene_index={k:sorted(set(v)) for k,v in sorted(gene_index.items(),key=lambda kv:kv[0].lower())}
    payload={
        'schema':'pAML-methylation-v2',
        'meta':{'cohortPatients':854,'associations':len(rows),'uniqueCpGs':len(probes),
          'geneLabels':len(gene_index),'eligibleDirectionsWithCpGs':len(dir_rows),
          'allClinicalDirectionalRecords':len(clin),'clinicalEligibleByPopulation':clin[clin.Status.eq('eligible')].groupby('Population').size().to_dict(),
          'regionAssociationCounts':ann.Region_GRCh38.value_counts().to_dict(),
          'uniqueCpGsByRegion':ann.groupby('Region_GRCh38').CpG.nunique().to_dict(),
          'strictPromoterCpGs':ann[ann.Region_GRCh38.eq('Promoter')].CpG.nunique(),
          'strictPromoterDirections':ann[ann.Region_GRCh38.eq('Promoter')].Direction.nunique(),
          'source':'11D + 12C + 13C; observational discovery suite only',
          'assembly':'GRCh38','node':'Mean beta of available unique CpGs in original primary NMF group; not differential',
          'arrow':'Signed CpG-wise effect for selected contrast; not change over time',
          'sexMethylationAvailable':False,'publicSuppressionThreshold':args.public_min_cell},
        'flags':FLAGS,'groups':group_array,'probes':probes,'directions':dir_rows,
        'rows':rows,'geneRows':gene_index,'groupMeans':dict(group_means),'clinical':clinical
    }
    output=Path(args.output);output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('w',encoding='utf-8') as f:
        f.write('/* pAML aggregate discovery CpG atlas v2; no patient-level records. */\nwindow.PAML_METHYLATION_ATLAS=')
        json.dump(payload,f,separators=(',',':'),ensure_ascii=False,allow_nan=False)
        f.write(';\n')
    audit={
        'schema':payload['schema'],'input_11d_12c_exact_order_match':verified,
        'meta':payload['meta'],'group_probe_original_mean_records':len(meanmap),
        'group_mean_consistency_cross_checks':reconstruction_checks,
        'eligible_clinical_rows_with_unredacted_counts':len(clin)-small_suppressed,
        'clinical_rows_with_suppressed_counts':small_suppressed,
        'group_mean_coverage_groups':len(group_means),
        'flag_true_counts':{f:int(ann[f].sum()) for f in FLAGS},
        'notes':[
          '12C contains all 11D selected associations; 11D is independently verified by keys if supplied.',
          'Original primary group CpG means: directly reported target means or reconstructed from observed source I and N partition.',
          'Promoters are TSS -2000 to +500 bp as defined by 12C; gene-body overlap and intergenic nearest-gene contexts are distinct.',
          'Nearest-gene assignment is distance-based, not a demonstrated regulatory connection.',
          'All 22 overlapping CpG provenance flags are retained, including categories with zero hits.',
          'Clinical rows are separate summary records. Female and Male clinical records are NOT sex-stratified methylation tests.',
          'Ineligible clinical records have subgroup counts and clinical top-hit metadata withheld; minimum public cell size applies.',
          'Group-level means and gene summaries are descriptive; no gene-level or cross-direction q-values are calculated.',
          'No patient IDs, individual methylation profiles or private H signatures exported.'
        ]}
    output.with_name('build-audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf-8')
    print('Created',output,'size',output.stat().st_size)
    print('Rows:',len(rows),'CpGs:',len(probes),'Gene labels:',len(gene_index),'Directions:',len(dir_rows))
    print('Region association counts:',payload['meta']['regionAssociationCounts'])
    print('Group probe means:',len(meanmap),'consistent duplicate checks:',reconstruction_checks)
    print('Clinical summaries:',len(clin),'redacted for small/ineligible:',small_suppressed)

if __name__=='__main__':main()
