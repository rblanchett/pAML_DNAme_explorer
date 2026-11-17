#!/usr/bin/env python3
"""Regenerate PUBLIC, aggregate-only CpG promoter explorer from Scripts 12C/13C.

Run: python scripts/build_promoter_explorer.py --annotation Table_12C_B_ALL_Directional_CpG_Annotations_GRCh38.csv --clinical Table_13C_02_Clinical_Directional_CpG_Integrated_Atlas.csv --output data/atlas.js --11d Table_11D_C_ALL_Selected_CpG_Direction_Associations.csv
Requires Python 3.9+ and pandas. NO patient-level values or identifiers are exported.
"""
import argparse, json, pathlib, math, collections
import pandas as pd

FLAGS = [
 'Primary_Restricted','Primary_Restricted_010','Neighborhood_Conserved','Neighborhood_Conserved_010',
 'Core_Target_06','Core_Source_06',
 'Prior_Bridge_Source_to_Target_07a','Prior_Bridge_Target_to_Source_07a',
 'Prior_11C_Target_Restricted','Prior_11C_Target_Conserved',
 'Prior_06b_A_Higher_015','Prior_06b_B_Higher_015',
 'Prior_06b_Shared_Higher_015','Prior_06b_Shared_Lower_015',
 'Prior_06b_Shared_Higher_Formal_015','Prior_06b_Shared_Lower_Formal_015',
 'Prior_06b_A_Defining_Concordant_Distinct_015','Prior_06b_B_Defining_Concordant_Distinct_015',
 'Prior_06b_Target_Higher_015','Prior_06b_Source_Higher_015',
 'Prior_06b_Target_Defining_Distinct_015','Prior_06b_Source_Defining_Distinct_015'
]
REQUIRED = ['CpG','Source','Target','Direction','Region_GRCh38','Promoter_Overlap_Genes','Mean_Beta_Primary','Mean_Beta_Incoming',
            'Mean_Beta_Source_Nonbridge','N_Primary','N_Incoming','N_Source_Nonbridge','BH_Incoming_vs_Nonbridge','BH_Primary_vs_Incoming',
            'BH_Primary_vs_Nonbridge','BH_TOST_Primary_Equiv_Incoming','W_Target_Loading_Percentile','W_Target_Component_Specificity',
            'W_Target_Column_Share','Delta_Incoming_Minus_Nonbridge','Chromosome','Start_hg38']+FLAGS

def num(v, decimals=7):
    try:
        f=float(v)
        return round(f,decimals) if math.isfinite(f) else None
    except (TypeError, ValueError):
        return None

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--annotation',required=True)
    ap.add_argument('--clinical',required=True)
    ap.add_argument('--11d',dest='d11',default=None,help='Optional unannotated 11D table for row/key integrity check')
    ap.add_argument('--output',required=True,help='Destination atlas.js')
    args=ap.parse_args()
    ann=pd.read_csv(args.annotation,low_memory=False)
    missing=set(REQUIRED)-set(ann.columns)
    if missing:raise ValueError('12C annotation missing: '+', '.join(sorted(missing)))
    clin=pd.read_csv(args.clinical,low_memory=False)
    for k in ('Population','Source','Target','N_Source','N_Bridge','N_NonBridge','Status','Has_11D_CpGs'):
        if k not in clin: raise ValueError('13C missing '+k)
    if args.d11:
        d11=pd.read_csv(args.d11,usecols=['CpG','Direction','Source','Target'],low_memory=False)
        if len(ann)!=len(d11) or ann[['CpG','Direction','Source','Target']].fillna('').astype(str).equals(d11[['CpG','Direction','Source','Target']].fillna('').astype(str)) is False:
            raise ValueError('11D and 12C association identity/order differ; do not merge blindly.')
        print('CHECK 11D: exactly matching all',len(d11),'associations')
    sel=ann.loc[ann.Region_GRCh38.eq('Promoter') & ann.Promoter_Overlap_Genes.notna()].copy()
    if sel.empty:raise ValueError('No promoter rows with real promoter-overlap gene annotations!')
    if not sel.Coordinate_Status.eq('GRCh38-mapped').all():raise ValueError('Unmapped promoter coordinates present')
    if sel[['CpG','Source','Target']].duplicated().any():raise ValueError('Duplicate CpG/Source/Target association')
    all_c=clin.loc[clin.Population.eq('All')].copy()
    if all_c[['Source','Target']].duplicated().any():raise ValueError('Duplicate 13C pooled pair')
    pair_lookup={(str(r.Source),str(r.Target)):r for r in all_c.itertuples(index=False)}
    dirs=sorted(set(zip(sel.Source.astype(str),sel.Target.astype(str))),key=lambda s:(int(s[0][3:]),int(s[1][3:])))
    directions=[];direction_idx={}
    for s,t in dirs:
        c=pair_lookup.get((s,t))
        if c is None or str(c.Status)!='eligible' or not bool(c.Has_11D_CpGs):
            raise ValueError(f'12C promoter direction {s}->{t} does not match eligible 13C All')
        direction_idx[(s,t)]=len(directions)
        directions.append([s,t,int(c.N_Source),int(c.N_Bridge),int(c.N_NonBridge)])
    probes=[];probe_idx={};probe_check={};gene_rows=collections.defaultdict(list)
    output_rows=[];groups_by_probe={}; conflict_check=0
    # Pool source bridge/nonbridge means using the CpG-specific numbers with observed beta.
    def update_mean(group,p,n,mean):
        nonlocal conflict_check
        if mean is None or not n or n<=0:return
        key=(group,p)
        if key not in groups_by_probe:
            groups_by_probe[key]=(mean,n)
        else:
            a,aa=groups_by_probe[key]
            if abs(a-mean)>1e-6 or aa!=n:
                conflict_check+=1
                if conflict_check<4:print('CONFLICT',group,probes[p][0],a,mean,aa,n)
    # source row iteration preserves upstream categories, associations, statistics.
    for tup in sel.itertuples(index=False,name=None):
        z=dict(zip(sel.columns,tup))
        cpg=str(z['CpG']); pp=probe_idx.get(cpg)
        if pp is None:
            pp=len(probes);probe_idx[cpg]=pp
            probes.append([cpg,str(z['Chromosome']),int(z['Start_hg38'])])
        else:
            if probes[pp][1:]!=[str(z['Chromosome']),int(z['Start_hg38'])]:
                raise ValueError(f'Inconsistent coordinates for {cpg}')
        src=str(z['Source']);tar=str(z['Target']);di=direction_idx[(src,tar)]
        mP=num(z['Mean_Beta_Primary'],12);mI=num(z['Mean_Beta_Incoming'],12);mN=num(z['Mean_Beta_Source_Nonbridge'],12)
        nP=int(z['N_Primary']);nI=int(z['N_Incoming']);nN=int(z['N_Source_Nonbridge'])
        mask=sum((1<<i) for i,k in enumerate(FLAGS) if bool(z[k]) and not pd.isna(z[k]))
        vals=[pp,di,num(mP),num(mI),num(mN),nP,nI,nN,
              num(z['BH_Incoming_vs_Nonbridge'],9),num(z['BH_Primary_vs_Incoming'],9),
              num(z['BH_Primary_vs_Nonbridge'],9),num(z['BH_TOST_Primary_Equiv_Incoming'],9),
              mask,num(z['W_Target_Loading_Percentile'],6),num(z['W_Target_Component_Specificity'],6),
              num(z['W_Target_Column_Share'],9),num(z['Delta_Incoming_Minus_Nonbridge'],7)]
        rid=len(output_rows);output_rows.append(vals)
        genes={s.strip() for s in str(z['Promoter_Overlap_Genes']).split(';') if s.strip()}
        for gene in genes:gene_rows[gene].append(rid)
        update_mean(tar,pp,nP,mP)
        if nI+nN>0 and mI is not None and mN is not None:
            update_mean(src,pp,nI+nN,(nI*mI+nN*mN)/(nI+nN))
    if conflict_check:raise ValueError(f'Inconsistent original group means ({conflict_check}); check data sources')
    for g,rr in gene_rows.items():gene_rows[g]=sorted(set(rr))
    groups=[];rawgroupcounts={}
    for z in all_c.itertuples(index=False):
        s=str(z.Source)
        if s in rawgroupcounts and rawgroupcounts[s]!=int(z.N_Source):raise ValueError(f'Group count inconsistent for {s}')
        rawgroupcounts[s]=int(z.N_Source)
    for i in range(1,24):
        group=f'nmf{i}'; groups.append([group,rawgroupcounts.get(group)])
    groupMeans=collections.defaultdict(list)
    for (g,pp),(val,n) in groups_by_probe.items():groupMeans[g].append([pp,round(val,7),n])
    for k in groupMeans:groupMeans[k].sort(key=lambda row:row[0])
    payload={'schema':'pAML-promoter-v1','meta':{'associations':len(output_rows),'uniqueCpGs':len(probes),
              'promoterGenes':len(gene_rows),'eligibleDirections':len(directions),'cohortPatients':854,
              'coordAssembly':'GRCh38','nodeMeaning':'Original primary NMF-group mean beta for available CpG sites (target-primary direct or reconstructed from source I+N)',
              'arrowMeaning':'Incoming bridge beta minus source nonbridge beta (default); statistical BH adjustment per direction and contrast'},
             'flags':FLAGS,'groups':groups,'probes':probes,'directions':directions,'rows':output_rows,
             'geneRows':dict(sorted(gene_rows.items(),key=lambda x:x[0].lower())),
             'groupMeans':dict(groupMeans)}
    output=pathlib.Path(args.output);output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('w',encoding='utf-8') as w:
        w.write('/* Generated aggregate-only promoter CpG associations. NO patient data. */\nwindow.PAML_PROMOTER_ATLAS=')
        json.dump(payload,w,separators=(',',':'),ensure_ascii=False,allow_nan=False)
        w.write(';\n')
    print('WROTE:',output,'bytes:',output.stat().st_size)
    print('STATS:',payload['meta'])
    print('AUDIT: reconstructed original group-CpG means:',len(groups_by_probe),'in',len(groupMeans),'groups, 0 inconsistencies')
    print('AUDIT: no patient IDs, individual beta values, USIs or private score tables exported')
    audit={'meta':payload['meta'],'original_group_mean_records':len(groups_by_probe),
           'groups_with_original_means':len(groupMeans),'11D_12C_row_match':bool(args.d11),
           'group_node_count':len(groups),'notes':['12C already carries 11D columns; 11D used for identity audit only.',
           'Only GRCh38 promoter overlaps are used. Semicolon-separated genes are expanded, preserving multi-gene overlaps.',
           'Means reconstructed by CpG-specific I/N observed counts for source groups; direct P means for target groups.',
           'No outcome data or private patient IDs exported.',
           'CpGs selected within the same cohort; not independent validation.',
           'Inference is pooled; 13C sex-specific clinical summaries are not sex-specific CpG statistics.']}
    auditp=output.with_name('build-audit.json');auditp.write_text(json.dumps(audit,indent=2),encoding='utf-8')
if __name__=='__main__':main()
