# Generates Data/rdkit_reference.json (RDKit 2026.09.1). Run: python rdkit_reference_gen.py
import json
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, rdMolDescriptors as D, Descriptors, QED, MACCSkeys, AllChem, Lipinski
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import rdFingerprintGenerator as G
RDLogger.DisableLog('rdApp.*')
smis = {
 'aspirin':'CC(=O)Oc1ccccc1C(=O)O','caffeine':'Cn1cnc2c1c(=O)n(C)c(=O)n2C','ibuprofen':'CC(C)Cc1ccc(cc1)C(C)C(=O)O',
 'paracetamol':'CC(=O)Nc1ccc(O)cc1','benzene':'c1ccccc1','pyridine':'c1ccncc1','pyrrole':'c1cc[nH]c1','furan':'c1ccoc1',
 'thiophene':'c1ccsc1','naphthalene':'c1ccc2ccccc2c1','indole':'c1ccc2[nH]ccc2c1','ethanol':'CCO','glycine':'NCC(=O)O',
 'nitrobenzene':'O=[N+]([O-])c1ccccc1','imatinib':'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1',
 'atorvastatin':'CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O',
 'diazepam':'CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21','morphine':'CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5',
 'cholesterol':'CC(C)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C',
 'celecoxib':'Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1','sildenafil':'CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12',
 'salicylic':'O=C(O)c1ccccc1O','pyridone':'O=c1cccc[nH]1','azulene':'c1ccc2cccc2cc1','acetonitrile':'CC#N','methylpyrrole':'Cn1cccc1',
 'chloroform':'ClC(Cl)Cl','glucose':'OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O','biphenyl':'c1ccc(-c2ccccc2)cc1','spiro':'C1CCC2(CC1)CCCC2',
 'norbornane':'C1CC2CCC1C2','adamantane':'C1C2CC3CC1CC(C2)C3','cyclohexanone':'O=C1CCCCC1','toluene':'Cc1ccccc1','phenol':'Oc1ccccc1',
 'anthracene':'c1ccc2cc3ccccc3cc2c1','pyrene':'c1cc2ccc3cccc4ccc(c1)c2c34','tetrazole':'c1nn[nH]n1','charged':'C[N+](C)(C)C.[Cl-]',
 'deuterated':'[2H]C([2H])([2H])O','c13':'[13CH4]','thiazole':'c1cscn1','quinoline':'c1ccc2ncccc2c1','purine':'c1ncc2[nH]cnc2n1',
 'methotrexate':'CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc1','penicillinG':'CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O',
}
out={}
for k,s in smis.items():
    m=Chem.MolFromSmiles(s)
    r={}
    r['smiles']=s
    r['logp']=Crippen.MolLogP(m); r['mr']=Crippen.MolMR(m); r['tpsa']=D.CalcTPSA(m)
    r['mw']=Descriptors.MolWt(m); r['exact']=Descriptors.ExactMolWt(m); r['formula']=D.CalcMolFormula(m)
    r['hbd']=Lipinski.NumHDonors(m); r['hba']=Lipinski.NumHAcceptors(m); r['rotb']=D.CalcNumRotatableBonds(m)
    r['rings']=D.CalcNumRings(m); r['arom']=D.CalcNumAromaticRings(m); r['qed']=QED.qed(m); r['fsp3']=D.CalcFractionCSP3(m)
    r['nhoh']=Lipinski.NHOHCount(m); r['no']=Lipinski.NOCount(m)
    r['aromatic_atoms']=[a.GetIdx() for a in m.GetAtoms() if a.GetIsAromatic()]
    r['morgan2']=list(AllChem.GetMorganFingerprintAsBitVect(m,2,2048).GetOnBits())
    r['morgan_counts']={str(k):v for k,v in AllChem.GetMorganFingerprint(m,2).GetNonzeroElements().items()}
    r['fcfp2']=list(AllChem.GetMorganFingerprintAsBitVect(m,2,2048,useFeatures=True).GetOnBits())
    r['maccs']=list(MACCSkeys.GenMACCSKeys(m).GetOnBits())
    r['ap']={str(k):v for k,v in D.GetAtomPairFingerprint(m).GetNonzeroElements().items()}
    r['tt']={str(k):v for k,v in D.GetTopologicalTorsionFingerprint(m).GetNonzeroElements().items()}
    r['murcko']=Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m))
    r['canon']=Chem.MolToSmiles(m)
    out[k]=r
pairs=[('aspirin','salicylic'),('aspirin','ibuprofen'),('caffeine','theophylline'),('benzene','toluene'),('morphine','imatinib')]
smis['theophylline']='Cn1c2nc[nH]c2c(=O)n(C)c1=O'
tani={}
for a,b in pairs:
    ma,mb=Chem.MolFromSmiles(smis[a]),Chem.MolFromSmiles(smis[b])
    from rdkit import DataStructs
    fa,fb=AllChem.GetMorganFingerprintAsBitVect(ma,2,2048),AllChem.GetMorganFingerprintAsBitVect(mb,2,2048)
    tani[a+'|'+b]={'ecfp4':DataStructs.TanimotoSimilarity(fa,fb),'maccs':DataStructs.TanimotoSimilarity(MACCSkeys.GenMACCSKeys(ma),MACCSkeys.GenMACCSKeys(mb)),'smiles':[smis[a],smis[b]]}
json.dump({'mols':out,'tanimoto':tani},open('rdkit_reference.json','w'),indent=0)
print('ok')
