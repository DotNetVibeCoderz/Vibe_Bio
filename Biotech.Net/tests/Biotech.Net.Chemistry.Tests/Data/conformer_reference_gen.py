# Generates Data/conformer_reference.json with RDKit (2026.09.1): explicit-H graphs, RDKit hybridisation,
# UFF parameters / energies / gradients at an ETKDGv3 geometry, the UFF-minimised energy from that geometry,
# and an ETKDGv3 conformer ensemble (heavy atoms) per molecule.
# Run: python conformer_reference_gen.py   (needs rdkit)
import json
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from rdkit.Chem import rdForceFieldHelpers as FFH
from rdkit.Chem import rdMolAlign
RDLogger.DisableLog('rdApp.*')

MOLS = {
    'aspirin': 'CC(=O)Oc1ccccc1C(=O)O',
    'caffeine': 'Cn1cnc2c1c(=O)n(C)c(=O)n2C',
    'ibuprofen': 'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O',
    'paracetamol': 'CC(=O)Nc1ccc(O)cc1',
    'imatinib': 'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1',
    'atorvastatin': 'CC(C)c1c(C(=O)Nc2ccccc2)c(-c2ccccc2)c(-c2ccc(F)cc2)n1CC[C@@H](O)C[C@@H](O)CC(=O)O',
    'diazepam': 'CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21',
    'morphine': 'CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5',
    'celecoxib': 'Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(N)(=O)=O)cc2)cc1',
    'sildenafil': 'CCCc1nn(C)c2c(=O)[nH]c(-c3cc(S(=O)(=O)N4CCN(C)CC4)ccc3OCC)nc12',
    'methotrexate': 'CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc1',
    'penicillinG': 'CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O',
    'glucose': 'OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O',
    'biotin': 'OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12',
    'benzamidine': 'NC(=N)c1ccccc1',
    'warfarin': 'CC(=O)C[C@@H](c1ccccc1)c1c(O)c2ccccc2oc1=O',
    'nicotine': 'CN1CCC[C@H]1c1cccnc1',
    'lidocaine': 'CCN(CC)CC(=O)Nc1c(C)cccc1C',
    'propranolol': 'CC(C)NC[C@@H](O)COc1cccc2ccccc12',
    'oseltamivir': 'CCOC(=O)C1=C[C@@H](OC(CC)CC)[C@H](NC(C)=O)[C@@H](N)C1',
    'ciprofloxacin': 'OC(=O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O',
    'captopril': 'C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O',
    'nitrobenzene': 'O=[N+]([O-])c1ccccc1',
    'acetate_ammonium': 'CC(=O)[O-].C[NH3+]',
    'phenylacetylene_nitrile': 'C#Cc1ccc(C#N)cc1',
    'dimethylsulfone_phosphate': 'CS(C)(=O)=O.COP(=O)(O)O',
    'thiophene_ester': 'COC(=O)c1cccs1',
    'cyclopropyl_cyclobutene': 'C1CC1C1=CCC1',
    'tamoxifen_like_EZ': r'CC/C(=C(\c1ccccc1)c1ccc(OCCN(C)C)cc1)c1ccccc1',
    'esomeprazole': 'COc1ccc2[nH]c([S@@](=O)Cc3ncc(C)c(OC)c3C)nc2c1',
}

def coords(m, cid=-1, heavy=False, nd=4):
    conf = m.GetConformer(cid)
    out = []
    for a in m.GetAtoms():
        if heavy and a.GetAtomicNum() == 1:
            continue
        p = conf.GetAtomPosition(a.GetIdx())
        out.append([round(p.x, nd), round(p.y, nd), round(p.z, nd)])
    return out

def paths(m):
    angles, torsions, inversions = [], [], []
    for j in m.GetAtoms():
        nb = [n.GetIdx() for n in j.GetNeighbors()]
        for a in range(len(nb)):
            for b in range(a + 1, len(nb)):
                angles.append([nb[a], j.GetIdx(), nb[b]])
        if len(nb) == 3 and j.GetAtomicNum() in (6, 7, 8, 15):
            inversions.append([nb[0], j.GetIdx(), nb[1], nb[2]])
    for bd in m.GetBonds():
        b, c = bd.GetBeginAtomIdx(), bd.GetEndAtomIdx()
        for a in m.GetAtomWithIdx(b).GetNeighbors():
            if a.GetIdx() == c: continue
            for d in m.GetAtomWithIdx(c).GetNeighbors():
                if d.GetIdx() in (b, a.GetIdx()): continue
                torsions.append([a.GetIdx(), b, c, d.GetIdx()])
    return angles, torsions, inversions

out = {}
for name, smi in MOLS.items():
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    r = {'smiles': smi}
    r['atoms'] = [a.GetAtomicNum() for a in m.GetAtoms()]
    r['bonds'] = [[b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()] for b in m.GetBonds()]
    r['hybridization'] = [str(a.GetHybridization()) for a in m.GetAtoms()]
    r['conjugated'] = [b.GetIsConjugated() for b in m.GetBonds()]
    ps = AllChem.ETKDGv3(); ps.randomSeed = 42
    AllChem.EmbedMolecule(m, ps)
    r['xyz'] = coords(m, nd=6)
    # UFF parameters (typing evidence)
    r['uff_bonds'] = [list(FFH.GetUFFBondStretchParams(m, b[0], b[1])) for b in r['bonds']]
    angles, torsions, inversions = paths(m)
    r['uff_angles'] = [[*a, *FFH.GetUFFAngleBendParams(m, *a)] for a in angles if FFH.GetUFFAngleBendParams(m, *a) is not None]
    tv = []
    for t in torsions:
        p = FFH.GetUFFTorsionParams(m, *t)
        if p is not None: tv.append([*t, p])
    r['uff_torsions'] = tv
    iv = []
    for t in inversions:
        p = FFH.GetUFFInversionParams(m, *t)
        if p is not None: iv.append([*t, p])
    r['uff_inversions'] = iv
    r['uff_vdw'] = [list(FFH.GetUFFVdWParams(m, i, i)) for i in range(m.GetNumAtoms())]
    ff = AllChem.UFFGetMoleculeForceField(m)
    r['uff_energy'] = ff.CalcEnergy()
    r['uff_grad'] = [round(g, 6) for g in ff.CalcGrad()]
    ff.Minimize(maxIts=20000, forceTol=1e-6, energyTol=1e-10)
    r['uff_min_energy'] = ff.CalcEnergy()
    r['uff_min_xyz'] = [[round(v, 6) for v in ff.Positions()[3 * i:3 * i + 3]] for i in range(m.GetNumAtoms())]
    # ETKDGv3 ensemble, heavy atoms only
    mh = Chem.AddHs(Chem.MolFromSmiles(smi))
    ps = AllChem.ETKDGv3(); ps.randomSeed = 42; ps.pruneRmsThresh = 0.5
    cids = AllChem.EmbedMultipleConfs(mh, numConfs=10, params=ps)
    r['ensemble'] = [coords(mh, c, heavy=True) for c in cids]
    # calibration: how well a second RDKit ensemble (other seed) covers this one (mean / max of the
    # symmetry-aware heavy-atom best-RMSD from each conformer above to the second ensemble)
    m2 = Chem.AddHs(Chem.MolFromSmiles(smi))
    ps2 = AllChem.ETKDGv3(); ps2.randomSeed = 7; ps2.pruneRmsThresh = 0.5
    c2 = AllChem.EmbedMultipleConfs(m2, numConfs=10, params=ps2)
    a = Chem.RemoveHs(mh); b = Chem.RemoveHs(m2)
    best = []
    for ci in cids:
        best.append(min(rdMolAlign.GetBestRMS(b, a, prbId=cj, refId=ci) for cj in c2))
    r['rdkit_self_coverage'] = [sum(best) / len(best), max(best)]
    # chirality: signed volume around each tagged centre using heavy-atom order (for checks independent of H)
    r['chiral_tags'] = [[a.GetIdx(), str(a.GetChiralTag())] for a in Chem.MolFromSmiles(smi).GetAtoms() if a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED]
    out[name] = r
json.dump({'rdkit': Chem.rdBase.rdkitVersion, 'mols': out}, open('conformer_reference.json', 'w'), separators=(',', ':'))
print('ok', len(out))
