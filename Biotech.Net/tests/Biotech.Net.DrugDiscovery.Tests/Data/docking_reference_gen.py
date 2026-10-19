"""Builds the redocking fixtures and the AutoDock Vina reference used by FlexibleDockingTests.

Inputs are downloaded from RCSB into a work directory; outputs are written next to this script:
  <id>_receptor.pdb   protein only, hydrogens added by PDBFixer (pH 7), cropped to residues with an atom
                      within 18 A of the ligand centre (a pocket-only receptor, used identically by Vina and Biotech.Net)
  <id>_ligand.sdf     crystal ligand with bond orders from the CCD SMILES (heavy atoms, crystal coordinates)
  vina_reference.json AutoDock Vina score_only terms for the crystal pose and a redocking run
                      (box 22.5 A centred on the ligand, exhaustiveness 8, seed 42)

Tools (versions used): RDKit 2026.09.1, PDBFixer/OpenMM 8.6.1, Meeko 0.8.0 (ligand PDBQT), AutoDock Vina 1.2.7.
Run: python docking_reference_gen.py --vina C:/path/to/vina.exe [--work tmpdir]
"""
import argparse
import json
import os
import re
import subprocess
import urllib.request

from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

COMPLEXES = {
    # id: (ligand residue name, chain, CCD SMILES)
    '3ptb': ('BEN', 'A', 'NC(=N)c1ccccc1'),
    '1stp': ('BTN', 'A', 'OC(=O)CCCC[C@@H]1SC[C@@H]2NC(=O)N[C@H]12'),
}
BOX = 22.5
CROP = 18.0


def fetch(pdb_id, work):
    path = os.path.join(work, pdb_id + '.pdb')
    if not os.path.exists(path):
        urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pdb_id.upper()}.pdb', path)
    return path


def ligand_block(raw, resn, chain):
    return [l for l in open(raw) if l.startswith('HETATM') and l[17:20] == resn and l[21] == chain]


def write_ligand(lines, smiles, out):
    pdb = Chem.MolFromPDBBlock(''.join(lines) + 'END\n', removeHs=True)
    tmpl = Chem.MolFromSmiles(smiles)
    mol = AllChem.AssignBondOrdersFromTemplate(tmpl, pdb)
    Chem.AssignStereochemistryFrom3D(mol)
    mol.SetProp('_Name', os.path.basename(out).split('_')[0] + ' crystal ligand')
    w = Chem.SDWriter(out)
    w.write(mol)
    w.close()
    return mol


def protonate_and_crop(raw, center, out, work):
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    fixer = PDBFixer(filename=raw)
    fixer.removeHeterogens(keepWater=False)
    fixer.findMissingResidues()
    fixer.missingResidues = {}
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    full = os.path.join(work, os.path.basename(out) + '.full.pdb')
    with open(full, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f, keepIds=True)
    atoms = [l for l in open(full) if l.startswith(('ATOM', 'HETATM'))]
    keep = set()
    for l in atoms:
        x, y, z = float(l[30:38]), float(l[38:46]), float(l[46:54])
        d2 = (x - center[0]) ** 2 + (y - center[1]) ** 2 + (z - center[2]) ** 2
        if l[76:78].strip() != 'H' and d2 < CROP ** 2:
            keep.add((l[21], l[22:27]))
    with open(out, 'w', newline='\n') as f:
        f.write('REMARK   1 PROTONATED WITH PDBFIXER (PH 7); RESIDUES WITHIN %.0f A OF THE LIGAND CENTRE\n' % CROP)
        serial = 1
        for l in atoms:
            if (l[21], l[22:27]) in keep:
                f.write(l[:6] + '%5d' % serial + l[11:].rstrip('\n') + '\n')
                serial += 1
        f.write('END\n')


def run(cmd):
    return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout


AROMATIC_C = {'PHE': {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'}, 'TYR': {'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'},
              'TRP': {'CG', 'CD1', 'CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'}, 'HIS': {'CG', 'CD2', 'CE1'}}


def receptor_pdbqt(rec_pdb, work):
    """AD4-typed receptor PDBQT (ADFR conventions): polar H kept as HD, non-polar H merged, aromatic ring
    carbons A, unprotonated histidine ring N as NA, O as OA, S as SA. Partial charges are zero (Vina ignores them)."""
    atoms = []
    for l in open(rec_pdb):
        if l.startswith(('ATOM', 'HETATM')):
            atoms.append((l, l[12:16].strip(), l[17:20], l[76:78].strip().upper(),
                          (float(l[30:38]), float(l[38:46]), float(l[46:54]))))
    def d2(a, b):
        return sum((a[k] - b[k]) ** 2 for k in range(3))
    heavy = [a for a in atoms if a[3] != 'H']
    hs = [a for a in atoms if a[3] == 'H']
    parent = {}
    for h in hs:
        best = min(heavy, key=lambda a: d2(a[4], h[4]))
        parent[id(h)] = best
    nh = {}
    for h in hs:
        p = parent[id(h)]
        nh[id(p)] = nh.get(id(p), 0) + 1
    out = []
    for a in atoms:
        line, name, res, el, xyz = a
        if el == 'H':
            if parent[id(a)][3] not in ('N', 'O'):
                continue
            t = 'HD'
        elif el == 'C':
            t = 'A' if name in AROMATIC_C.get(res, ()) else 'C'
        elif el == 'N':
            t = 'NA' if res == 'HIS' and name in ('ND1', 'NE2') and nh.get(id(a), 0) == 0 else 'N'
        elif el == 'O':
            t = 'OA'
        elif el == 'S':
            t = 'SA'
        else:
            t = el
        out.append(line[:54] + '  1.00  0.00    +0.000 %-2s' % t + chr(10))
    path = os.path.join(work, os.path.basename(rec_pdb)[:-4] + '.pdbqt')
    open(path, 'w').write(''.join(out))
    return path


def ligand_pdbqt(mol, work, name):
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    mh = Chem.AddHs(mol, addCoords=True)
    setup = MoleculePreparation().prepare(mh)[0]
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setup)
    path = os.path.join(work, name + '_ligand.pdbqt')
    open(path, 'w').write(pdbqt)
    return path, pdbqt


def parse_score(text):
    terms = {}
    for key, pat in [('affinity', r'Estimated Free Energy of Binding\s*:\s*(-?[\d.]+)'),
                     ('intermolecular', r'\(1\) Final Intermolecular Energy\s*:\s*(-?[\d.]+)'),
                     ('intramolecular', r'\(2\) Final Total Internal Energy\s*:\s*(-?[\d.]+)'),
                     ('torsional', r'\(3\) Torsional Free Energy\s*:\s*(-?[\d.]+)')]:
        m = re.search(pat, text)
        terms[key] = float(m.group(1)) if m else None
    return terms


def docked_poses(out_pdbqt, ref):
    from meeko import PDBQTMolecule, RDKitMolCreate
    scores = [float(l.split()[3]) for l in open(out_pdbqt) if l.startswith('REMARK VINA RESULT')]
    mol = Chem.RemoveHs(RDKitMolCreate.from_pdbqt_mol(PDBQTMolecule.from_file(out_pdbqt, skip_typing=True))[0])
    res = []
    for k, conf in enumerate(mol.GetConformers()):
        rms = rdMolAlign.CalcRMS(mol, ref, prbId=conf.GetId())  # symmetry-aware, in place (no superposition)
        res.append({'affinity': scores[k], 'rmsd_to_crystal': rms})
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--vina', required=True)
    ap.add_argument('--work', default='docking_work')
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True)
    here = os.path.dirname(os.path.abspath(__file__))
    ref = {'vina': run([a.vina, '--version']).strip(), 'box_size': BOX, 'crop_radius': CROP, 'complexes': {}}
    for pid, (resn, chain, smi) in COMPLEXES.items():
        raw = fetch(pid, a.work)
        lig = write_ligand(ligand_block(raw, resn, chain), smi, os.path.join(here, f'{pid}_ligand.sdf'))
        conf = lig.GetConformer()
        center = [sum(conf.GetAtomPosition(i)[k] for i in range(lig.GetNumAtoms())) / lig.GetNumAtoms() for k in range(3)]
        rec = os.path.join(here, f'{pid}_receptor.pdb')
        protonate_and_crop(raw, center, rec, a.work)
        rq = receptor_pdbqt(rec, a.work)
        lq, lq_text = ligand_pdbqt(lig, a.work, pid)
        score = parse_score(run([a.vina, '--receptor', rq, '--ligand', lq, '--score_only', '--autobox']))
        out = os.path.join(a.work, pid + '_out.pdbqt')
        box = ['--center_x', '%.3f' % center[0], '--center_y', '%.3f' % center[1], '--center_z', '%.3f' % center[2],
               '--size_x', str(BOX), '--size_y', str(BOX), '--size_z', str(BOX)]
        run([a.vina, '--receptor', rq, '--ligand', lq, *box, '--exhaustiveness', '8', '--seed', '42',
             '--num_modes', '9', '--out', out])
        poses = docked_poses(out, lig)
        torsdof = int(re.search(r'TORSDOF (\d+)', lq_text).group(1))
        ref['complexes'][pid] = {'ligand': resn, 'smiles': smi, 'center': center, 'crystal_score': score,
                                 'torsdof': torsdof, 'redock': poses}
        print(pid, score, poses[0])
    json.dump(ref, open(os.path.join(here, 'vina_reference.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
