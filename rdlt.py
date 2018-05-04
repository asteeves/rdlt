#! /usr/bin/env python

from __future__ import print_function
import sys, pickle, argparse

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def writeHeader(molname, loplsflag):
    print("""import "oplsaa.lt"    # <-- defines the standard "OPLSAA" force field""")
    if loplsflag:
        print("""import "loplsaa.lt"   # <-- custom parameters for long alkane chains taken from
                      #     Sui et al. J.Chem.Theory.Comp (2012), 8, 1459
                      #     To use the ordinary OPLSAA force field parameters,
                      #     (instead of the Sui et al. parameters), change the
                      #     atom types below from "@atom:81L","@atom:85LCH2" to
                      #     "@atom:81" and "@atom:85"  (defined in "oplsaa.lt")""")
    print("""{0} inherits OPLSAA {{""".format(molname))

def writeFooter(molname):
    print("""}} # {0}



# Note: You don't need to supply the partial partial charges of the atoms.
#       If you like, just fill the fourth column with zeros ("0.000").
#       Moltemplate and LAMMPS will automatically assign the charge later""".format(molname))

def writeAtoms(m):
    print("# atom-id  mol-id  atom-type charge      X         Y        Z")
    print("  write('Data Atoms') {")
    conf = m.GetConformer(0)
    for at in m.GetAtoms():
        point = conf.GetAtomPosition(at.GetIdx())
        print("\t{0}\t$mol\t{1}\t0\t{2:8.3f}\t{3:8.3f}\t{4:8.3f}".format(
                                    '$atom:'+at.GetSymbol()+str(at.GetIdx()),
                                    at.GetProp('AtomType'),
                                    point.x, point.y, point.z
                                    ))
    print("  }")

def writeBonds(m):
    bonds = m.GetBonds()
    print("  write('Data Bond List') {")
    for bond in bonds:
        b = bond.GetBeginAtom()
        e = bond.GetEndAtom()
        bname = b.GetSymbol()+str(b.GetIdx())
        ename = e.GetSymbol()+str(e.GetIdx())
        print("\t$bond:{0}\t$atom:{1}\t$atom:{2}".format(bname+ename,
                                                       bname, ename))
    print("  }")

def lt_to_molecule(ltfn):
    """Reads a moltemplate .lt file and returns an RDKit molecule for
    comparison purposes. Only works on .lt files with specific formatting.
    Doesn't have bond type perception, so doesn't generate useful smiles.
    ** Don't use this for anything. **
    """
    ltmol = Chem.RWMol()
    with open(ltfn,'r') as infile:
        for line in [line.strip() for line in infile if line.strip()]:
            # for this to work, the name of the atom must contain the
            #atomic symbol at the beginning. Otherwise would need to
            # look up based on type or mass.
            if line.split()[0][:5] == "$atom":
                label = line.split(':')[1].split()[0]
                #print(label)
                #construct a new atom by passing the atomic symbol
                #filter removes numbers
                newatom = Chem.Atom(''.join(filter(str.isalpha, label)))
                atomid = ltmol.AddAtom(newatom)
            elif line.split()[0][:5] == "$bond":
                #assumes bond - atom - atom style entries with atom id
                # in the field
                id1str = line.split()[1].split(':')[1]
                id1 = int(''.join(filter(str.isdigit, id1str)))
                id2str = line.split()[2].split(':')[1]
                id2 = int(''.join(filter(str.isdigit, id2str)))
                #this makes everything a single bond, so won't allow building
                # of a valid smiles from the incomplete graph
                ltmol.AddBond(id1,id2)
                #print(id1,id2)
    newmol = ltmol.GetMol()
    Chem.SanitizeMol(newmol)
    AllChem.EmbedMolecule(newmol,AllChem.ETKDG())
    print(Chem.MolToSmiles(newmol))


def read_cdict(cdictin):
    with open(cdictin, 'rb') as f:
        cdict = pickle.load(f)
    return cdict

def sum_of_charges(m, cdict):
    test_sum = 0
    print("# Given the current charge dictionary, the atoms will have the following charges:")
    for atom in m.GetAtoms():
        atype = atom.GetProp('AtomType')
        print("# Atom {0} is type {1} with charge {2}".format(atom.GetIdx(),atype,cdict[atype]))
        test_sum += cdict[atype]
    print("# The sum of the atomic charges is: {:.2f}".format(test_sum))
    if abs(test_sum) > 0.001:
        print("""
            # WARNING: The net charge appears to be non-zero! This may indicate
            incompatible atom types.
            """)

def generateFeatureDefn(fpath, fdefout, cdictout):
    """Write a feature definition file in RDKit style from the moltemplate
    conversion document. Only need to run this function if the conversion
    document has been changed.

    fpath -- file path of the moltemplate conversion doc
    fdefout -- file path to write the feature definition file
    cdictout -- file path to create a dictionary of atom types to charges
    """
    with open(fpath,'r') as infile, open(fdefout,'w') as outfile:
        feat_index = 0
        cdict = {}
        for line in [line.strip() for line in infile if line.strip()]:
            if line[0]!='*':
                el, atomname, typename, patt, lttype, chg, desc = [el.strip() for el in line.split("|")]
                # write lttype, SMARTS to feature definintion file
                # NOTE: using feature family to store the atom names is dangerous
                # because rdkit won't assign mutliple features in same family.
                # So I had to assign an index to make unique names [AHS]
                fdefn = \
"""
DefineFeature {0} {1}
Family {2}{3}
EndFeature""".format(lttype, patt, feat_index, atomname)
                # add new charge dictionary entry
                cdict[lttype]=float(chg)
                feat_index+=1
                outfile.write(fdefn)

    with open(cdictout,'wb') as f:
        pickle.dump(cdict,f, protocol=2)


def main():
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", "--smi",
                        help="Smiles string of the molecule to be atom-typed",
                        required=True)
    parser.add_argument("-n", "--name",
                        help="Name of the molecule to be typed",
                        default="LIG")
    parser.add_argument("-l", "--loplsflag",
                        help="Use the lopls atom atype definitions",
                        action="store_true")
    parser.add_argument("-f", "--fdef",
                        help="OPLS feature definition file path",
                        default="./opls_lt.fdefn")
    parser.add_argument("-x", "--lfdef",
                        help="LOPLS feature definition file path",
                        default="./lopls_lt.fdefn")
    parser.add_argument("-r", "--refresh",
                        help="""
                        Overwrite/make new feature defintion files in the current directory
                        given a file path containing a moltemplate conversion document.
                        Use caution, as this will overwrite files.
                        """)
    parser.add_argument("-c", "--charge",
                        help="Check the net charge of the molecule based on a charge dictionary",
                        action="store_true")
    args = parser.parse_args()

    #Build rdkit molecule from smiles and generate a conformer
    m = AllChem.AddHs(Chem.MolFromSmiles(args.smi))
    AllChem.EmbedMolecule(m,AllChem.ETKDG())

    # WARNING: This part is dumb. Will update the lopls definitions ONLY
    # if the lopls flag is used. If a path is passed with the refresh command
    #
    if args.refresh and args.loplsflag:
        generateFeatureDefn(args.refresh,'./lopls_lt.fdefn','./lopls_lt_dict.pkl')
    elif args.refresh:
        generateFeatureDefn(args.refresh,'./opls_lt.fdefn','./opls_lt_dict.pkl')

    #Build a feature factory from the defintion file and assign all features
    factory = Chem.ChemicalFeatures.BuildFeatureFactory(args.fdef)
    features = factory.GetFeaturesForMol(m)

    #Use the features to assign an atom type property
    [m.GetAtomWithIdx(f.GetAtomIds()[0]).SetProp('AtomType',f.GetType()) for f in features];

    #if lopls defitions are desired, redo the feature process
    # overwrite atomtypes
    if args.loplsflag:
        #print('loplsflag is {}'.format(loplsflag) )
        lfactory = Chem.ChemicalFeatures.BuildFeatureFactory(args.lfdef)
        lfeatures = lfactory.GetFeaturesForMol(m)
        #print(len(lfeatures))
        #for f in lfeatures:
        #    print(f.GetId(), f.GetFamily(), f.GetType(), f.GetAtomIds())
        [m.GetAtomWithIdx(f.GetAtomIds()[0]).SetProp('AtomType',f.GetType()) for f in lfeatures];
        #[print(at.GetProp('AtomType')) for at in m.GetAtoms()]

    #find untyped atoms
    #
    failure = False
    for at in m.GetAtoms():
        try:
            at.GetProp('AtomType')
        except KeyError:
            print("Atom {0} does not have an assigned atom type!".format(at.GetIdx()))
            failure = True
    #if any failed to type, quit
    if failure:
        sys.exit("""Refusing to write a .lt file without type assignments.
Check the SMARTS pattern that defines the expected atom type.""")


    #basic output
    writeHeader(args.name,args.loplsflag)
    writeAtoms(m)
    writeBonds(m)
    writeFooter(args.name)

    if args.charge:
        # Read charge dictionaries for testing
        opls_cdict = read_cdict('./opls_lt_dict.pkl')
        if args.loplsflag:
            lopls_cdict = read_cdict('./lopls_lt_dict.pkl')
            opls_cdict.update(lopls_cdict)

        sum_of_charges(m,opls_cdict)


if __name__=='__main__':
    main()
